
############################################################################################
## get.largest.interval
############################################################################################

#' Get the largest interval for each gene, given multiple TSS and TTS annotations
#'
#' The input bed6 file can be derived from a gencode annotation file,
#' as described in the vignette
#' @param bed A bed6 frame with comprehensive gene annotations, defaults to NULL
#' @return A bed6 frame with the largest annotation for each gene
#' @export
#' @examples
#' # get intervals for furthest TSS and TTS +/- interval
#' bed.long = get.largest.interval(bed=dat0)
get.largest.interval = function(bed=NULL){

  # get most extreme start and end sites
  ends0 = bed[with(bed, order(gene, -end)), ]
  starts0 = bed[with(bed, order(gene, start)), ]
  ends = ends0[!duplicated(ends0$gene),]
  starts = starts0[!duplicated(starts0$gene),]
  bed.out = ends
  bed.out$start = starts$start

  return(bed.out)

} # get.largest.interval


############################################################################################
## read.count.transcript
############################################################################################

#' Evaluate reads and compute densities for each transcript of each gene.
#'
#' Input coordinates for evaluation of counts and densities.
#' Return the maximal values across all transcripts of a given gene.
#' For by="cnt" (default), we return the density associated with the max count for a specific gene.
#' For by="den", we return the density associated with the max count for a specific gene.
#' @param bed A bed6 frame
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param by Either "cnt" (default) or "den" for selecting transcripts based on either count or density
#' @return
#' A list with counts and desity.
#' counts = vector with end of gene counts for each gene.
#' density = vector with end of gene densities for each gene.
#' @export
#' @examples
#' # get read counts and densities for the max transcript of each annotated gene
#' transcript.reads = read.count.transcript(bed=gencode.transcript, bw.plus=bw.plus, bw.minus=bw.minus)
#' hist(log(transcript.reads$density), breaks=200, col="black",xlab="log read density",main="")
#' hist(log(transcript.reads$counts), breaks=200, col="black",xlab="log read count",main="")

read.count.transcript = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, by="cnt"){

  genes = unique(bed$gene)
  res = c()
  for(gene in genes){

    # transcript coords
    ind = which(bed$gene == gene)
    dat = bed[ind,]
    strand = dat$strand[1]

    # count reads and compute densities for all transcripts
    # sort upstream to downstream
    if(strand=="-"){
      dat.bed = dat[order(dat$end,decreasing=T),]
      counts = bed.region.bpQuery.bigWig(bw.minus, dat.bed)
    }
    if(strand=="+"){
      dat.bed = dat[order(dat$start),]
      counts = bed.region.bpQuery.bigWig(bw.plus, dat.bed)
    }
    dens = counts / (dat$end - dat$start)

    # tabulate the results, take the upstream transcript for ties
    if(by == "cnt"){ind.max = which(counts == max(counts))[1]}
    if(by == "den"){ind.max = which(dens == max(dens))[1]}
    out = data.frame(gene=gene, count=counts[ind.max], density=dens[ind.max])
    res = rbind(res, out)

  } # gene loop

  # return results
  counts = res$count
  density = res$density
  names(counts) = names(density) = genes
  return(list(counts=counts, density=density))

} # read.count.transcript


############################################################################################
## read.count.body
############################################################################################

#' Count gene body reads and compute density
#'
#' Input coordinates for evaluation of counts and densities.
#' This function works with input that includes one annotation per gene.
#' @param bed A bed6 frame
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param bp.start Optional number of bases to omit from the upstream gene boundary
#' @return
#' A list with counts and desity.
#' counts = vector with end of gene counts for each gene.
#' density = vector with end of gene densities for each gene.
#' @export
#' @examples
#' # get read counts and densities at the body of each annotated gene
#'
#' body.reads = read.count.body(bed=largest.interval.bed, bw.plus=bw.plus,
#' bw.minus=bw.minus, bp.start=bp.start)
#' hist(log(body.reads$density), breaks=200, col="black",xlab="log read density",main="")
#' hist(log(body.reads$counts), breaks=200, col="black",xlab="log read count",main="")
read.count.body = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, bp.start=0){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")

  # modify coords to only look at the body of the gene
  bed.plus$start = bed.plus$start + bp.start
  bed.minus$end = bed.minus$end - bp.start

  # get read counts
  counts.plus = bed.region.bpQuery.bigWig(bw.plus, bed.plus)
  counts.minus = bed.region.bpQuery.bigWig(bw.minus, bed.minus)
  names(counts.plus) = bed.plus$gene
  names(counts.minus) = bed.minus$gene

  # aggregate end count information, compute densities, return results
  counts = c(counts.plus, counts.minus)
  plus.denom = bed.plus$end-bed.plus$start
  minus.denom = bed.minus$end-bed.minus$start
  density.minus = counts.minus/minus.denom
  density.plus = counts.plus/plus.denom
  count.density = c(density.plus, density.minus)
  return(list(counts=counts, density=count.density))

} # read.count.body


############################################################################################
## read.count.end
############################################################################################

#' Select a regions containing the gene ends
#'
#' Select a fraction of the annotated gene end and consider an additional number of base pairs beyond
#' the gene end within which to count reads.
#' This function works with input that includes one annotation per gene.
#' @param bed A bed6 frame
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param fraction.end Fraction of the annotated gene end (0,1)
#' @param add.to.end Number of base pairs beyond the gene end within which to count reads
#' @return
#' A list with counts and desity.
#' counts = vector with end of gene counts for each gene.
#' density = vector with end of gene densities for each gene.
#' @export
#' @examples
#' # get read counts and densities at the end of each annotated gene
#' fraction.end = 0.1
#' add.to.end = 0
#' end.reads = read.count.end(bed=bed.long, bw.plus=bw.plus, bw.minus=bw.minus,
#'                           fraction.end=fraction.end, add.to.end=add.to.end)
#' hist(log(end.reads$density))
#' hist(log(end.reads$counts))
read.count.end = function(bed=NULL, bw.plus=NULL, bw.minus=NULL,
                          fraction.end=NULL, add.to.end=NULL){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")

  # modify plus strand to only look at the end of the gene
  end1 = bed.plus$end
  bed.plus$end = end1 + add.to.end
  bed.plus$start = end1 - round( fraction.end * (end1 - bed.plus$start) )

  # modify minus strand to only look at the end of the gene
  str1 = bed.minus$start
  bed.minus$start = str1 - add.to.end
  bed.minus$end = str1 + round( fraction.end * (bed.minus$end - str1) )

  # get read counts
  counts.plus = bed.region.bpQuery.bigWig(bw.plus, bed.plus)
  counts.minus = bed.region.bpQuery.bigWig(bw.minus, bed.minus)
  names(counts.plus) = bed.plus$gene
  names(counts.minus) = bed.minus$gene

  # aggregate end count information, compute densities, return results
  counts = c(counts.plus, counts.minus)
  plus.denom = bed.plus$end-bed.plus$start
  minus.denom = bed.minus$end-bed.minus$start
  density.minus = counts.minus/minus.denom
  density.plus = counts.plus/plus.denom
  count.density = c(density.plus, density.minus)
  return(list(counts=counts, density=count.density))

} # read.count.end


############################################################################################
## get.TSS
############################################################################################

#' Select an optimal transcription start site (TSS) for each gene
#'
#' We set a range around each annotated exon 1 within which to search for a read density.
#' For each gene, we select the upstream coordinate of the first exon with the greatest read density, 
#' in the specified region (bp.range),
#' out of all first exons annotated for the given gene.
#' Note that input and output genes must be matched (bed.in and bed.out).
#' @param bed.in A bed6 frame with first exon gene annotations
#' @param bed.out A bed6 frame for output (same genes as bed.in)
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param bp.range The range in bp for identifying the exon 1 with the greatest density (e.g., c(20,120))
#' @param cnt.thresh Read count number below which the TSS is not evaluated, 
#' by default the gene is removed from the analysis for max(counts) < cnt.thresh (low.read=FALSE).
#' Indicate low.read=TRUE to use the upstream-most TSS for max(counts) < cnt.thresh.
#' @param low.read Logical for whether to pick the upstream-most TSS if max(counts) < cnt.thresh (default FALSE)
#' @param by indicate whether to select the TSS based on max read count (by="cnt", default) or density (cnt="den")
#' @return
#' A list with a bed6 frame with inferred TSSs (bed) and a list of problematic genes (issues).
#' Problems documented include start > end and start < 0. For start > end,
#' the original input corrdinates are retained (bed.out).
#' For start < 0, start = 0 is applied.
#' @export
#' @examples
#' # select the TSS for each gene and incorporate these TSSs into the largest interval coordinates
#' bp.range = c(20,120)
#' cnt.thresh = 5
#' bed.out = largest.interval.expr.bed
#' bed.in = gencode.firstExon[gencode.firstExon$gene %in% bed.out$gene,]
#' TSS.gene = get.TSS(bed.in=bed.in, bed.out=bed.out,
#'                    bw.plus=bw.plus, bw.minus=bw.minus,
#'                    bp.range=bp.range, cnt.thresh=cnt.thresh)
get.TSS = function(bed.in=NULL, bed.out=NULL, bw.plus=NULL, bw.minus=NULL,
                   bp.range=NULL, cnt.thresh=NULL, low.read=FALSE, by="cnt"){

  # error handling. note that input and output genes must be matched
  if( length(unique(bed.in$gene)) != length(unique(bed.out$gene)) ){stop("input genes do not match output genes")}
  if( setequal(unique(bed.in$gene),unique(bed.out$gene)) == FALSE ){stop("input genes do not match output genes")}

  # loop through each unique gene
  # map reads a range for each annotated TSS
  bed0 = bed.out
  bed = bed.in
  genes = unique(bed$gene)
  tss = rep(0,length(genes))
  for(ii in 1:length(genes)){

    # get gene of interest and strand
    ind.gene = which(bed$gene == genes[ii])
    strand.gene = unique(bed$strand[ind.gene])[1]
    chr.gene = unique(bed$chr[ind.gene])[1]

    # get all sets of intervals around each annotated TSS
    # sort upstream to downstream
    if(strand.gene == "-"){
      starts = unique(bed$end[ind.gene])
      ranges = data.frame(chr=chr.gene, start=starts-bp.range[2],
                  end=starts-bp.range[1], stringsAsFactors=FALSE)
      ranges = ranges[order(ranges$end,decreasing=TRUE),]
    }
    if(strand.gene == "+"){
      starts = unique(bed$start[ind.gene])
      ranges = data.frame(chr=chr.gene, start=starts+bp.range[1],
                  end=starts+bp.range[2], stringsAsFactors=FALSE)
      ranges = ranges[order(ranges$start),]
    }

    # get counts and densities for each interval of interest
    if(strand.gene == "-"){ counts = bed.region.bpQuery.bigWig(bw.minus, ranges) }
    if(strand.gene == "+"){ counts = bed.region.bpQuery.bigWig(bw.plus, ranges) }
    dens = counts / (ranges$end - ranges$start)

    # select the best TSS, select the most upstream for ties
    if(by == "cnt"){ind.max = which(counts == max(counts))[1]}
    if(by == "den"){ind.max = which(dens == max(dens))[1]}
    
    if(counts[ind.max] < cnt.thresh){
      if(low.read==TRUE){
        if(strand.gene == "-"){ tss[ii] = ranges$end[1] + bp.range[1] }
        if(strand.gene == "+"){ tss[ii] = ranges$start[1] - bp.range[1] }
      }
      next
    }
    if(strand.gene == "-"){ tss[ii] = ranges$end[ind.max] + bp.range[1] }
    if(strand.gene == "+"){ tss[ii] = ranges$start[ind.max] - bp.range[1] }

  } # ii
  names(tss) = genes

  # note that genes with less counts than cnt.thresh in the max count region
  # and low.read=FALSE
  # remove such genes here
  in.tss.remove = which(tss == 0)
  tss.remove = names(tss)[in.tss.remove]
  bed.out = bed.out[!(bed.out$gene %in% tss.remove),]

  # apply TSS coordinates to the expression-filtered long gene annotations
  # output results
  issues = c()
  out = apply.TSS.coords(bed=bed.out, tss=tss)
  ind = which(out$end < out$start)
  if(length(ind) > 0){
    inds = sapply(out$gene[ind],function(x){which(bed0$gene==x)}) %>% unlist
    out[ind,] = bed0[inds,]
    issues = bed0$gene[inds]
  }
  ind = which(out$start < 0)
  if(length(ind) > 0){
    out$start[ind] = 0
    issues = out$gene[ind]
  }
  return(list(bed=out, issues=issues))

} # get.TSS


############################################################################################
## apply.TSS.coords
############################################################################################

#' Apply previously identified TSSs to an existing gene annotation
#'
#' This function will apply transcription start site (TSS) coordinates to a bed6 frame.
#' This function is embedded within get.TSS().
#' @param bed A bed6 frame with gene annotations
#' @param tss TSS estimates, from inside get.TSS().
#' Note that the TSS genes must be a subset of the genes in the bed6 data.
#' @return
#' A bed6 frame with estimated TSSs incorporated
#' @export
#' @examples
#' out = apply.TSS.coords(bed=bed.out, tss=tss)
apply.TSS.coords = function(bed=NULL, tss=NULL){

  # separate data by strand
  bed.plus = bed %>% filter(strand == "+")
  bed.minus = bed %>% filter(strand == "-")

  # get indices
  inds.plus = sapply(bed.plus$gene,function(x){which(names(tss)==x)}) %>% unlist
  inds.minus = sapply(bed.minus$gene,function(x){which(names(tss)==x)}) %>% unlist

  # set TSS and generate output
  bed.plus$start = tss[inds.plus]
  bed.minus$end = tss[inds.minus]
  out = rbind(bed.plus, bed.minus)
  return(out)

} # apply.TSS.coords


############################################################################################
## get.dups
############################################################################################

#' Get duplicates in which start sites match or end sites match
#'
#' This function identifies an extreme form of gene overlaps in which multiple identifiers
#' correspond to identical start or end coordinates.
#' @param bed A bed frame
#' @return
#' A bed frame with duplicate entries. Col 7 indicates ordered overlap loci.
#' @export
#' @examples
#' dups0 = get.dups(bed = TSS.gene)

get.dups = function(bed=NULL){

  # identical start coordinates
  dup.start = bed[(duplicated(bed[,1:2]) |
                     duplicated(bed[,1:2],
                                fromLast=TRUE)),]
  dup.start = dup.start[order(dup.start$chr,dup.start$start,dup.start$end),]

  # identical end coordinates
  dup.end = bed[(duplicated(bed[,c(1,3)]) |
                   duplicated(bed[,c(1,3)],
                              fromLast=TRUE)),]
  dup.end = dup.end[order(dup.end$chr,dup.end$start,dup.end$end),]

  # determine the count of duplicate cases and output
  out = unique(rbind(dup.start,dup.end)) %>% mutate(cases=0)
  out = out[order(out$chr, out$start, out$end),]
  if(nrow(out)==0){return(NULL)}
  track = 1
  for(ii in 1:nrow(out)){
    if(out$cases[ii] == 0){
      ind.start = which(out$start == out$start[ii])
      ind.end = which(out$end == out$end[ii])
      if(length(ind.start)>1){out$cases[ind.start] = track}
      if(length(ind.end)>1){out$cases[ind.end] = track}
      track = track + 1
    }
  }
  out = out[order(out$cases, out$chr, out$start, out$end),]
  return(out)

} # get.dups


############################################################################################
## remove.overlaps
############################################################################################

#' Filter genes with overlaps
#'
#' For a given case in which multiple genes overlap, this function evaluates read density
#' for all transcripts cooresponding to each gene.
#' The gene with the largest read count (by="cnt", default) or density (by="den") across all transcripts is retained in the input annotation.
#' The related functions get.dups() and gene.overlaps() can be used to identify overlap cases.
#' @param bed A bed6 frame with gene annotations
#' @param overlaps Frame with 6 columns in the bed6 format and a 7th column indicating integer overlap cases
#' @param by Indicate whether to select transcripts based on read counts (by="cnt") or densities (by="den",default)
#' @return
#' A bed6 frame with overlaps removed. That is, a modified version of the input bed.
#' @export
#' @examples
#' # run overlap analysis
#' overlap.data = gene.overlaps( bed = bed.long.filtered2.tss )
#' has.start.inside = overlap.data$has.start.inside
#' is.a.start.inside = overlap.data$is.a.start.inside
#' dim(has.start.inside)
#' dim(is.a.start.inside)
remove.overlaps = function(bed=NULL, overlaps=NULL, transcripts=NULL,
                           bw.plus=NULL, bw.minus=NULL, by="den"){

  # loop through each overlap case and identify genes to keep
  keep = rep(NA, length(unique(overlaps$cases)))
  for(ii in unique(overlaps$cases)){

    # get the genes for a given overlap locus
    dat = overlaps[overlaps$cases==ii,]
    if(nrow(dat)<2){stop("too few entries, <2 overlaps for a case")}
    genes = dat$gene

    # loop through genes, get the max density across all transcripts
    densit = rep(0,length(genes))
    count = rep(0,length(genes))
    names(densit) = names(count) = genes
    for(jj in genes){
      transcr = transcripts[transcripts$gene==jj,]
      strand = transcr$strand[1]
      if(strand=="+"){counts = bed.region.bpQuery.bigWig(bw.plus, transcr)}
      if(strand=="-"){counts = bed.region.bpQuery.bigWig(bw.minus, transcr)}
      den = counts / (transcr$end - transcr$start)
      densit[jj] = max(den)
      count[jj] = max(counts)
    } # gene loop

    if(by == "den"){
      keep[ii] = names(densit)[which(densit == max(densit))[1]]
    }
    if(by == "cnt"){
      keep[ii] = names(count)[which(count == max(count))[1]]
    }

  } # overlaps loop

  # remove overlapping genes and return the output
  remove.genes = overlaps$gene[!(overlaps$gene %in% keep)]
  out = bed[!(bed$gene %in% remove.genes),]
  return(out)

} # remove.overlaps



############################################################################################
## gene.overlaps
############################################################################################

#' Documentation of gene overlaps
#'
#' This function identifies gene overlaps.
#' @param bed A bed6 frame with processed gene annotations
#' @return
#' A list with has.start.inside, is.a.start.inside, and cases. The object has.start.inside is a bed6 frame that documents genes in which there are starts
#' of other genes. The object is.a.start.inside is a bed6 frame that documents genes that start inside other genes. The object
#' cases documents loci in which 2 or more genes overlap (col 7 of a bed format frame).
#' @export
#' @examples
#' # run overlap analysis
#' overlap.data = gene.overlaps( bed = TSS.gene.filtered1 )
#' has.start.inside = overlap.data$has.start.inside
#' is.a.start.inside = overlap.data$is.a.start.inside
#' head(has.start.inside)
#' head(is.a.start.inside)
#' head(overlap.data$cases)
#' length(unique(overlap.data$cases$gene))
#' 
gene.overlaps = function(bed=NULL){

  # separate by strand
  bed.plus = bed %>% filter(strand == "+")
  bed.minus = bed %>% filter(strand == "-")

  # data for documenting overlaps
  has.start.inside = c()
  is.a.start.inside = c()
  overlap.dat = rbind(bed.plus, bed.minus) %>% mutate(cases=0)
  case.ind = 1

  # filter plus by start sites to document end/start overlaps
  # loop through each chromosome
  chr.plus = unique(bed.plus$chr)
  for(ii in 1:length(chr.plus)){
    inds.chr = which(bed.plus$chr == chr.plus[ii])

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      start = bed.plus$start[inds.chr][jj]
      end = bed.plus$end[inds.chr][jj]
      ind.overlap = which(bed.plus$start[inds.chr][-jj]>=start &
                            bed.plus$start[inds.chr][-jj]<=end)
      if(length(ind.overlap) > 0){
        outer.gene = bed.plus$gene[inds.chr][jj]
        inner.genes = bed.plus$gene[inds.chr][-jj][ind.overlap]
        overlap.genes = paste0(inner.genes, collapse=",")

        # this gene has a start inside
        bed.orig = bed.plus[inds.chr[jj],]
        bed.orig$xy = paste0(overlap.genes, " starts inside this gene")
        has.start.inside = rbind(has.start.inside, bed.orig)

        # this gene starts inside another gene
        bed.orig = bed.plus[inds.chr[-jj][ind.overlap],]
        bed.orig$xy = paste0("this gene starts inside ", outer.gene)
        is.a.start.inside = rbind(is.a.start.inside, bed.orig)

        # document cases
        inds = sapply(c(outer.gene,inner.genes),function(x){which(overlap.dat$gene==x)}) %>% unlist
        if(max(unique(overlap.dat$cases[inds]))>0){
          uni = unique(overlap.dat$cases[inds])
          caseind = min(uni[which(uni>0)])
          overlap.dat$cases[inds] = caseind
          overlap.dat$cases[overlap.dat$cases == max(uni[which(uni>0)])] = caseind
        } else {
          overlap.dat$cases[inds] = case.ind
          case.ind = case.ind + 1
        }
      } # if overlaps...
    } ## jj

  } # ii, plus

  # filter minus by start sites to document end/start overlaps
  # loop through each chromosome
  chr.minus = unique(bed.minus$chr)
  for(ii in 1:length(chr.minus)){
    inds.chr = which(bed.minus$chr == chr.minus[ii])

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      start = bed.minus$end[inds.chr][jj]
      end = bed.minus$start[inds.chr][jj]
      ind.overlap = which(bed.minus$end[inds.chr][-jj]<=start &
                            bed.minus$end[inds.chr][-jj]>=end)
      if(length(ind.overlap) > 0){
        outer.gene = bed.minus$gene[inds.chr][jj]
        inner.genes = bed.minus$gene[inds.chr][-jj][ind.overlap]
        overlap.genes = paste0(inner.genes, collapse=",")

        # this gene has a start inside
        bed.orig = bed.minus[inds.chr[jj],]
        bed.orig$xy = paste0(overlap.genes, " starts inside this gene")
        has.start.inside = rbind(has.start.inside, bed.orig)

        # this gene starts inside another gene
        bed.orig = bed.minus[inds.chr[-jj][ind.overlap],]
        bed.orig$xy = paste0("this gene starts inside ", outer.gene)
        is.a.start.inside = rbind(is.a.start.inside, bed.orig)

        # document cases
        overlapers = c(outer.gene,inner.genes)
        inds = sapply(overlapers,function(x){which(overlap.dat$gene==x)}) %>% unlist
        if(max(unique(overlap.dat$cases[inds]))>0){
          uni = unique(overlap.dat$cases[inds])
          caseind = min(uni[which(uni>0)])
          overlap.dat$cases[inds] = caseind
          overlap.dat$cases[overlap.dat$cases == max(uni[which(uni>0)])] = caseind
        } else {
          overlap.dat$cases[inds] = case.ind
          case.ind = case.ind + 1
        }
      } # if overlaps...
    } ## jj

  } # ii, minus

  # format data and generate output
  # reset cases counts to interval 1 if necessary
  # e.g., length(unique(cases1$cases)) == max(unique(cases1$cases))
  cases0 = overlap.dat[which(overlap.dat$cases>0),]
  cases0 = cases0[order(cases0$cases),]
  cases1 = cases0
  for(ii in 1:length(unique(cases0$cases))){
    cases1$cases[cases0$cases == unique(cases0$cases)[ii]] = ii
  }
  return( list(has.start.inside=has.start.inside,
               is.a.start.inside=is.a.start.inside,
               cases=cases1) )

} # gene.overlaps



############################################################################################
## inside.starts
############################################################################################

#' Identify genes with multiple starts inside
#'
#' This function identifies genes within which >1 start is observed.
#' That is, relatively 'large' genes in which 2 or more genes originate.
#' Data generated from gene.overlaps() serve as an input to this function.
#' @param vec Input the xy column from the $is.a.start.inside output of gene.overlaps()
#' @return
#' A vector of genes within which multiple starts are observed.
#' @export
#' @examples
#' overlap.data = gene.overlaps( bed = TSS.gene.filtered1 )
#' is.a.start.inside = overlap.data$is.a.start.inside
#' mult.inside.starts = inside.starts(vec = is.a.start.inside$xy)
inside.starts = function(vec=NULL){
  mult.inside.starts = vec[duplicated(vec) | duplicated(vec, fromLast=TRUE)] %>% unique
  mult.inside.starts = sapply(mult.inside.starts,function(x){
    strsplit(x,"inside ")[[1]][2]}) %>% unname
  return(mult.inside.starts)
} # inside.starts



############################################################################################
## adjacent.gene.coords
############################################################################################

#' Function to identify coordinates for adjacent gene pairs
#'
#' Adjacent gene pairs can be identified based on manual analyses, 
#' despite overlaps in the largest interval coordinates.
#' We empirically define pause sites for these genes by binning the region spanned by both genes,
#' fitting smooth spline curves to the binned read counts, and identifying the two largest peaks separated
#' by a specified distance. We set a bin size and apply the constraint that the identified
#' pause peaks must be a given distance apart. 
#' For the spline fits, we set the
#' number of knots to the number of bins divided by specified setting.
#' We identify the TSSs as the exon 1 coordinates closest to the pause site peaks.
#' @param fix.genes A frame with upstream and downstream genes
#' @param exon1 A bed6 frame with first exon gene annotations
#' @param bed.long Long gene annotations, see get.largest.interval()
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param bp.bin The interval will be separated into adjacent bins of this size
#' @param shift.up Look upstream of the long gene start by this amount to search for a viable TSS
#' @param delta.tss Amount by which to shift the identied TSS upstream of the identified interval
#' @param knot.div The number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number
#' @param knot.thresh Minimum number of knots
#' @param diff.tss Minimum distance between TSSs
#' @param pause.bp Number of bp by which a pause site should be considered downstrean of a TSS
#' @param dist.from.start Distance between the TTS of the upstream gene and the TSS of the downstream gene
#' @param fname File name (.pdf) for output plots
#' @return
#' A frame with the adjusted coordinates for the input gene set, along with a plot. 
#' The titles for the plots indicate the upstream gene, the downstream gene, and the strand.
#' For genes on the plus strand, the upstream gene is on the left. 
#' For the minus strand, the upstream gene is on the right.
#' The red line denotes the spline fit and the vertical lines indicate the pause peaks.
#' The plot is intended for diagnostic purposes.
#' @export
#' @examples
#' # get the start sites for adjacent gene pairs
#' fix.genes = rbind(fix.genes.id, fix.genes.ov)
#' bp.bin = 5
#' knot.div = 40
#' shift.up = 100
#' delta.tss = 50
#' diff.tss = 1000
#' dist.from.start  = 50
#' adjacent.coords = adjacent.gene.coords(fix.genes=fix.genes, bed.long=TSS.gene,
#'                                        exon1=gencode.firstExon,
#'                                        bw.plus=bw.plus, bw.minus=bw.minus,
#'                                        knot.div=knot.div, bp.bin=bp.bin,
#'                                        shift.up=shift.up, delta.tss=delta.tss,
#'                                        dist.from.start=dist.from.start,
#'                                        diff.tss=diff.tss, fname="adjacentSplines.pdf")
#'
adjacent.gene.coords = function(fix.genes=NULL, bed.long=NULL, exon1=NULL, bw.plus=NULL,
                             bw.minus=NULL, bp.bin=NULL, shift.up=NULL,
                             dist.from.start = 50,
                             delta.tss=NULL, knot.div=4, knot.thresh=5,
                             diff.tss=1000, pause.bp=120, fname="adjacentSplines.pdf"){

  TSS = c()
  pdf(fname)
  par(mfrow=c(3,3))

  # loop through each gene set
  for(ii in 1:nrow(fix.genes)){

    # get counts over binned region spanning both genes
    first = fix.genes$upstream[ii]
    second = fix.genes$downstream[ii]
    ind.first = which(bed.long$gene==first)
    ind.second = which(bed.long$gene==second)
    bed.pair = bed.long[c(ind.first,ind.second),]
    range = c(min(bed.pair$start), max(bed.pair$end))
    strand = bed.long$strand[ind.first]
    main = paste0(first,", ",second,", ",strand)

    if(strand=="+"){
      range[1] = range[1] - shift.up
      starts = seq(range[1], range[2]-bp.bin, by=bp.bin)
      ends = seq(range[1]+bp.bin, range[2], by=bp.bin)
      len = min(length(starts),length(ends))
      segments = data.frame(bed.pair$chr[1], starts[1:len], ends[1:len],
                            stringsAsFactors=FALSE)
      names(segments) = c("chr","start","end")
      segments$end[nrow(segments)] = range[2]
      counts = bed.region.bpQuery.bigWig(bw.plus, segments)
    }
    if(strand=="-"){
      range[2] = range[2] + shift.up
      ends = seq(range[2], range[1]+bp.bin, by=-bp.bin)
      starts = seq(range[2]-bp.bin, range[1], by=-bp.bin)
      len = min(length(starts),length(ends))
      segments = data.frame(bed.pair$chr[1], starts[1:len], ends[1:len],
                            stringsAsFactors=F)
      names(segments) = c("chr","start","end")
      segments$start[nrow(segments)] = range[1]
      counts = bed.region.bpQuery.bigWig(bw.minus, segments)
    }

    # fit splines and identify peaks indicative of TSSs
    nknots = round( length(counts) / knot.div )
    if(nknots < knot.thresh){nknots = knot.thresh}
    if(strand=="+"){
      spl = smooth.spline(x=segments$start, y=counts, all.knots=FALSE, nknots=nknots)
      pks = findpeaks(spl$y, nups=0)
      pks = pks[order(-pks[,1]),]
      pk.coords = segments$start[pks[,2]]
      ind.tss2 = which(abs(diff(pk.coords)) > diff.tss)[1] + 1
      tss1 = pk.coords[1]
      tss2 = pk.coords[ind.tss2]
      tss = c(tss1,tss2) - delta.tss
      spl = smooth.spline(x=segments$start, y=counts, all.knots=FALSE, nknots=nknots)
      pred = predict(spl,segments$start)
      plot(segments$start,counts,main=main,xlab="coords (bp)")
      points(pred$x,pred$y,type="l",col="red", xlab="coords (bp)");
      abline(v=tss)
      new = cbind(c(first,second), c(min(tss),max(tss)))
    }
    if(strand=="-"){
      spl = smooth.spline(x=segments$end, y=counts, all.knots=FALSE, nknots=nknots)
      pks = findpeaks(spl$y, nups=0)
      pks = pks[order(-pks[,1]),]
      pk.coords = rev(segments$end)[pks[,2]]
      ind.tss2 = which(abs(diff(pk.coords)) > diff.tss)[1] + 1
      tss1 = pk.coords[1]
      tss2 = pk.coords[ind.tss2]
      tss = c(tss1,tss2) + delta.tss
      spl = smooth.spline(x=segments$end, y=counts, all.knots=FALSE, nknots=nknots)
      pred = predict(spl,segments$end)
      plot(segments$end,counts,main=main, xlab="coords (bp)")
      points(pred$x,pred$y,type="l",col="red", xlab="coords (bp)")
      abline(v=tss)
      new = cbind(c(second,first), c(min(tss),max(tss)))
    }
    new =as.data.frame(new, stringsAsFactors=FALSE); names(new) = c("gene","TSS")
    TSS = rbind(TSS,new)

  } # ii
  dev.off()
  TSS$TSS = as.numeric(TSS$TSS)

  # get the TSS based on the closest exon 1
  new.bed = bed.long[bed.long$gene %in% c(fix.genes$upstream,fix.genes$downstream),]
  for(ii in 1:nrow(fix.genes)){

    # specify gene, strand
    upgene = fix.genes$upstream[ii]
    dngene = fix.genes$downstream[ii]
    strand = bed.long$strand[bed.long$gene == upgene]
    pause.up = TSS$TSS[TSS$gene == upgene]
    pause.dn = TSS$TSS[TSS$gene == dngene]

    if(strand == "+"){
      # upstream ordered exon 1
      tss.up = exon1$start[exon1$gene == upgene] %>% unique
      tss.dn = exon1$start[exon1$gene == dngene] %>% unique
      tss.up = tss.up[order(tss.up)]
      tss.dn = tss.dn[order(tss.dn)]

      # downstream closest exon 1
      dp_tss = pause.dn - tss.dn
      ind.dn = which(abs(dp_tss) == min(abs(dp_tss)))[1]
      tss.dn = tss.dn[ind.dn]

      # upstream closest exon 1
      tss.up = tss.up[tss.up < tss.dn]
      dp_tss = pause.up - tss.up
      ind.up = which(abs(dp_tss) == min(abs(dp_tss)))[1]
      tss.up = tss.up[ind.up]

      # apply inferred TSS coords
      new.bed$start[new.bed$gene == upgene] = tss.up
      new.bed$start[new.bed$gene == dngene] = tss.dn
    }

    if(strand == "-"){
      # upstream ordered exon 1
      tss.up = exon1$end[exon1$gene == upgene] %>% unique
      tss.dn = exon1$end[exon1$gene == dngene] %>% unique
      tss.up = tss.up[order(tss.up,decreasing=TRUE)]
      tss.dn = tss.dn[order(tss.dn,decreasing=TRUE)]

      # downstream closest exon 1
      dp_tss = tss.dn - pause.dn
      ind.dn = which(abs(dp_tss) == min(abs(dp_tss)))[1]
      tss.dn = tss.dn[ind.dn]

      # upstream closest exon 1
      tss.up = tss.up[tss.up > tss.dn]
      dp_tss = tss.up - pause.up
      ind.up = which(abs(dp_tss) == min(abs(dp_tss)))[1]
      tss.up = tss.up[ind.up]

      # apply inferred TSS coords
      new.bed$end[new.bed$gene == upgene] = tss.up
      new.bed$end[new.bed$gene == dngene] = tss.dn
    }

  } # exon 1 fix.genes loop, tss

  # apply TTSs, correct the upstream TTS
  for(ii in 1:nrow(fix.genes)){

    # specify gene, strand
    upgene = fix.genes$upstream[ii]
    dngene = fix.genes$downstream[ii]
    strand = bed.long$strand[bed.long$gene == upgene]

    if(strand=="+"){
      tss.dn = new.bed$start[new.bed$gene == dngene]
      new.bed$end[new.bed$gene == upgene] = tss.dn - dist.from.start
    }

    if(strand=="-"){
      tss.dn = new.bed$end[new.bed$gene == dngene]
      new.bed$start[new.bed$gene == upgene] = tss.dn + dist.from.start
    }

  } # exon 1 fix.genes loop, tts

  return(new.bed)
} # adjacent.gene.coords


############################################################################################
## adjacent.coords.plot
############################################################################################


#' Plot adjacent gene pair coordinates
#'
#' Plot the coordinates, along with read data, for overlapping genes 
#' that were identified as an adjacent pair, 
#' with coordinates determined using adjacent.gene.coords().
#' @param adjacent.coords bed6 frame with coordinates for adjacent genes
#' @param pair Frame with one row, upstream and downstream genes
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param xper Fraction of the x-axis before the gene starts
#' @param yper Fraction of the y-axis with open space above the read data
#' @return
#' A plot with read data and gene coordinate annotations.
#' @export
#' @examples
#' # visualize coordinates for a specific pair
#' adjacent.coords.plot(adjacent.coords=adjacent.coords,
#'                   pair=fix.genes[4,],
#'                   bw.plus=bw.plus, bw.minus=bw.minus)
#'
adjacent.coords.plot = function(adjacent.coords=NULL, pair=NULL,
                                bw.plus=NULL, bw.minus=NULL, bp.bin=5,
                                xper=0.2, yper=0.3) {

  # coords for adjacent genes
  yper = yper + 1
  genes = c(pair$upstream, pair$downstream)
  strand = adjacent.coords$strand[adjacent.coords$gene %in% genes[1]]
  chr = adjacent.coords$chr[adjacent.coords$gene %in% genes[1]]
  coords0 = adjacent.coords[adjacent.coords$gene %in% genes,]
  min.str = min(coords0$start)
  max.end = max(coords0$end)
  totdist = max.end - min.str
  str = min.str - xper * totdist
  end = max.end + xper * totdist

  # map reads to the interval
  starts = seq(str, end-bp.bin, by=bp.bin)
  ends = seq(str+bp.bin, end, by=bp.bin) 
  len = min(length(starts),length(ends))
  segments = data.frame(chr, starts[1:len], ends[1:len],
                        stringsAsFactors=FALSE)
  names(segments) = c("chr","start","end")
  if(strand=="+"){
    counts = bed.region.bpQuery.bigWig(bw.plus, segments)
    bp.axis = segments$start
    main = paste0(pair$upstream, ", ", pair$downstream)
  }
  if(strand=="-"){
    counts = bed.region.bpQuery.bigWig(bw.minus, segments)
    bp.axis = segments$end
    main = paste0(pair$downstream, ", ", pair$upstream)
  }

  # plot read data
  max.reads = max(counts)
  ymax = yper * max.reads
  plot(bp.axis,counts,type="l",ylim=c(0,ymax), main=main,
       xlab=paste0("coordinates (bp), ",chr,", strand ",strand))

  # upstream coords
  dopen = ymax - max.reads
  coord.lev = max.reads + dopen/2
  upcoords = adjacent.coords[adjacent.coords$gene==pair$upstream,]
  x = seq(upcoords$start, upcoords$end, 1)
  y = rep(coord.lev, length(x))
  points(x,y,type="l",lwd=4)

  # downstream coords
  dopen = ymax - max.reads
  coord.lev = max.reads + dopen/3
  upcoords = adjacent.coords[adjacent.coords$gene==pair$downstream,]
  x = seq(upcoords$start, upcoords$end, 1)
  y = rep(coord.lev, length(x))
  points(x,y,type="l",lwd=4)

} # adjacent.coords.plot


############################################################################################
## get.end.intervals
############################################################################################

#' Get intervals at the gene ends for transcription termination site (TTS) estimation
#'
#' We define regions at the gene ends to examine read counts for TTS identification.
#' Note that transcription frequently extends beyond the poly-A site of a gene. To capture the end of
#' transcription, it is critical to examine regions beyond annotated gene boundries. We
#' evaluate evidence of transcriptional termination in regions extending from a 3’ subset of the gene
#' to a selected number of base pairs downstream of the most distal annotated gene end. We initiate
#' the search region with the end of the largest gene annotation (see get.largest.interval()).
#' We extend the search region up to a given distance past the annotated gene end. We also apply the constraint that a
#' TTS cannot be identified closer than a specified distance to the previously identified TSS of a downstream gene.
#' This analysis also incorporates the constraint that a gene end region identified
#' cannot cross the TSS of a downstream gene, thereby preventing gene overlaps on a given strand.
#' Thus, we clip the amount of bases added on to the gene end as necessary to avoid overlaps.
#' The output of this function is the input for get.TTS().
#' @param bed Processed bed6 frame with one entry per gene and identified TSSs
#' @param add.to.end Number of bases to add to the end of each gene to define the search region
#' @param fraction.end Fraction of the gene annotation to consider at the end of the gene
#' @param dist.from.start The maximal allowable distance between the end of one gene and the start of the next
#' @return
#' A bed6 frame for TTS evaluation, see get.TTS().
#' @export
#' @examples
#' # get intervals for TTS evaluation
#' add.to.end = 100000
#' fraction.end = 0.2
#' dist.from.start = 50
#' bed.for.tts.eval = get.end.intervals(bed=bed.long.filtered4.tss,
#'                                  add.to.end=add.to.end,
#'                                  fraction.end=fraction.end,
#'                                  dist.from.start=dist.from.start)
get.end.intervals = function(bed=NULL, add.to.end=NULL, fraction.end=NULL,
                             dist.from.start=NULL){

  # separate by strand
  bed.plus.orig = bed %>% filter(strand == "+")
  bed.minus.orig = bed %>% filter(strand == "-")
  bed.plus.new = bed.plus.orig
  bed.minus.new = bed.minus.orig

  # modify plus strand to only look at the end of the gene
  end1 = bed.plus.orig$end
  bed.plus.new$end = end1 + add.to.end
  bed.plus.new$start = end1 - round( fraction.end * (end1 - bed.plus.orig$start) )
  bed.plus.big = bed.plus.new

  # modify minus strand to only look at the end of the gene
  str1 = bed.minus.orig$start
  bed.minus.new$start = str1 - add.to.end
  bed.minus.new$end = str1 + round( fraction.end * (bed.minus.orig$end - str1) )
  bed.minus.big = bed.minus.new

  # filter plus by start sites to remove overlaps
  # loop through each chromosome
  chr.plus = unique(bed.plus.new$chr)
  for(ii in 1:length(chr.plus)){

    inds.chr = which(bed.plus.new$chr == chr.plus[ii])
    old.start = bed.plus.orig$start[inds.chr]

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      new.start = bed.plus.new$start[inds.chr][jj]
      new.end = bed.plus.new$end[inds.chr][jj]
      ind.overlap = which(old.start>=new.start & old.start<=new.end)
      if(length(ind.overlap) > 0){
        ind.closest = which(old.start[ind.overlap] == min(old.start[ind.overlap]))[1]
        old = bed.plus.big$end[inds.chr][jj]
        new = old.start[ind.overlap][ind.closest] - dist.from.start
        bed.plus.new$end[inds.chr][jj] = new
        bed.plus.new$xy[inds.chr][jj] = old - new # clip distance
      }
    } ## jj

    # make sure the start of the closest gene is not too close
    closest.start = sapply(bed.plus.new$end[inds.chr],function(end){
      deltas = bed.plus.orig$start[inds.chr] - end
      ind.close = which(deltas < dist.from.start & deltas > 0)
      if(length(ind.close) > 0){
        ind.cl = which(bed.plus.orig$start[inds.chr][ind.close] ==
                         min(bed.plus.orig$start[inds.chr][ind.close]))[1]
        return( bed.plus.orig$start[inds.chr][ind.close][ind.cl] )
      } else {
        return(end + dist.from.start + 100)
      }
    }) %>% unlist()

    # modify gene ends to be shifted from start sites from adjacent genes
    dists.plus = closest.start - bed.plus.new$end[inds.chr]
    ind.modify = which(dists.plus < dist.from.start)
    if(length(ind.modify)>0){
      old = bed.plus.big$end[inds.chr][ind.modify]
      new = closest.start[ind.modify] - dist.from.start
      bed.plus.new$end[inds.chr][ind.modify] = new
      bed.plus.new$xy[inds.chr][ind.modify] = old - new
    }

  } # ii, plus

  # filter minus by start sites to remove end/start overlaps
  # loop through each chromosome
  chr.minus = unique(bed.minus.new$chr)
  for(ii in 1:length(chr.minus)){

    inds.chr = which(bed.minus.new$chr == chr.minus[ii])
    old.start = bed.minus.orig$end[inds.chr]

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      new.start = bed.minus.new$end[inds.chr][jj]
      new.end = bed.minus.new$start[inds.chr][jj]
      ind.overlap = which(old.start>=new.end & old.start<=new.start)
      if(length(ind.overlap) > 0){
        ind.closest = which(old.start[ind.overlap] == max(old.start[ind.overlap]))[1]
        old = bed.minus.big$start[inds.chr][jj]
        new = old.start[ind.overlap][ind.closest] + dist.from.start
        bed.minus.new$start[inds.chr][jj] = new
        bed.minus.new$xy[inds.chr][jj] = new - old
      }
    } ## jj

    # make sure the start of the closest gene is not too close
    closest.start = sapply(bed.minus.new$start[inds.chr],function(end){
      deltas = end - bed.minus.orig$end[inds.chr]
      ind.close = which(deltas < dist.from.start & deltas > 0)
      if(length(ind.close) > 0){
        ind.cl = which(bed.minus.orig$end[inds.chr][ind.close] ==
                         max(bed.minus.orig$end[inds.chr][ind.close]))
        return( bed.minus.orig$end[inds.chr][ind.close][ind.cl] )
      } else {
        return(end - dist.from.start - 100)
      }
    }) %>% unlist()

    # modify gene ends to be shifted from start sites from adjacent genes
    dists.minus = bed.minus.new$start[inds.chr] - closest.start
    ind.modify = which(dists.minus < dist.from.start)
    if(length(ind.modify)>0){
      old = bed.minus.big$start[inds.chr][ind.modify]
      new = closest.start[ind.modify] + dist.from.start
      bed.minus.new$start[inds.chr][ind.modify] = new
      bed.minus.new$xy[inds.chr][ind.modify] = new - old
    }
  } # ii, minus

  # output
  # set xy column to zero for unclipped genes
  bed.out = rbind(bed.plus.new, bed.minus.new)
  bed.out$xy[bed.out$xy == "na"] = 0
  bed.out$xy = as.numeric(bed.out$xy)
  ind.prob = which(bed.out$start < 0)
  if(length(ind.prob) > 0){bed.out$start[ind.prob] = 0}
  return(bed.out)

} # get.end.intervals


############################################################################################
## get.TTS
############################################################################################

#' Defining gene ends
#'
#' Given regions within which to search for TTSs, we opperationally define the TTSs by binning the TTS search
#' regions, counting reads within the bins, fitting smooth spline curves to the bin counts, and
#' detecting points at which the curves peak and then decay towards zero. We applied the
#' constraint that there must be a specified number of bases in the TTS search region, otherwise the
#' TTS analysis is not applied. Similarly, if the number of knots identified is too low, then we set the
#' number of knots to a specified threshold. For this analysis, we set a peak search region at the beginning of
#' the TTS search region and identify the maximal peak from the spline fit. Then we identify the point
#' at which the spline fit decays to a threshold relative to the peak.
#' We reasoned that the peak search region should be largest, as a fraction of the TTS search region, for genes with
#' the greatest numbers of clipped bases. Such cases occur when the conventional
#' gene ends are proximal to identified TSSs, and we should include these entire regions for analysis of
#' the TTS. Similarly, we reasoned that for genes with substantially less clipped bases, and
#' correspondingly larger TTS search regions with greater potential for observing enhancers or divergent transcripts,
#' the peak search regions should be smaller fractions of the TTS search regions. We use
#' an exponential model to define the peak search regions.
#' The user should use the output of get.end.intervals() as an input to this function.
#' Note that gene ends will be set to zeros if there are no counts or if the interval is too small. 
#' The gene ends in tss will be the default if there are problems with estimating the TTS.
#' This function calls get.gene.end().
#' @param bed A bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param tss A bed6 gene coordinates containing estimated TSSs and largest interval gene ends 
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param bp.bin The interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh The TTS is defined as pk.thresh percent of max peak of the spline fit
#' @param knot.div the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
#' @param cnt.thresh Read count number below which the TTS is not evaluated
#' @param knot.thresh Minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end The maximal length of the search region (bp)
#' @param tau.dist Distance constant for the exponential defining the region for peak detection
#' @param frac.max Maximal fraction of gene end region for peal detection
#' @param frac.min Minimal fraction of gene end region for peal detection
#' @param sum.thresh Threshold for the min sum of read counts, below this default to the input TTS (default 20)
#' @return
#' A bed6 file with TTS estimates incorporated, original TTSs are present if the TTS could not be estimated.
#' The output is a list including the bed frame along with several metrics (see example below and vignette).
#' @export
#' @examples
#' # identify gene ends
#' add.to.end = max(bed.for.tts.eval$xy)
#' knot.div = 40
#' pk.thresh = 0.05
#' bp.bin = 50
#' knot.thresh = 5
#' cnt.thresh = 5
#' tau.dist = 10000
#' frac.max = 1
#' frac.min = 0.2
#' inferred.coords = get.TTS(bed=bed.for.tts.eval, tss=TSS.gene.filtered3,
#'                           bw.plus=bw.plus, bw.minus=bw.minus,
#'                           bp.bin=bp.bin, add.to.end=add.to.end,
#'                           knot.div=knot.div,
#'                           pk.thresh=pk.thresh, knot.thresh=knot.thresh,
#'                           cnt.thresh=cnt.thresh, tau.dist=tau.dist,
#'                           frac.max=frac.max, frac.min=frac.min)
#'
get.TTS = function(bed=NULL, tss=NULL, bw.plus=NULL, bw.minus=NULL, bp.bin=NULL,
                        pk.thresh=0.1, knot.div=40, knot.thresh=5,
                        cnt.thresh=5, add.to.end=50000, sum.thresh=20,
                        tau.dist=10000, frac.max=1, frac.min=0.2){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")

  # get TTS for each gene
  TTS.plus = get.gene.end(bed=bed.plus,
                          bw=bw.plus, bp.bin=bp.bin,
                          pk.thresh=pk.thresh,
                          knot.div=knot.div,
                          cnt.thresh=cnt.thresh,
                          knot.thresh=knot.thresh,
                          add.to.end=add.to.end,
                          tau.dist=tau.dist,
                          frac.max=frac.max,
                          frac.min=frac.min,
                          sum.thresh=sum.thresh)
  TTS.minus = get.gene.end(bed=bed.minus,
                            bw=bw.minus, bp.bin=bp.bin,
                            pk.thresh=pk.thresh,
                            knot.div=knot.div,
                            cnt.thresh=cnt.thresh,
                            knot.thresh=knot.thresh,
                            add.to.end=add.to.end,
                            tau.dist=tau.dist,
                            frac.max=frac.max,
                            frac.min=frac.min,
                            sum.thresh=sum.thresh)

  # apply identified TTS to the strand-separated search region frames 
  bed.plus.out = bed.plus
  bed.minus.out = bed.minus
  bed.plus.out$end = TTS.plus$tts
  bed.minus.out$start = TTS.minus$tts
  
  # address problematic cases in which the TTS cound not be identified 
  # such entries are zeros, see logic cases in get.gene.end()
  ind.pos = which(bed.plus.out$end == 0)
  if(length(ind.pos) > 0){
    genes = bed.plus.out$gene[ind.pos]
    inds = sapply(genes,function(x){which(tss$gene==x)})
    bed.plus.out$end[ind.pos] = tss$end[inds]
  }
  ind.neg = which(bed.minus.out$start == 0)
  if(length(ind.neg) > 0){
    genes = bed.minus.out$gene[ind.neg]
    inds = sapply(genes,function(x){which(tss$gene==x)})
    bed.minus.out$start[ind.neg] = tss$start[inds]
  }
  
  # aggregate data and address start<0 cases if necessary
  bed.out = rbind(bed.plus.out, bed.minus.out)
  ind.prob = which(bed.out$start < 0)
  if(length(ind.prob) > 0){bed.out$start[ind.prob] = 0}
  
  # apply TSSs to the frame with TTSs
  tss.plus = tss$start[tss$strand=="+"]
  tss.minus = tss$end[tss$strand=="-"]
  names(tss.plus) = tss$gene[tss$strand=="+"]
  names(tss.minus) = tss$gene[tss$strand=="-"]
  out = apply.TSS.coords(bed=bed.out, tss=c(tss.plus, tss.minus))
  
  # return the frame with TSSs and TTSs along with info on spline fitting
  return(list(bed=out,
              minus.lowcount=TTS.minus$lowcount,
              plus.lowcount=TTS.plus$lowcount,
              minus.knotmod=TTS.minus$knotmod,
              plus.knotmod=TTS.plus$knotmod))

} # get.TTS


############################################################################################
## get.gene.end
############################################################################################

#' Estimate the transcription termination sites (TTSs)
#'
#' This function estimates the TTSs. This function is called by get.TTS() 
#' and is not designed to be run on its own. TTSs that cannot be estimated are assigned zeros.
#' @param bed A bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param bw Plus or minus strand bigWig data
#' @param bp.bin The interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh The TTS is defined as pk.thresh percent of max peak of the spline fit (identified within the peak search region)
#' @param knot.div The number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number
#' @param cnt.thresh Read count number below which the TTS is not evaluated
#' @param knot.thresh Minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end The maximal length of the search region (bp)
#' @param tau.dist Distance constant for the exponential function defining the region for peak detection
#' @param frac.max Maximal fraction of gene end region for peak detection
#' @param frac.min Minimal fraction of gene end region for peak detection
#' @param sum.thresh Threshold for the min sum of read counts
#' @return
#' A vector of TTS coordinates
#' @export
#' @examples
#' See get.TTS()
get.gene.end = function(bed=NULL, bw=NULL, bp.bin=NULL, pk.thresh=NULL, knot.div=NULL,
                   cnt.thresh=NULL, knot.thresh=NULL, add.to.end=NULL,
                   tau.dist=NULL, frac.max=NULL, frac.min=NULL, sum.thresh=NULL){

  # check strand
  strand = unique(bed$strand)
  if(length(strand) > 1){stop("issue with strand separation in get.TTS()")}

  # loop through each gene and select a TTS coordinate
  tts = rep(0,nrow(bed)) # all zeros to start (default)
  lowcount = c()
  knotmod = c()
  for(ii in 1:nrow(bed)){

    # set region for peak analysis
    clip = bed$xy[ii]
    x = c(-add.to.end:0)
    y0 = exp(-x/tau.dist)
    y = ( (y0-min(y0))/(max(y0)-min(y0)) ) * (frac.max-frac.min) + frac.min
    frac.dist = y[add.to.end - clip + 1]

    # divide the gene into bins and obtain counts
    gene = bed$gene[ii]
    if( bed$end[ii] - bed$start[ii] < bp.bin ){next}
    if(strand=="+"){
      starts = seq(bed$start[ii], (bed$end[ii]-bp.bin), by=bp.bin)
      ends = seq((bed$start[ii]+bp.bin), bed$end[ii], by=bp.bin)
      lim = min( length(starts), length(ends) )
      bed.map = data.frame(bed$chr[ii], starts[1:lim], ends[1:lim],
                           stringsAsFactors=F)
      names(bed.map) = c("chr","start","end")
      bed.map$end[nrow(bed.map)] = bed$end[ii]
    }
    if(strand=="-"){
      ends = seq(bed$end[ii], (bed$start[ii]+bp.bin), by=-bp.bin)
      starts = seq((bed$end[ii]-bp.bin), bed$start[ii], by=-bp.bin)
      lim = min( length(starts), length(ends) )
      bed.map = data.frame(bed$chr[ii], starts[1:lim], ends[1:lim],
                           stringsAsFactors=F)
      names(bed.map) = c("chr","start","end")
      bed.map$start[nrow(bed.map)] = bed$start[ii]
    }
    counts = bed.region.bpQuery.bigWig(bw, bed.map)
    if(length(counts) < cnt.thresh){lowcount = c(lowcount,ii); next}
    if(sum(counts) < sum.thresh){next}

    # fit a spline to the count profile
    # get peaks from the spline fit and estimate the best TTS
    # TTS defined by pk.thresh percent of max peak of the spline fit
    nknots = round( length(bed.map$end) / knot.div )
    if(nknots < knot.thresh){knotmod = c(knotmod,ii); nknots = knot.thresh}
    if(strand=="+"){
      spl = smooth.spline(x=bed.map$end, y=counts, all.knots=FALSE, nknots=nknots)
      
      srch = c(min(bed.map$start), max(bed.map$end))
      dsrch = srch[2] - srch[1]
      fsrch = srch[1] + dsrch * frac.dist
      dd = spl$x - fsrch
      isearch = which(abs(dd) == min(abs(dd)))[1]
      
      pks = findpeaks(spl$y[1:isearch], nups=0)
      indpk = which( pks[,1] == max(pks[,1]) )[1]
      ind.tts = which(spl$y[pks[indpk,2]:length(spl$y)] < pk.thresh * max(pks[,1]))
      if(length(ind.tts) == 0){
        TTS = spl$x[length(spl$x)]
      } else {
        ind.tts = min( which(spl$y[pks[indpk,2]:length(spl$y)] <
                               pk.thresh * max(pks[,1])) ) + pks[indpk,2] - 1
        TTS = spl$x[ind.tts]
      }
    } # plus
    
    if(strand=="-"){
      spl = smooth.spline(x=-bed.map$start, y=counts, all.knots=FALSE, nknots=nknots)
      
      srch = c(min(bed.map$start), max(bed.map$end))
      dsrch = srch[2] - srch[1]
      fsrch = srch[2] - dsrch * frac.dist
      dd = spl$x + fsrch
      isearch = which(abs(dd) == min(abs(dd)))[1]
      
      pks = findpeaks(spl$y[1:isearch], nups=0)
      indpk = which( pks[,1] == max(pks[,1]) )[1]
      ind.tts = which(spl$y[pks[indpk,2]:length(spl$y)] < pk.thresh * max(pks[,1]))
      if(length(ind.tts)==0){
        TTS = -spl$x[length(spl$x)]
      } else {
        ind.tts = min( which(spl$y[pks[indpk,2]:length(spl$y)] <
                               pk.thresh * max(pks[,1])) ) + pks[indpk,2] - 1
        TTS = -spl$x[ind.tts]
      }
    } # minus

    tts[ii] = TTS
  } # ii
  return(list(tts=tts, lowcount=lowcount, knotmod=knotmod))
} # get.gene.end


############################################################################################
## tts.plot
############################################################################################

#' Plot inferred and largest interval coordinates
#'
#' This function plots the inferred and largest interval coordinates along with read data.
#' The plot annotates the peak search region and the total search region.
#' @param coords A bed6 frame with inferred coordinates
#' @param gene.end A bed6 frame with peak search regions
#' @param long.gene A bed6 frame with largest interval coordinates
#' @param gene A gene for plotting
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param add.to.end The maximal length of the search region (bp)
#' @param tau.dist Distance constant for the exponential defining the region for peak detection
#' @param frac.max Maximal fraction of gene end region for peak detection
#' @param frac.min Minimal fraction of gene end region for peak detection
#' @param xper Fraction of the x-axis before the gene starts
#' @param yper Fraction of the y-axis with open space above the read data
#' @return
#' A plot with reads, the interred annotation (black), the largest interval annotation (gray),
#' the TTS search region (solid, vertical), and the TTS peak region (dashed, vertical).
#' @export
#' @examples
#' # first run inferred.coords = get.TTS(...)
#' coords = inferred.coords$bed
#' gene.end = bed.for.tts.eval
#' long.gene = largest.interval.bed
#' tts.plot(coords=coords, gene.end=gene.end, long.gene=long.gene,
#'          gene = "Pparg", xper=0.2, yper=0.3,
#'          frac.min=frac.min, frac.max=frac.max,
#'          bw.plus=bw.plus, bw.minus=bw.minus, bp.bin=5,
#'          add.to.end=add.to.end, tau.dist=tau.dist)
#'
tts.plot = function(coords=NULL, gene.end=NULL, long.gene=NULL,
                    bw.plus=NULL, bw.minus=NULL, bp.bin=5, gene=NULL,
                    xper=0.2, yper=1.3, add.to.end=NULL, tau.dist=NULL,
                    frac.max=NULL, frac.min=NULL) {

  # coords for adjacent genes
  strand = coords$strand[coords$gene %in% gene]
  chr = coords$chr[coords$gene %in% gene]
  coord = coords[coords$gene == gene,]
  long = long.gene[long.gene$gene == gene,]

  # total search window, x-axis boundaries
  search = gene.end[gene.end$gene == gene,]
  min.str = min(c(coord$start, search$start, long$start))
  max.end = max(c(coord$end, search$end, long$end))
  totdist = max.end - min.str
  str = min.str - xper * totdist
  end = max.end + xper * totdist

  # map reads to the interval
  starts = seq(str, end-bp.bin, by=bp.bin)
  ends = seq(str+bp.bin, end, by=bp.bin)
  len = min(length(starts),length(ends))
  segments = data.frame(chr, starts[1:len], ends[1:len],
                        stringsAsFactors=FALSE)
  names(segments) = c("chr","start","end")
  if(strand=="+"){
    counts = bed.region.bpQuery.bigWig(bw.plus, segments)
    bp.axis = segments$start
    main = gene
  }
  if(strand=="-"){
    counts = bed.region.bpQuery.bigWig(bw.minus, segments)
    bp.axis = segments$end
    main = gene
  }

  # plot read data
  max.reads = max(counts)
  ymax = (1+yper) * max.reads
  plot(bp.axis,counts,type="l",ylim=c(0,ymax), main=main,
       xlab=paste0("coordinates (bp), ",chr,", strand ",strand))

  # inferred coords
  dopen = ymax - max.reads
  coord.lev = max.reads + dopen/2
  x = seq(coord$start, coord$end, 1)
  y = rep(coord.lev, length(x))
  points(x,y,type="l",lwd=4,col="darkgray")

  # largest interval coords
  c2 = long.gene[long.gene$gene == gene,]
  coord.lev = max.reads + dopen/3
  x = seq(c2$start, c2$end, 1)
  y = rep(coord.lev, length(x))
  points(x,y,type="l",lwd=4,col="black")

  # plot search window
  abline(v=search$start)
  abline(v=search$end)

  # set region for peak analysis
  clip = search$xy
  x = c(-add.to.end:0)
  y0 = exp(-x/tau.dist)
  y = ( (y0-min(y0))/(max(y0)-min(y0)) ) * (frac.max-frac.min) + frac.min
  frac.dist = y[add.to.end - clip + 1]
  ds = frac.dist * (search$end - search$start)
  if(strand == "+"){
    abline(v=search$start+ds, lty=2)
  }
  if(strand == "-"){
    abline(v=search$end-ds, lty=2)
  }

} # tts.plot



############################################################################################
## gene.end.plot
############################################################################################

#' Plot the results of transcription terminination site identification
#'
#' This function can be used to evaluate the performance of get.TTS().
#' We also recommend visualizing the results of data-driven gene annotations using a genome browser.
#' @param bed Bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param gene A specific gene for plotting
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param bp.bin The interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh The TTS is defined as pk.thresh percent of max peak of the spline fit
#' @param knot.div The number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
#' @param cnt.thresh Read count number below which the TTS is not evaluated
#' @param knot.thresh Minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end The maximal length of the search region (bp)
#' @param tau.dist Distance constant for the exponential defining the region for peak detection
#' @param frac.max Maximal fraction of gene end region for peal detection
#' @param frac.min Minimal fraction of gene end region for peal detection
#' @return
#' A single plot is generated with the binned reads (circles), a fitted curve (red), 
#' and the location of the identified TTS (vertical line).
#' @export
#' @examples
#' gene.end.plot(bed=bed.for.tts.eval, gene="Pparg",
#'               bw.plus=bw.plus, bw.minus=bw.minus,
#'               bp.bin=bp.bin, add.to.end=add.to.end, knot.div=knot.div,
#'               pk.thresh=pk.thresh, knot.thresh=knot.thresh,
#'               cnt.thresh=cnt.thresh, tau.dist=tau.dist,
#'               frac.max=frac.max, frac.min=frac.min)
gene.end.plot = function(bed=NULL, gene=NULL, bw.plus=NULL, bw.minus=NULL, bp.bin=NULL,
                         pk.thresh=0.1, knot.div=40, knot.thresh=5,
                         cnt.thresh=5, add.to.end=50000,
                         tau.dist=10000, frac.max=1, frac.min=0.2){

  # gene information
  dat = bed[bed$gene==gene[1],]
  strand = dat$strand

  # set region for peak analysis
  clip = dat$xy
  x = c(-add.to.end:0)
  y0 = exp(-x/tau.dist)
  y = ( (y0-min(y0))/(max(y0)-min(y0)) ) * (frac.max-frac.min) + frac.min
  frac.dist = y[add.to.end - clip + 1]

  # divide the gene into bins and obtain counts
  if( dat$end - dat$start < bp.bin ){print("coords error")}
  if(strand=="+"){
    bw = bw.plus
    starts = seq(dat$start, (dat$end-bp.bin), by=bp.bin)
    ends = seq((dat$start+bp.bin), dat$end, by=bp.bin)
    lim = min( length(starts), length(ends) )
    bed.map = data.frame(dat$chr, starts[1:lim], ends[1:lim],
                         stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    bed.map$end[nrow(bed.map)] = dat$end
  }
  if(strand=="-"){
    bw = bw.minus
    ends = seq(dat$end, (dat$start+bp.bin), by=-bp.bin)
    starts = seq((dat$end-bp.bin), dat$start, by=-bp.bin)
    lim = min( length(starts), length(ends) )
    bed.map = data.frame(dat$chr, starts[1:lim], ends[1:lim],
                         stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    bed.map$start[nrow(bed.map)] = dat$start
  }
  counts = bed.region.bpQuery.bigWig(bw, bed.map)
  if(length(counts) < cnt.thresh){print("counts are too low")}


  # fit a spline to the count profile
  # get peaks from the spline fit and estimate the best TTS
  # TTS defined by pk.thresh percent of max peak of the spline fit
  nknots = round( length(bed.map$end) / knot.div )
  if(nknots < knot.thresh){nknots = knot.thresh}
  if(strand=="+"){
    spl = smooth.spline(x=bed.map$end, y=counts, all.knots=FALSE, nknots=nknots)
    
    srch = c(min(bed.map$start), max(bed.map$end))
    dsrch = srch[2] - srch[1]
    fsrch = srch[1] + dsrch * frac.dist
    dd = spl$x - fsrch
    isearch = which(abs(dd) == min(abs(dd)))[1]
    
    pks = findpeaks(spl$y[1:isearch], nups=0)
    indpk = which( pks[,1] == max(pks[,1]) )
    ind.tts = which(spl$y[pks[indpk,2]:length(spl$y)] < pk.thresh * max(pks[,1]))
    if(length(ind.tts) == 0){
      TTS = spl$x[length(spl$x)]
    } else {
      ind.tts = min( which(spl$y[pks[indpk,2]:length(spl$y)] <
                             pk.thresh * max(pks[,1])) ) + pks[indpk,2] - 1
      TTS = spl$x[ind.tts]
    }
    pred = predict(spl,bed.map$end)
    plot(bed.map$end,counts, ylab=paste0(gene," counts"),xlab="coord (bp)");
    points(pred$x,pred$y,type="l",col="red"); abline(v=TTS)
  } # plus
  if(strand=="-"){
    spl = smooth.spline(x=-bed.map$start, y=counts, all.knots=FALSE, nknots=nknots)
    
    srch = c(min(bed.map$start), max(bed.map$end))
    dsrch = srch[2] - srch[1]
    fsrch = srch[2] - dsrch * frac.dist
    dd = spl$x + fsrch
    isearch = which(abs(dd) == min(abs(dd)))[1]
    
    pks = findpeaks(spl$y[1:isearch], nups=0)
    indpk = which( pks[,1] == max(pks[,1]) )
    ind.tts = which(spl$y[pks[indpk,2]:length(spl$y)] < pk.thresh * max(pks[,1]))
    if(length(ind.tts)==0){
      TTS = -spl$x[length(spl$x)]
    } else {
      ind.tts = min( which(spl$y[pks[indpk,2]:length(spl$y)] <
                             pk.thresh * max(pks[,1])) ) + pks[indpk,2] - 1
      TTS = -spl$x[ind.tts]
    }
    pred = predict(spl,-bed.map$start)
    plot(-bed.map$start,counts,ylab=paste0(gene," counts"), xlab="coord (bp)");
    points(pred$x,pred$y,type="l",col="red"); abline(v=-TTS)
  } # minus

} # gene.end.plot



############################################################################################
## TTS.boundry.match
############################################################################################

#' Check for the fraction of identified TTSs that match the search boundry
#'
#' This function considers information in both coords and bed.for.tts.eval (see vignette)
#' to identify the fraction of TTSs that match the boundary of the search region.
#' @param coords A bed6 frame with inferred coordinates (from get.TTS())
#' @param bed.for.tts.eval A bed6 frame with TTS search regions (from get.end.intervals())
#' @return
#' A fraction [0,1].
#' @export
#' @examples
#' frac.match = TTS.boundary.match(coords=coords, bed.for.tts.eval=bed.for.tts.eval)
TTS.boundary.match = function(coords=NULL, bed.for.tts.eval=NULL){
  gene.ends.plus = coords %>% filter(strand=="+") %>% select(end) %>% data.matrix %>% as.numeric
  gene.ends.minus = coords %>% filter(strand=="-") %>% select(start) %>% data.matrix %>% as.numeric
  tts.eval.plus = bed.for.tts.eval %>% filter(strand=="+") %>% select(end) %>% data.matrix %>% as.numeric
  tts.eval.minus = bed.for.tts.eval %>% filter(strand=="-") %>% select(start) %>% data.matrix %>% as.numeric
  plus.match = gene.ends.plus[gene.ends.plus %in% tts.eval.plus]
  minus.match = gene.ends.minus[gene.ends.minus %in% tts.eval.minus]
  pos = length(plus.match)
  neg = length(minus.match)
  frac.match = (pos+neg)/nrow(bed.for.tts.eval)
  return(frac.match)
} # TTS.boundry.match




############################################################################################
## TTS.boundary.clip
############################################################################################

#' Look at whether TTSs identified at the search boundry were clipped.
#'
#' This function considers information in both coords and bed.for.tts.eval (see vignette)
#' to identify the fractions of search boundary TTSs that were clipped and unclipped.
#' @param coords A bed6 frame with inferred coordinates (from get.TTS())
#' @param bed.for.tts.eval A bed6 frame with TTS search regions (from get.end.intervals())
#' @return
#' A list with two fractions [0,1], frac.clip and frac.noclip,
#' and a vector of clip distances (clip.tts) for genes with boundary TTSs.
#' @export
#' @examples
#' # look at whether TTSs identified at the search boundry had clipped boundries
#' frac.bound.clip = TTS.boundry.clip(coords=coords, bed.for.tts.eval=bed.for.tts.eval)
#' frac.bound.clip$frac.clip
#' frac.bound.clip$frac.noclip
#'
#' # plot didtribution of clip distances for genes with boundary TTSs
#' hist(frac.bound.clip$clip.tts,main="",xlab="bp clipped for boundry genes")
#'
TTS.boundary.clip = function(coords=NULL, bed.for.tts.eval=NULL){
  clip.tts.plus = merge(coords, bed.for.tts.eval,by="gene") %>%
    filter(strand.x=="+",end.x==end.y)
  clip.tts.minus = merge(coords, bed.for.tts.eval,by="gene") %>%
    filter(strand.x=="-",start.x==start.y)
  clip.tts = rbind(clip.tts.plus, clip.tts.minus)
  clip.tts0 = clip.tts %>% filter(xy.x==0)
  n.boundry.genes = nrow(clip.tts0)
  frac.bound.noclip = length(which(clip.tts$xy.x == 0)) / nrow(clip.tts)
  frac.bound.clip = length(which(clip.tts$xy.x > 0)) / nrow(clip.tts)
  return(list(frac.clip=frac.bound.clip, frac.noclip=frac.bound.noclip, clip.tts=clip.tts$xy.x))
} # TTS.boundary.clip



############################################################################################
## chrom.size.filter
############################################################################################

#' Filter identified gene coordinates for chromosome sizes.
#'
#' First acquire a chrom.sizes file (e.g., http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes; see vignette).
#' If any ends are greater than the species chromosome size, reduce the end to the size limit.
#' Any negative starts will be set to 0.
#' @param coords A bed6 frame with inferred coordinates (from get.TTS())
#' @param chrom.sizes A two column frame with chromosomes and their respective sizes in bp
#' @return
#' A list with a bed6 frame with corrected coordinates (bed)
#' and a log documenting which corrections, if any, were made (log).
#' @export
#' @examples
#' coords.filt = chrom.size.filter(coords=coords, chrom.sizes=chrom.sizes)
#' head(coords.filt$bed)
#' coords.filt$log
#'
chrom.size.filter = function(coords=NULL, chrom.sizes=NULL){
  log.fix = c()
  for(ii in 1:length(unique(coords$chr))){
    ind.chr = which(coords$chr == unique(coords$chr)[ii])
    ind.fix = which(coords$end[ind.chr] >
                      chrom.sizes$size[chrom.sizes$chr == unique(coords$chr)[ii]])
    if(length(ind.fix) > 0){
      coords$end[ind.chr][ind.fix] =
        chrom.sizes$size[chrom.sizes$chr == unique(coords$chr)[ii]]
      log.fix = c(log.fix, rep(unique(coords$chr)[ii],length(ind.fix)))
    }
  }
  return(list(bed=coords, log=log.fix))
} # chrom.size.filter


############################################################################################
## eval.tss
############################################################################################

#' Evaluate TSS inference performance
#'
#' This function was designed to evaluate the results of out TSS identification analysis.
#' The user should input inferred and largest interval coordinates (bed1 and bed2, respectively).
#' We specify a region centered on the TSS. We obtain read counts in bins that span this window. 
#' We sort the genes based on the bin with the
#' maximal reads and we scale the data to the interval (0,1) for visualization. 
#' We compute the distances between the TSS and
#' the bin with the max reads within the specified window. Note that this function calls
#' TSS.count.dist() for the TSS distance analysis.
#' @param bed1 A bed frame with inferred TSSs
#' @param bed2 A bed frame with reference TSSs
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param window Region size, centered on the TSS for analysis
#' @param bp.bin The interval will be separated into adjacent bins of this size
#' @param bin.thresh Minimal value of reads in the bin with the max count, below this value the transcript will be removed
#' @param fname A .pdf file name for output plots (if NULL, no plot)
#' @return
#' A list of two lists and a plot (optional, specifiy fname for the plot to be printed to file).
#' Each list contains three vector elements: raw, scaled, and dist.
#' $tss.dists.inf cooresponds to the TSS data (bed1) and
#' $tss.dists.lng cooresponds to the largest interval data (bed2).
#' All outputs are organized based on TSS position (left-upstream).
#' raw = binned raw counts. scaled = binned scaled counts (0,1),
#' with genes sorted based on the distance between the TSS and the max read interval.
#' dist = upper bound distances ( min(|dists|) = bp.bin ).
#' See examples and vignette for more details. The first row of plots includes dirstibutions
#' of distances between the TSS and the maximal reads for the inferred TSSs (left) and largest interval TSSs (right).
#' The second row of plots shows the relationship between the distances from TSSs to max reads and the read depth
#' in the region of read count evaluation (left inferred, right largest interval).
#' The heatmaps show read profiles (horizontal axis, distances centered in the middle of the window; vertical axis, distinct genes).
#' The genes are sorted based on the inferred TSSs to the max reads (left), and the largest interval distances (right) are organized
#' based on the gene order for the inferred TSSs.
#' @export
#' @examples
#' #
eval.tss = function(bed1=NULL, bed2=NULL, bw.plus=NULL, bw.minus=NULL,
                    window=NULL, bp.bin=NULL, fname=NULL, bin.thresh=10){

  # look at read distribution around identified TSSs
  tss.dists.inf = TSS.count.dist(bed=bed1, bw.plus=bw.plus, bw.minus=bw.minus,
                             window=window, bp.bin=bp.bin, bin.thresh=bin.thresh)


  # look at read distribution around 'long gene' annotation TSSs
  tss.dists.lng = TSS.count.dist(bed=bed2, bw.plus=bw.plus, bw.minus=bw.minus,
                                 window=window, bp.bin=bp.bin, bin.thresh=bin.thresh)
  
  # filter to get common genes in both sets
  all(rownames(tss.dists.inf$raw) == rownames(tss.dists.inf$scaled))
  all(rownames(tss.dists.inf$raw) == names(tss.dists.inf$dist))
  common.genes = intersect(names(tss.dists.inf$dist),names(tss.dists.lng$dist))
  
  tss.dists.inf$raw = tss.dists.inf$raw[rownames(tss.dists.inf$raw) %in% common.genes,]
  tss.dists.inf$scaled = tss.dists.inf$scaled[rownames(tss.dists.inf$scaled) %in% common.genes,]
  tss.dists.inf$dist = tss.dists.inf$dist[names(tss.dists.inf$dist) %in% common.genes]
  
  tss.dists.lng$raw = tss.dists.lng$raw[rownames(tss.dists.lng$raw) %in% common.genes,]
  tss.dists.lng$scaled = tss.dists.lng$scaled[rownames(tss.dists.lng$scaled) %in% common.genes,]
  tss.dists.lng$dist = tss.dists.lng$dist[names(tss.dists.lng$dist) %in% common.genes]

  if(is.null(fname)==FALSE){

    id.inds = sapply(names(tss.dists.inf$dist),function(x){which(names(tss.dists.lng$dist)==x)})
    bk = seq(-window/2, window/2, bp.bin)
    pdf(fname); par(mfrow=c(2,2))
    hist(tss.dists.inf$dist,main="inferred",xlab="dist from TSS to read max (bp)",
         breaks=bk,col="black")
    hist(tss.dists.lng$dist,main="largest",xlab="dist from TSS to read max (bp)",
         breaks=bk,col="black")
    plot(rowSums(tss.dists.inf$raw), tss.dists.inf$dist, xlab="window read depth",
         ylab="dist from TSS to read max (bp)")
    plot(rowSums(tss.dists.lng$raw), tss.dists.lng$dist, xlab="window read depth",
         ylab="dist from TSS to read max (bp)")
    aheatmap(tss.dists.inf$scaled,Colv=NA,Rowv=NA,color=brewer.pal(9,"Greys"))
    aheatmap(tss.dists.lng$scaled[id.inds,],Colv=NA,Rowv=NA,color=brewer.pal(9,"Greys"))
    dev.off()

  } # plots

  return(list(tss.dists.inf=tss.dists.inf, tss.dists.lng=tss.dists.lng))

} # eval.tss


############################################################################################
## TSS.count.dist
############################################################################################

#' Get distances from identified transcription start sites (TSSs) to the nearest regions with peak reads
#'
#' This function was designed to evaluate the results of out TSS identification analysis. We specify a region centered
#' on the TSS. We obtain read counts in bins that span this window. We sort the genes based on the bin with the
#' maximal reads and we scale the data to the interval (0,1) for visualization. We compute the distances between the TSS and
#' the bin with the max reads within the specified window.
#' @param bed A bed frame for TSS evaluation
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param window Region size, centered on the TSS for analysis
#' @param bp.bin The interval will be separated into adjacent bins of this size
#' @param bin.thresh Minimal value of reads in the bin with the max count, below this value the transcript will be removed
#' @return
#' A list with three vector elements: raw, scaled, and dist.
#' All outputs are organized based on TSS position (left-upstream).
#' raw = binned raw counts. scaled = binned scaled counts (0,1).
#' dist = upper bound distances ( min(|dists|) = bp.bin ).
#' See examples and vignette for more details.
#' @export
#' @examples
#' see eval.tss()
TSS.count.dist = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, window=NULL, bp.bin=NULL, bin.thresh=NULL){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")
  bed.plus.window = bed.plus
  bed.minus.window = bed.minus

  # set windows around TSS
  bed.plus.window$start = bed.plus$start - window/2
  bed.plus.window$end = bed.plus$start + window/2
  bed.minus.window$start = bed.minus$end - window/2
  bed.minus.window$end = bed.minus$end + window/2

  # set read data frames
  bp.relative = c(seq(-window/2,-1,by=bp.bin), seq(bp.bin, window/2,by=bp.bin))
  counts.plus = matrix(0,nrow(bed.plus.window),window/bp.bin) %>%
    as.data.frame(stringsAsFactors=F)
  counts.minus = matrix(0,nrow(bed.minus.window),window/bp.bin) %>%
    as.data.frame(stringsAsFactors=F)
  rownames(counts.plus) = bed.plus.window$gene
  rownames(counts.minus) = bed.minus.window$gene
  names(counts.plus) = names(counts.minus) = bp.relative

  # plus counts
  for(ii in 1:nrow(bed.plus.window)){
    start = seq(bed.plus.window$start[ii],
                (bed.plus.window$end[ii]-bp.bin), by=bp.bin)
    end = seq((bed.plus.window$start[ii]+bp.bin),
              bed.plus.window$end[ii], by=bp.bin)
    bed.map = data.frame(bed.plus.window$chr[ii], start, end,
                         stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    counts.plus[ii,] = bed.region.bpQuery.bigWig(bw.plus, bed.map)
  }

  # minus counts
  for(ii in 1:nrow(bed.minus.window)){
    end = seq(bed.minus.window$end[ii],
              (bed.minus.window$start[ii]+bp.bin), by=-bp.bin)
    start = seq((bed.minus.window$end[ii]-bp.bin),
                bed.minus.window$start[ii], by=-bp.bin)
    bed.map = data.frame(bed.minus.window$chr[ii],
                         start, end, stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    counts.minus[ii,] = bed.region.bpQuery.bigWig(bw.minus, bed.map)
  }

  # aggregate data and sort by max bin reads
  counts = rbind(counts.plus, counts.minus)
  bin.max = apply(counts,1,function(x){which(x==max(x))[1]}) %>% unlist
  
  # filter data for low read counts
  # bin.thresh = 10
  ind.rem = which(counts[bin.max] < bin.thresh)
  counts = counts[-ind.rem,]
  bin.max = bin.max[-ind.rem]
  
  # sort the data
  indices = order(bin.max,decreasing=FALSE)
  counts.raw = counts[indices,]

  # scale each row to (0,1)
  counts.scl = apply(counts.raw,1,function(x){
    (x-min(x))/(max(x)-min(x))}) %>% t
  
  # dist from TSS to peak, upper bound
  dist0 = (bin.max[indices] - (window/bp.bin)/2) * bp.bin
  dist0[dist0 <= 0] = dist0[dist0 <= 0] - bp.bin

  # return data
  return(list(raw=counts.raw, scaled=counts.scl, dist=dist0))

} # TSS.count.dist


############################################################################################
## TTS.count.dist
############################################################################################

#' Get read counts around transcription termination sites (TTSs)
#'
#' We set a window around the TTS and segment the window into bins. To sort the read data,
#' we compute cumulative counts of reads, and we sort based on the bin at which a specified percentage of
#' the reads are found. We also take a ratio of gene counts downstream / upstream of the
#' TTS with regions within an interval.
#' @param bed A bed frame for TSS evaluation
#' @param bw.plus Plus strand bigWig data
#' @param bw.minus Minus strand bigWig data
#' @param window Region size, centered on the TTS for analysis
#' @param bp.bin The interval will be separated into adjacent bins of this size
#' @param frac.max Fraction of cumulative distribution for sorting entries are sorted by indices min( cumsum(x)/sum(x) ) < frac.max
#' @param ratio.region Number of bp on either side of the TTS to evaluate the ratio entries
#' @return
#' A list with three elements: raw, scaled, and ratio.
#' All outputs are organized based on TTS position (left-upstream, see frac.max).
#' raw = binned raw counts.
#' scaled = binned scaled counts (0,1).
#' ratio = downstream(ratio.region) / upstream(ratio.region).
#' @export
#' @examples
#' # look at read distribution around identified TTSs
#' window = 1000
#' bp.bin = 10
#' frac.max = 0.8
#' ratio.region = 300
#' tts.dists = TTS.count.dist(bed=coords, bw.plus=bw.plus, bw.minus=bw.minus,
#'      window=window, bp.bin=bp.bin, frac.max=frac.max,
#'      ratio.region=ratio.region)
TTS.count.dist = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, window=NULL,
                          bp.bin=NULL, frac.max=NULL, ratio.region=NULL){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")
  bed.plus.window = bed.plus
  bed.minus.window = bed.minus

  # set windows around TTS
  bed.plus.window$start = bed.plus$end - window/2
  bed.plus.window$end = bed.plus$end + window/2
  bed.minus.window$start = bed.minus$start - window/2
  bed.minus.window$end = bed.minus$start + window/2

  # set read data frames
  bp.relative = c(seq(-window/2,-1,by=bp.bin), seq(bp.bin, window/2,by=bp.bin))
  counts.plus = matrix(0,nrow(bed.plus.window),window/bp.bin) %>%
    as.data.frame(stringsAsFactors=F)
  counts.minus = matrix(0,nrow(bed.minus.window),window/bp.bin) %>%
    as.data.frame(stringsAsFactors=F)
  rownames(counts.plus) = bed.plus.window$gene
  rownames(counts.minus) = bed.minus.window$gene
  names(counts.plus) = names(counts.minus) = bp.relative

  # plus counts
  for(ii in 1:nrow(bed.plus.window)){
    start = seq(bed.plus.window$start[ii],
                (bed.plus.window$end[ii]-bp.bin), by=bp.bin)
    end = seq((bed.plus.window$start[ii]+bp.bin),
              bed.plus.window$end[ii], by=bp.bin)
    bed.map = data.frame(bed.plus.window$chr[ii],
                         start, end, stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    counts.plus[ii,] = bed.region.bpQuery.bigWig(bw.plus, bed.map)
  }

  # minus counts
  for(ii in 1:nrow(bed.minus.window)){
    end = seq(bed.minus.window$end[ii],
              (bed.minus.window$start[ii]+bp.bin), by=-bp.bin)
    start = seq((bed.minus.window$end[ii]-bp.bin),
                bed.minus.window$start[ii], by=-bp.bin)
    bed.map = data.frame(bed.minus.window$chr[ii],
                         start, end, stringsAsFactors=F)
    names(bed.map) = c("chr","start","end")
    counts.minus[ii,] = bed.region.bpQuery.bigWig(bw.minus, bed.map)
  }

  # aggregate data, sort by cumulative sum, and scale (zscore)
  counts = rbind(counts.plus, counts.minus)
  cumulatives = apply(counts,1,function(x){cumsum(x)/sum(x)}) %>% t
  cum.cut = apply(cumulatives,1,function(x){
    ind = which(x>frac.max)
    if(length(ind)>0){return(min(ind))}else{return(window/bp.bin)}
  })
  indices = order(cum.cut,decreasing=FALSE)
  counts.raw = counts[indices,]
  counts.scl = apply(counts.raw,1,function(x){
    (x-min(x))/(max(x)-min(x))}) %>% t

  # get counts ratio - after TTS / before TTS
  reg.bins = ratio.region / bp.bin
  cntr.bin = (window/bp.bin)/2
  bins.up = (cntr.bin-reg.bins+1):cntr.bin
  bins.dn = (cntr.bin+1):(cntr.bin+reg.bins)
  ratio = apply(counts.raw,1,function(x){
    sum(x[bins.dn])/sum(x[bins.up]) })

  # return data
  return(list(raw=counts.raw, scaled=counts.scl, ratio=ratio))

} # TTS.count.dist
