
############################################################################################
## get.largest.interval
############################################################################################

#' Get the largest interval for each gene, given multiple TSS and TTS annotations
#'
#' The input bed6 file can be derived from a gencode annotation file,
#' as described in the vignette
#' @param bed bed6 frame with comprehensive gene annotations, defaults to NULL
#' @return a bed6 frame with the largest annotation for each gene
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
## read.count.end
############################################################################################

#' Select a regions containing the gene ends
#'
#' Select a fraction of the annotated gene end and consider an additional number of base pairs beyond
#' the gene end within which to count reads
#' @param bed a bed6 frame
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param fraction.end fraction of the annotated gene end (0,1)
#' @param add.to.end number of base pairs beyond the gene end within which to count reads
#' @return
#' a list with counts and desity.
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
#' We set a range around each annotated gene start within which to search for a region of peak read
#' density. Such regions of peak read density occur at ’pause sites’ that are typically 20-80 bp from
#' the TSS. Within each region around an annotated gene start, we translate a sliding
#' window throughout to find the sub-region with maximal density. We determine the window with the
#' maximal density out of all regions considered and define the TSS as the upstream boundry of this
#' window, minus a shift of a specified amount. The output of this function should be processed using the function apply.TSS.coords().
#' @param bed bed6 file with comprehensive gene annotations
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param bp.range range to be centered on each gene start for evaluating read counts.
#'        Note that the value should be divisible by 2
#' @param bp.delta number of reads to move incrementally - sliding window interval
#' @param bp.bin bin size for evaluating read counts within a given region
#' @param delta.tss amount by which to shift the identied TSS upstream of the identified interval
#' @param cnt.thresh read count number below which the TTS is not evaluated
#' @return
#' A vector of indentified transcription start site coordinates. Note that strand and chromosome information are not provided here.
#' @export
#' @examples
#' # select the TSS for each gene
#' bp.range = 1200
#' bp.delta = 10
#' bp.bin = 50
#' delta.tss = 50
#' cnt.thresh = 5
#' TSS.gene = get.TSS(bed=dat0, bw.plus=bw.plus, bw.minus=bw.minus,
#'                   bp.range=bp.range, bp.delta=bp.delta, bp.bin=bp.bin,
#'                   delta.tss=delta.tss, cnt.thresh=cnt.thresh)
get.TSS = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, bp.range=NULL,
                   bp.delta=NULL, bp.bin=NULL, delta.tss=NULL, cnt.thresh=NULL){

  if(bp.range %% 2 != 0){stop("select a value for bp.range that is divisible by 2")}

  # loop through each unique gene
  # map reads a range for each annotated TSS
  genes = unique(bed$gene)
  tss = rep(0,length(genes))
  for(ii in 1:length(genes)){

    # get gene of interest and strand
    ind.gene = which(bed$gene == genes[ii])
    strand.gene = unique(bed$strand[ind.gene])[1]

    # get all sets of intervals around each annotated TSS
    if(strand.gene == "-"){ starts = unique(bed$end[ind.gene]) }
    if(strand.gene == "+"){ starts = unique(bed$start[ind.gene]) }
    ranges = sapply(starts,function(x){ c(x-bp.range/2, x+bp.range/2) })
    intervals = c()
    for(jj in 1:ncol(ranges)){
      x = ranges[,jj]
      starts = seq(x[1], x[2]-bp.bin, by=bp.delta)
      ends = seq(x[1]+bp.bin, x[2], by=bp.delta)
      intervals = rbind(intervals, cbind(starts,ends))
    } # jj

    # get counts for each interval of interest
    bed.map = cbind(unique(bed$chr[ind.gene]), intervals) %>%
      as.data.frame(stringsAsFactors=FALSE)
    bed.map[,2:3] = apply(bed.map[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
    names(bed.map) = c("chr","start","end")
    if(strand.gene == "-"){ counts = bed.region.bpQuery.bigWig(bw.minus, bed.map) }
    if(strand.gene == "+"){ counts = bed.region.bpQuery.bigWig(bw.plus, bed.map) }

    # select the best TSS
    ind.max.cnt = which(counts == max(counts))[1]
    if(counts[ind.max.cnt] < cnt.thresh){next}
    if(strand.gene == "-"){ tss[ii] = bed.map$end[ind.max.cnt] + delta.tss + 1 }
    if(strand.gene == "+"){ tss[ii] = bed.map$start[ind.max.cnt] - delta.tss }

  } # ii

  names(tss) = genes
  return(tss)

} # get.TSS


############################################################################################
## apply.TSS.coords
############################################################################################

#' Apply previously identified TSSs to an existing gene annotation
#'
#' This function will apply transcription start site (TSS) coordinates to a bed6 frame to incorporate the output from
#' get.TSS() into an existing bed frame.
#' @param bed bed6 file with comprehensive gene annotations
#' @param tss TSS estimates, output from get.TSS().
#'        Note that the TSS genes must be a subset of the genes in the bed6 data.
#' @return
#' A bed6 file with estimated TSSs incorporated
#' @export
#' @examples
#' # apply TSS coordinates to the expression-filtered long gene annotations
#' bed.long.filtered.tss = apply.TSS.coords(bed=bed.long.filtered, tss=TSS.gene)
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
## gene.overlaps
############################################################################################

#' Documentation of gene overlaps
#'
#' This function identifies gene overlaps that do not involve identical coordinates.
#' @param bed bed6 file with processed gene annotations
#' @return
#' A list with has.start.inside and is.a.start.inside. has.start.inside = bed6 file that documents genes in which there are starts
#'         of other genes. is.a.start.inside = bed6 file that documents genes that start inside other
#'         genes.
#' @export
#' @examples
#' # run overlap analysis
#' overlap.data = gene.overlaps( bed = bed.long.filtered2.tss )
#' has.start.inside = overlap.data$has.start.inside
#' is.a.start.inside = overlap.data$is.a.start.inside
#' dim(has.start.inside)
#' dim(is.a.start.inside)
gene.overlaps = function(bed=NULL){

  # separate by strand
  bed.plus = bed %>% filter(strand == "+")
  bed.minus = bed %>% filter(strand == "-")

  # data for documenting overlaps
  has.start.inside = c()
  is.a.start.inside = c()

  # filter plus by start sites to document end/start overlaps
  # loop through each chromosome
  chr.plus = unique(bed.plus$chr)
  for(ii in 1:length(chr.plus)){
    inds.chr = which(bed.plus$chr == chr.plus[ii])

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      start = bed.plus$start[inds.chr][jj]
      end = bed.plus$end[inds.chr][jj]
      ind.overlap = which(bed.plus$start[inds.chr]>start &
                            bed.plus$start[inds.chr]<end)
      if(length(ind.overlap) > 0){
        overlap.genes = paste0(bed.plus$gene[inds.chr][ind.overlap], collapse=",")

        # this gene has a start inside
        bed.orig = bed.plus[inds.chr[jj],]
        bed.orig$xy = paste0(overlap.genes, " starts inside this gene")
        has.start.inside = rbind(has.start.inside, bed.orig)

        # this gene starts inside another gene
        bed.orig = bed.plus[inds.chr[ind.overlap],]
        bed.orig$xy = paste0("this gene starts inside ", bed.plus$gene[inds.chr][jj])
        is.a.start.inside = rbind(is.a.start.inside, bed.orig)
      }
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
      ind.overlap = which(bed.minus$end[inds.chr]<start & bed.minus$end[inds.chr]>end)
      if(length(ind.overlap) > 0){
        overlap.genes = paste0(bed.minus$gene[inds.chr][ind.overlap], collapse=",")

        # this gene has a start inside
        bed.orig = bed.minus[inds.chr[jj],]
        bed.orig$xy = paste0(overlap.genes, " starts inside this gene")
        has.start.inside = rbind(has.start.inside, bed.orig)

        # this gene starts inside another gene
        bed.orig = bed.minus[inds.chr[ind.overlap],]
        bed.orig$xy = paste0("this gene starts inside ", bed.minus$gene[inds.chr][jj])
        is.a.start.inside = rbind(is.a.start.inside, bed.orig)
      }
    } ## jj

  } # ii, minus

  # generate output
  return( list(has.start.inside=has.start.inside,
               is.a.start.inside=is.a.start.inside) )

} # gene.overlaps


############################################################################################
## adjacent.gene.tss
############################################################################################


#' Function to identify transcription start sites (TSSs) for adjacent gene pairs
#'
#' We have identified adjacent gene pairs with overlaps based on manual analyses.
#' We empirically define TSSs for these genes by binning the region spanned by both genes,
#' fitting smooth spline curves to the binned read counts, and identifying the two largest peaks separated
#' by a specified distance. We set a bin size and apply the constraint that the identified
#' 'pause' peaks must be a given distance apart. The search region includes a specified distance upstream of
#' the upstream-most gene. We shift the identified peaks upstream. For the spline fits, we set the
#' number of knots to the number of bins divided by specified setting. See also get.TSS().
#' @param fix.genes frame with upstream and downstream genes
#' @param bed.long long gene annotations, see get.largest.interval()
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param bp.bin the interval will be separated into adjacent bins of this size
#' @param shift.up look upstream of the long gene start by this amount to search for a viable TSS
#' @param delta.tss amount by which to shift the identied TSS upstream of the identified interval
#' @param knot.div the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number
#' @param knot.thresh minimum number of knots
#' @param diff.tss minimum distance between TSSs
#' @param fname file name (.pdf) for output plots
#' @return
#' A vector of TSSs, along with a plot. The titles for the plots indicate the upstream gene, the downstream gene, and the strand.
#' For genes on the plus strand, the upstream gene is on the left. For the minus strand, the upstream gene is on the right.
#' The red line denotes the spline fit and the vertical lines indicate the TSSs.
#' The plot is intended for diagnostic purposes.
#' @export
#' @examples
#' # get the start sites for adjacent gene pairs
#' bp.bin = 10
#' knot.div = 40
#' shift.up = 100
#' delta.tss = 50
#' diff.tss = 2000
#' TSS.adjacent = adjacent.gene.tss(fix.genes=fix.genes, bed.long=bed.long,
#'                               bw.plus=bw.plus, bw.minus=bw.minus,
#'                               knot.div=knot.div, bp.bin=bp.bin,
#'                               shift.up=shift.up, delta.tss=delta.tss,
#'                               diff.tss=diff.tss, fname="adjacentTSS.pdf")
adjacent.gene.tss = function(fix.genes=NULL, bed.long=NULL, bw.plus=NULL,
                             bw.minus=NULL, bp.bin=NULL, shift.up=NULL,
                             delta.tss=NULL, knot.div=4, knot.thresh=5,
                             diff.tss=1000, fname="adjacentTSS.pdf"){

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
  return(TSS)
} # adjacent.gene.tss



############################################################################################
## adjacent.gene.tts
############################################################################################

#' Function to identify gene ends (TTSs) for adjacent gene pairs
#'
#' We set the TTSs for the upstream genes to a given distance from the starts of the downstream genes for the
#' adjacent gene pairs. See also adjacent.gene.tss().
#' @param fix.genes frame with upstream and downstream genes for adjacent gene pairs
#' @param bed.long comprehensive gene annotations with large intervals defined by get.largest.interval()
#' @param TSS.adjacent TSSs for adjancent genes, identified from adjacent.gene.tss()
#' @param dist.from.start distance between the TTS of the upstream gene and the TSS of the downstream gene
#' @return
#' A bed6 frame with TTSs for upstream genes in adjacent gene pairs
#' @export
#' @examples
#' # get the end sites for adjacent upstream genes, integrate with TSSs
#' dist.from.start = 50
#' adjacent.tts.up = adjacent.gene.tts(fix.genes=fix.genes, bed.long=bed.long,
#'                                  TSS.adjacent=TSS.adjacent, dist.from.start=dist.from.start)
#' tss = TSS.adjacent$TSS
#' names(tss) = TSS.adjacent$gene
#' adjacent.tts.up = apply.TSS.coords(bed=adjacent.tts.up, tss=tss)
adjacent.gene.tts = function(fix.genes=NULL, bed.long=NULL, TSS.adjacent=NULL,
                             dist.from.start=NULL){

  # get upstream TTSs
  gene.ends.up = bed.long[bed.long$gene %in% fix.genes$upstream,]
  for(ii in 1:nrow(fix.genes)){
    upstr = TSS.adjacent$TSS[TSS.adjacent$gene == fix.genes$upstream[ii]]
    downstr = TSS.adjacent$TSS[TSS.adjacent$gene == fix.genes$downstream[ii]]
    strand = bed.long$strand[bed.long$gene == fix.genes$upstream[ii]]
    if(strand == "+"){gene.ends.up$end[gene.ends.up$gene==fix.genes$upstream[ii]] =
      as.numeric(downstr) - dist.from.start}
    if(strand == "-"){gene.ends.up$start[gene.ends.up$gene==fix.genes$upstream[ii]] =
      as.numeric(downstr) + dist.from.start}
  }

  # return results
  return(gene.ends.up)

} # adjacent.gene.tts


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
#' @param bed processed bed6 file with one entry per gene and identified TSSs
#' @param add.to.end number of bases to add to the end of each gene to define the search region
#' @param fraction.end fraction of the gene annotation to consider at the end of the gene
#' @param dist.from.start the maximal allowable distance between the end of one gene and the start of the next
#' @return
#' A bed6 for gene end evaluation
#' @export
#' @examples
#' # get intervals for TTS evaluation
#' add.to.end = 100000
#' fraction.end=0.1
#' dist.from.start=50
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

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      old.start = bed.plus.orig$start[inds.chr]
      new.start = bed.plus.new$start[inds.chr][jj]
      new.end = bed.plus.new$end[inds.chr][jj]
      ind.overlap = which(old.start>new.start & old.start<new.end)
      if(length(ind.overlap) > 0){
        ind.closest = which(old.start[ind.overlap] == min(old.start[ind.overlap]))[1]
        old = bed.plus.big$end[inds.chr][jj]
        new = old.start[ind.overlap][ind.closest] - dist.from.start
        bed.plus.new$end[inds.chr][jj] = new
        bed.plus.new$xy[inds.chr][jj] = old - new
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

    # check for start sites within the annotoation of a given gene
    for(jj in 1:length(inds.chr)){
      old.start = bed.minus.orig$end[inds.chr]
      new.start = bed.minus.new$end[inds.chr][jj]
      new.end = bed.minus.new$start[inds.chr][jj]
      ind.overlap = which(old.start>new.end & old.start<new.start)
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
  bed.out = rbind(bed.plus.new, bed.minus.new)
  return(bed.out)

} # get.end.intervals


############################################################################################
## get.gene.end
############################################################################################

#' Defining gene ends
#'
#'Given regions within which to search for TTSs, we opperationally define the TTSs by binning the gene
#' end regions, counting reads within the bins, fitting smooth spline curves to the bin counts, and
#' detecting points at which the curves decay towards zero. We applied the
#' constraint that there must be a specified number of bases in the gene end interval, otherwise the
#' TTS analysis is not applied. Similarly, if the number of knots identified is too low, then we set the
#' number of knots to a specified threshold. For this analysis, we set a sub-region at the beginning of
#' the gene end region and identify the maximal peak from the spline fit. Then we identify the point
#' at which the spline fit decays to a threshold level of of the peak level.
#' We reasoned that the sub-region should be largest for genes with
#' the greatest numbers of clipped bases, because such cases occur when the conventional
#' gene ends are proximal to identified TSSs, and we should include these entire regions for analysis of
#' the TTS. Similarly, we reasoned that for genes with substantially less clipped bases, and
#' correspondingly larger gene end regions with greater potential for observing enhancers or divergent transcripts,
#' the sub-regions should be smaller sections of the upstream-most gene end region. We use
#' an exponential model to define the sub-regions.
#' The user should use the output of get.end.intervals() as an input to this function.
#' Note that gene ends will be set to zeros if there are no counts or if the interval is too small.
#' This function calls get.TTS().
#' @param bed bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param bp.bin the interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh the TTS is defined as pk.thresh percent of max peak of the spline fit
#' @param knot.div the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
#' @param cnt.thresh read count number below which the TTS is not evaluated
#' @param knot.thresh minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end the maximal length of the search region (bp)
#' @param tau.dist distance constant for the exponential defining the region for peak detection
#' @param frac.max maximal fraction of gene end region for peal detection
#' @param frac.min minimal fraction of gene end region for peal detection
#' @return
#' A bed6 file with TTS estimates incorporated, zeros are present if the TTS could not be estimated.
#' The output is a list including the bed frame along with several metrics (see example below and vignette).
#' @export
#' @examples
#' # identify gene ends
#' add.to.end = max(bed.for.tts.eval$xy)
#' knot.div = 40
#' pk.thresh = 0.02
#' bp.bin = 50
#' knot.thresh = 5
#' cnt.thresh = 5
#' tau.dist = 10000
#' frac.max = 1
#' frac.min = 0.2
#' gene.ends = get.gene.end(bed=bed.for.tts.eval, bw.plus=bw.plus, bw.minus=bw.minus,
#'                      bp.bin=bp.bin, add.to.end=add.to.end, knot.div=knot.div,
#'                      pk.thresh=pk.thresh, knot.thresh=knot.thresh,
#'                      cnt.thresh=cnt.thresh, tau.dist=tau.dist,
#'
#' # get metrics
#' minus.lowcount = length(gene.ends$minus.lowcount)
#' plus.lowcount = length(gene.ends$plus.lowcount)
#' minus.knotmod = length(gene.ends$minus.knotmod)
#' plus.knotmod = length(gene.ends$plus.knotmod)
#' ends = gene.ends$bed
get.gene.end = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, bp.bin=NULL,
                        pk.thresh=0.1, knot.div=40, knot.thresh=5,
                        cnt.thresh=5, add.to.end=50000,
                        tau.dist=10000, frac.max=1, frac.min=0.2){

  # separate genes by strand
  bed.plus = bed %>% filter(strand=="+")
  bed.minus = bed %>% filter(strand=="-")

  # get TTS for each gene
  TTS.plus = get.TTS(bed=bed.plus,
                     bw=bw.plus, bp.bin=bp.bin,
                     pk.thresh=pk.thresh,
                     knot.div=knot.div,
                     cnt.thresh=cnt.thresh,
                     knot.thresh=knot.thresh,
                     add.to.end=add.to.end,
                     tau.dist=tau.dist,
                     frac.max=frac.max,
                     frac.min=frac.min)
  TTS.minus = get.TTS(bed=bed.minus,
                      bw=bw.minus, bp.bin=bp.bin,
                      pk.thresh=pk.thresh,
                      knot.div=knot.div,
                      cnt.thresh=cnt.thresh,
                      knot.thresh=knot.thresh,
                      add.to.end=add.to.end,
                      tau.dist=tau.dist,
                      frac.max=frac.max,
                      frac.min=frac.min)

  # aggregate data and return
  bed.plus$end = TTS.plus$tts
  bed.minus$start = TTS.minus$tts
  bed.out = rbind(bed.plus, bed.minus)
  return(list(bed=bed.out,
              minus.lowcount=TTS.minus$lowcount,
              plus.lowcount=TTS.plus$lowcount,
              minus.knotmod=TTS.minus$knotmod,
              plus.knotmod=TTS.plus$knotmod))

} # get.gene.end


############################################################################################
## get.TTS
############################################################################################

#' Estimate the transcription termination sites (TTSs)
#'
#' This function estimates the TTSs, this function is called by get.gene.end() and is not designed to be run on its own.
#' @param bed bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param bw plus or minus strand bigWig data
#' @param bp.bin the interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh the TTS is defined as pk.thresh percent of max peak of the spline fit
#' @param knot.div the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number
#' @param cnt.thresh read count number below which the TTS is not evaluated
#' @param knot.thresh minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end the maximal length of the search region (bp)
#' @param tau.dist distance constant for the exponential defining the region for peak detection
#' @param frac.max maximal fraction of gene end region for peal detection
#' @param frac.min minimal fraction of gene end region for peal detection
#' @return
#' A vector of TTS coordinates
#' @export
#' @examples
#' See get.gene.end()
get.TTS = function(bed=NULL, bw=NULL, bp.bin=NULL, pk.thresh=NULL, knot.div=NULL,
                   cnt.thresh=NULL, knot.thresh=NULL, add.to.end=NULL,
                   tau.dist=NULL, frac.max=NULL, frac.min=NULL){

  # check strand
  strand = unique(bed$strand)
  if(length(strand) > 1){stop("issue with strand separation in get.TTS()")}

  # loop through each gene and select a TTS coordinate
  tts = rep(0,length(bed$end))
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

    # fit a spline to the count profile
    # get peaks from the spline fit and estimate the best TTS
    # TTS defined by pk.thresh percent of max peak of the spline fit
    nknots = round( length(bed.map$end) / knot.div )
    if(nknots < knot.thresh){knotmod = c(knotmod,ii); nknots = knot.thresh}
    if(strand=="+"){
      spl = smooth.spline(x=bed.map$end, y=counts, all.knots=FALSE, nknots=nknots)
      pks = findpeaks(spl$y[1:round(frac.dist*length(spl$y))], nups=0)
      indpk = which( pks[,1] == max(pks[,1]) )
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
      pks = findpeaks(spl$y[1:round(frac.dist*length(spl$y))], nups=0)
      indpk = which( pks[,1] == max(pks[,1]) )
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
} # get.TTS


############################################################################################
## gene.end.plot
############################################################################################

#' Plot the results of transcription terminination site identification
#'
#' This function can be used to evaluate the performance of get.gene.end() and get.TTS().
#' We also recommend visualizing the results of data-driven gene annotations using a genome browser.
#' @param bed bed6 gene coordinates with gene end intervals for estimation of the TTS
#' @param gene a specific gene for plotting
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param bp.bin the interval at the gene end will be separated into adjacent bins of this size
#' @param pk.thresh the TTS is defined as pk.thresh percent of max peak of the spline fit
#' @param knot.div the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
#' @param cnt.thresh read count number below which the TTS is not evaluated
#' @param knot.thresh minimum number of knots, if knot.div gives a lower number, use this
#' @param add.to.end the maximal length of the search region (bp)
#' @param tau.dist distance constant for the exponential defining the region for peak detection
#' @param frac.max maximal fraction of gene end region for peal detection
#' @param frac.min minimal fraction of gene end region for peal detection
#' @return
#' A single plot is generated
#' @export
#' @examples
#' gene.end.plot(bed=bed.for.tts.eval, gene="Aamp",
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
    pks = findpeaks(spl$y[1:round(frac.dist*length(spl$y))], nups=0)
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
    pks = findpeaks(spl$y[1:round(frac.dist*length(spl$y))], nups=0)
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
## TSS.count.dist
############################################################################################

#' Get distances from identified transcription start sites (TSSs) to the nearest regions with peak reads
#'
#' This function was designed to evaluate the results of out TSS identification analysis. We specify a region centered
#' on the TSS. We obtain read counts in bins that span this window. We sort the genes based on the bin with the
#' maximal reads and we scale the data to the interval (0,1) for visualization. We compute the distances between the TSS and
#' the bin with the max reads within the specified window.
#' @param bed = bed file for TSS evaluation
#' @param bw.plus = plus strand bigWig data
#' @param bw.minus = minus strand bigWig data
#' @param window = region size, centered on the TSS for analysis
#' @param bp.bin = the interval will be separated into adjacent bins of this size
#' @return
#' A list with three vector elements: raw, scaled, and dist.
#' All outputs are organized based on TSS position (left-upstream).
#' raw = binned raw counts. scaled = binned scaled counts (0,1).
#' dist = upper bound distances ( min(|dists|) = bp.bin ).
#' See examples and vignette for more details.
#' @export
#' @examples
#' # look at read distribution around identified TSSs
#' bp.bin = 10
#' window = 1000
#' tss.dists = TSS.count.dist(bed=bed.tss.tts, bw.plus=bw.plus, bw.minus=bw.minus,
#'                            window=window, bp.bin=bp.bin)
#'
#' # look at read distribution around 'long gene' annotation TSSs
#' tss.dists.lng = TSS.count.dist(bed=bed.long[bed.long$gene %in% names(tss.dists$dist),],
#'                                bw.plus=bw.plus, bw.minus=bw.minus,
#'                                window=window, bp.bin=bp.bin)
TSS.count.dist = function(bed=NULL, bw.plus=NULL, bw.minus=NULL, window=NULL,
                          bp.bin=NULL){

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
#' @param bed bed frame for TSS evaluation
#' @param bw.plus plus strand bigWig data
#' @param bw.minus minus strand bigWig data
#' @param window region size, centered on the TTS for analysis
#' @param bp.bin the interval will be separated into adjacent bins of this size
#' @param frac.max fraction of cumulative distribution for sorting entries are sorted by indices min( cumsum(x)/sum(x) ) < frac.max
#' @param ratio.region number of bp on either side of the TTS to evaluate the ratio entries
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
#' tts.dists = TTS.count.dist(bed=bed.tss.tts, bw.plus=bw.plus, bw.minus=bw.minus,
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
