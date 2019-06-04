
############################################################################################
## single.overlaps
############################################################################################

#' Assign identifiers to TUs with single gene overlaps
#'
#' This function will assign identifiers to TUs that overlapped with single genes.
#' Trusted gene annotations are considered in terms of their degree of overlap with
#' TUs identified in an unbiased manner. Marginal regions of TUs outside of gene overlaps
#' are given generic identifiers.
#' @param overlaps frame with bed intersect of inferred TUs (col1-6)
#' with annotated TUs (col7-12) and overlaping bps (col14)
#' @param tss.thresh = number of bp a TU beginning can be off from an annotation in
#' order to be assigned that annotation
#' @param delta.tss max distance between an upstream gene end and downstream gene start
#' @param delta.tts max difference distance between and annotated gene end and the start
#' of a downstream gene before an intermediate TU id is assigned
#' @return
#' A list with bed data (bed) and counts for each class (cnt1-4).
#' bed = data with chr, start, end, id, class, and strand.
#' cnt1 = count of TU from class 1 overlaps
#' @export
#' @examples
#' sing.hmm.rows = overlap.tu[!duplicated(overlap.tu[,c(1:3,6)]) &
#'      !duplicated(overlap.tu[,c(1:3,6)], fromLast=TRUE),]
#' overlaps = sing.hmm.rows
#' names(overlaps)[1:6] = c("infr.chr","infr.start","infr.end",
#'                          "infr.gene","infr.xy","infr.strand")
#' tss.thresh = 200
#' delta.tss = 50
#' delta.tts = 1000
#' class1234 = single.overlaps(overlaps=overlaps, tss.thresh=tss.thresh,
#'                             delta.tss=delta.tss, delta.tts=delta.tts)
#' new.ann.sing = class1234$bed
single.overlaps = function(overlaps=NULL, tss.thresh=NULL,
                           delta.tss=NULL, delta.tts=NULL){

  if(nrow(overlaps) != nrow(unique(overlaps))){stop("all input rows are not unique")}

  # identify unique TUs and set data objects
  class1.counts = class2.counts = class3.counts = class4.counts = 0
  class1 = class2 = class3 = class4 = c()

  # loop through each unique TU and address each class
  for(ii in 1:nrow(overlaps)){


    # class 1 - the annotation is within the TU
    if(overlaps$infr.start[ii] < overlaps$ann.start[ii] &
       overlaps$infr.end[ii] > overlaps$ann.end[ii]){
      fr = overlaps[ii,]
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=1)
      class1 = rbind(class1, new)
      class1.counts = class1.counts + 1
    } # class 1


    # class 2 - the TU is within the annotation
    if(overlaps$infr.start[ii] > overlaps$ann.start[ii] &
       overlaps$infr.end[ii] < overlaps$ann.end[ii]){
      fr = overlaps[ii,]
      fr$infr.start = fr$ann.start
      fr$infr.end = fr$ann.end
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=2)
      class2 = rbind(class2, new)
      class2.counts = class2.counts + 1
    } # class 2


    # class 3 - TU upstream overlap
    if(overlaps$infr.start[ii] < overlaps$ann.start[ii] &
       overlaps$infr.end[ii] < overlaps$ann.end[ii] &
       overlaps$infr.strand[ii] == "+"){
      fr = overlaps[ii,]
      fr$infr.end = fr$ann.end
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=3)
      class3 = rbind(class3, new)
      class3.counts = class3.counts + 1
    } # class 3 - plus
    if(overlaps$infr.start[ii] > overlaps$ann.start[ii] &
       overlaps$infr.end[ii] > overlaps$ann.end[ii] &
       overlaps$infr.strand[ii] == "-"){
      fr = overlaps[ii,]
      fr$infr.start = fr$ann.start
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=3)
      class3 = rbind(class3, new)
      class3.counts = class3.counts + 1
    } # class 3 - minus


    # class 4 - TU downstream overlap
    if(overlaps$infr.start[ii] > overlaps$ann.start[ii] &
       overlaps$infr.end[ii] > overlaps$ann.end[ii] &
       overlaps$infr.strand[ii] == "+"){
      fr = overlaps[ii,]
      #if(fr$infr.strand=="+"){stop(ii)}
      fr$infr.start = fr$ann.start
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=4)
      class4 = rbind(class4, new)
      class4.counts = class4.counts + 1
    } # class 4 - plus
    if(overlaps$infr.start[ii] < overlaps$ann.start[ii] &
       overlaps$infr.end[ii] < overlaps$ann.end[ii] &
       overlaps$infr.strand[ii] == "-"){
      fr = overlaps[ii,]
      #if(fr$infr.strand=="-"){stop(ii)}
      fr$infr.end = fr$ann.end
      new = single.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                  delta.tss=delta.tss,
                                  delta.tts=delta.tts, cl=4)
      class4 = rbind(class4, new)
      class4.counts = class4.counts + 1
    } # class 4 - minus

  } # ii

  # rename unassigned TUs and return data:
  class1$id[class1$id==0] = paste0("tu_class1_",c(1:length(which(class1$id==0))))
  class2$id[class2$id==0] = paste0("tu_class2_",c(1:length(which(class2$id==0))))
  class3$id[class3$id==0] = paste0("tu_class3_",c(1:length(which(class3$id==0))))
  class4$id[class4$id==0] = paste0("tu_class4_",c(1:length(which(class4$id==0))))
  out = rbind(class1, class2, class3, class4)
  class.cnt = (class1.counts+class2.counts+class3.counts+class4.counts)
  if(nrow(overlaps) == class.cnt){
    return(list(bed=out, cnt1=class1.counts, cnt2=class2.counts,
                cnt3=class3.counts, cnt4=class4.counts))
  } else {
    return(list(bed=out, cnt1=class1.counts, cnt2=class2.counts,
                cnt3=class3.counts, cnt4=class4.counts))
    stop("some entries could no be assigned")
  }

} # single.overlaps


############################################################################################
## single.overlap.assign
############################################################################################

#' Function to assign identifiers for class1 overlap TUs
#'
#' This function is called inside of single.overlaps()
#' and is not recommended to be used on its own.
#' @param fr frame with full bedtools overlap data for a class1
#' TU with multiple internal annotations
#' @param tss.thresh number of bp a TU beginning can be off
#' from an annotation in order to be assigned that annotation
#' @param delta.tss max distance between an upstream gene end and downstream gene start
#' @param cl multi-overlap class
#' @return
#' A bed data frame with chr, start, end, id, class, and strand.
#' Note that id = 0 for regions of large TUs that do no match input annotations.
#' @export
#' @examples
#' See single.overlaps()
single.overlap.assign = function(fr=NULL, tss.thresh=NULL, delta.tss=NULL,
                                 delta.tts=NULL, cl=NULL){

  # determine strand and set data object
  ov.class = cl
  chr = fr$infr.chr[1]
  strand = fr$infr.strand[1]
  out = data.frame(matrix(0,1,6),stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","id","class","strand")

  # assing TU identifiers for the plus strand
  if(strand == "+"){

    # assign the first subset
    tss.dist = fr$ann.start - fr$infr.start
    if(tss.dist > tss.thresh){
      out$chr[1] = chr
      out$start[1] = fr$infr.start
      out$end[1] = fr$ann.start - delta.tss
      out$class[1] = ov.class
      out$strand[1] = strand
      end = fr$ann.end
      start = fr$ann.start
      id = fr$ann.gene
      out = rbind(out, c(chr, start, end, id, ov.class, strand))
    } else {
      out$chr[1] = chr
      out$start[1] = fr$infr.start
      out$end[1] = fr$ann.end
      out$id[1] = fr$ann.gene
      out$class[1] = ov.class
      out$strand[1] = strand
    }

    # assign the last subset
    last.dist = fr$infr.end - as.numeric(out$end[nrow(out)])
    if(last.dist > delta.tts){
      start = as.numeric(out$end[nrow(out)]) + delta.tss
      end = fr$infr.end
      out = rbind(out, c(chr, start, end, 0, ov.class, strand))
    } else {
      out$end[nrow(out)] = fr$infr.end
    } # last subset

  } # plus

  # assing TU identifiers for the minus strand
  if(strand == "-"){

    # assign the first subset
    tss.dist = fr$infr.end - fr$ann.end
    if(tss.dist > tss.thresh){
      out$chr[1] = chr
      out$end[1] = fr$infr.end[1]
      out$start[1] = fr$ann.end + delta.tss
      out$class[1] = ov.class
      out$strand[1] = strand
      end = fr$ann.end
      start = fr$ann.start
      id = fr$ann.gene
      out = rbind(out, c(chr, start, end, id, ov.class, strand))
    } else {
      out$chr[1] = chr
      out$end[1] = fr$infr.end[1]
      out$start[1] = fr$ann.start
      out$id[1] = fr$ann.gene
      out$class[1] = ov.class
      out$strand[1] = strand
    } # first subset

    # assign the last subset
    last.dist = as.numeric(out$start[nrow(out)]) - fr$infr.start
    if(last.dist > delta.tts){
      end = as.numeric(out$start[nrow(out)]) - delta.tss
      start = fr$infr.start[1]
      out = rbind(out, c(chr, start, end, 0, ov.class, strand))
    } else {
      out$start[nrow(out)] = fr$infr.start
    } # last subset

  } # minus

  # return data:
  return(out)

} # single.overlap.assign


############################################################################################
## multi.overlaps
############################################################################################

#' Assign identifiers to TUs with multiple gene overlaps
#'
#' This function will assign identifiers to TUs that overlapped with multiple genes.
#' Trusted gene annotations are considered in terms of their degree of overlap with
#' TUs identified in an unbiased manner. Marginal regions of TUs outside of gene overlaps
#' are given generic identifiers.
#' @param overlaps frame with bed intersect of inferred TUs (col1-6) with annotated TUs (col7-12) and overlaping bps (col14)
#' @param tss.thresh number of bp a TU beginning can be off from an annotation in order to be assigned that annotation
#' @param delta.tss max distance between an upstream gene end and downstream gene start
#' @param delta.tts max difference distance between and annotated gene end and the start of a downstream gene before an intermediate TU id is assigned
#' @return
#' A list with bed data (bed) and counts for each class (cnt5-8).
#' bed = data with chr, start, end, id, class, and strand.
#' cnt5 = count of TU from class 5 overlaps
#' @export
#' @examples
#' dup.hmm.rows = overlap.tu[duplicated(overlap.tu[,c(1:3,6)]) |
#'     duplicated(overlap.tu[,c(1:3,6)], fromLast=TRUE),]
#' dup.full = dup.hmm.rows
#' names(dup.full)[1:6] = c("infr.chr","infr.start","infr.end",
#'                          "infr.gene","infr.xy","infr.strand")
#' nrow(dup.full[,1:3] %>% unique)
#' nrow(dup.full %>% select(ann.gene) %>% unique)
#' tss.thresh = 200
#' delta.tss = 50
#' delta.tts = 1000
#' class5678 = multi.overlaps(overlaps=dup.full, tss.thresh=tss.thresh,
#'                            delta.tss=delta.tss, delta.tts=delta.tts)
#' new.ann.mult = class5678$bed
multi.overlaps = function(overlaps=NULL, tss.thresh=NULL,
                          delta.tss=NULL, delta.tts=NULL){

  # identify unique TUs and set data objects
  dup.uniq = unique(overlaps[,c(1:3,6)])
  class5.counts = class6.counts = class7.counts = class8.counts = 0
  class5 = class6 = class7 = class8 = c()

  # loop through each unique TU and address each class
  for(ii in 1:nrow(dup.uniq)){

    # get the frame with all overlaps for a given TU
    inds = which(dup.full$infr.chr == dup.uniq$infr.chr[ii] &
                   dup.full$infr.start == dup.uniq$infr.start[ii] &
                   dup.full$infr.end == dup.uniq$infr.end[ii] &
                   dup.full$infr.strand == dup.uniq$infr.strand[ii])
    fr = dup.full[inds,]

    # get ranges for overlapping genes to identify overlap class
    gene.range = c(min(fr$ann.start), max(fr$ann.end))

    # class 5 - all annotations internal
    if(gene.range[1] > dup.uniq$infr.start[ii] &
       gene.range[2] < dup.uniq$infr.end[ii]){
      new = multi.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                 delta.tss=delta.tss,
                                 delta.tts=delta.tts, cl=5)
      class5 = rbind(class5, new)
      class5.counts = class5.counts + 1
    } # class 5

    # class 6 - upstream overhang
    if( (gene.range[1] < dup.uniq$infr.start[ii] &
         gene.range[2] < dup.uniq$infr.end[ii] &
         dup.uniq$infr.strand[ii]=="+") |
        (gene.range[2] > dup.uniq$infr.end[ii] &
         gene.range[1] > dup.uniq$infr.start[ii] &
         dup.uniq$infr.strand[ii]=="-") ){
      if(dup.uniq$infr.strand[ii]=="+"){
        ind.upstream = which(fr$ann.start == min(fr$ann.start))
        fr$infr.start = fr$ann.start[ind.upstream]
      }
      if(dup.uniq$infr.strand[ii]=="-"){
        ind.upstream = which(fr$ann.end == max(fr$ann.end))
        fr$infr.end = fr$ann.end[ind.upstream]
      }
      new = multi.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                 delta.tss=delta.tss,
                                 delta.tts=delta.tts, cl=6)
      class6 = rbind(class6, new)
      class6.counts = class6.counts + 1
    } # class6

    # class 7 - downstream overhang
    if( (gene.range[1] > dup.uniq$infr.start[ii] &
         gene.range[2] > dup.uniq$infr.end[ii] &
         dup.uniq$infr.strand[ii]=="+") |
        (gene.range[2] < dup.uniq$infr.end[ii] &
         gene.range[1] < dup.uniq$infr.start[ii] &
         dup.uniq$infr.strand[ii]=="-") ){
      if(dup.uniq$infr.strand[ii]=="+"){
        ind.dnstream = which(fr$ann.end == max(fr$ann.end))
        fr$infr.end = fr$ann.end[ind.dnstream]
      }
      if(dup.uniq$infr.strand[ii]=="-"){
        ind.dnstream = which(fr$ann.start == min(fr$ann.start))
        fr$infr.start = fr$ann.start[ind.dnstream]
      }
      new = multi.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                 delta.tss=delta.tss,
                                 delta.tts=delta.tts, cl=7)
      class7 = rbind(class7, new)
      class7.counts = class7.counts + 1

    } # class 7

    # class 8 - up and downstream overhangs
    if( gene.range[1] < dup.uniq$infr.start[ii] &
        gene.range[2] > dup.uniq$infr.end[ii] ){
      fr$infr.start = fr$ann.start[ which(fr$ann.start == min(fr$ann.start)) ]
      fr$infr.end = fr$ann.end[ which(fr$ann.end == max(fr$ann.end)) ]
      new = multi.overlap.assign(fr=fr, tss.thresh=tss.thresh,
                                 delta.tss=delta.tss,
                                 delta.tts=delta.tts, cl=8)
      class8 = rbind(class8, new)
      class8.counts = class8.counts + 1
    } # class 8

  } # ii

  # rename unassigned TUs and return data:
  class5$id[class5$id==0] = paste0("tu_class5_",c(1:length(which(class5$id==0))))
  class6$id[class6$id==0] = paste0("tu_class6_",c(1:length(which(class6$id==0))))
  class7$id[class7$id==0] = paste0("tu_class7_",c(1:length(which(class7$id==0))))
  class8$id[class8$id==0] = paste0("tu_class8_",c(1:length(which(class8$id==0))))
  out = rbind(class5, class6, class7, class8)
  class.cnt = (class5.counts+class6.counts+class7.counts+class8.counts)
  if(nrow(dup.uniq) == class.cnt){
    return(list(bed=out, cnt5=class5.counts, cnt6=class6.counts,
                cnt7=class7.counts, cnt8=class8.counts))
  } else {
    return(list(bed=out, cnt5=class5.counts, cnt6=class6.counts,
                cnt7=class7.counts, cnt8=class8.counts))
    stop("some entries could no be assigned")
  }

} # multi.overlaps


############################################################################################
## multi.overlap.assign
############################################################################################

#' Function to assign identifiers for class5 overlap TUs
#'
#' This function is called inside of multiple.overlaps()
#' and is not recommended to be used on its own.
#' @param fr = frame with full bedtools overlap data for a class5 TU
#' with multiple internal annotations
#' @param tss.thresh = number of bp a TU beginning can be off
#' from an annotation in order to be assigned that annotation
#' @param delta.tss = max distance between an upstream gene end and
#' downstream gene start cl = multi-overlap class
#' @param delta.tts max difference distance between and annotated gene
#' end and the start of a downstream gene before an intermediate TU id is assigned
#' @param cl multi-overlap class
#' @return
#' A bed data with chr, start, end, id, class, and strand.
#' Note that id = 0 for regions of large TUs that do no match input annotations.
#' @export
#' @examples
#' See multiple.overlaps()
multi.overlap.assign = function(fr=NULL, tss.thresh=NULL, delta.tss=NULL,
                                delta.tts=NULL, cl=NULL){

  # determine strand and set data object
  ov.class = cl
  chr = fr$infr.chr[1]
  strand = fr$infr.strand[1]
  out = data.frame(matrix(0,1,6),stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","id","class","strand")

  # assing TU identifiers for the plus strand
  if(strand == "+"){

    # assign the first subset
    ind.upstream = which(fr$ann.start == min(fr$ann.start))
    tss.dist = fr$ann.start[ind.upstream] - fr$infr.start[1]
    if(tss.dist > tss.thresh){
      out$chr[1] = chr
      out$start[1] = fr$infr.start[1]
      out$end[1] = fr$ann.start[ind.upstream] - delta.tss
      out$class[1] = ov.class
      out$strand[1] = strand
      end = fr$ann.end[ind.upstream]
      start = fr$ann.start[ind.upstream]
      id = fr$ann.gene[ind.upstream]
      out = rbind(out, c(chr, start, end, id, ov.class, strand))
    } else {
      out$chr[1] = chr
      out$start[1] = fr$infr.start[1]
      out$end[1] = fr$ann.end[ind.upstream]
      out$id[1] = fr$ann.gene[ind.upstream]
      out$class[1] = ov.class
      out$strand[1] = strand
    } # first subset

    # assign the intermediate subsets
    fr.new = fr[-ind.upstream,]
    for(jj in 1:(nrow(fr)-1)){
      ind.up = which(fr.new$ann.start == min(fr.new$ann.start))
      gene.dist = fr.new$ann.start[ind.up] - as.numeric(out$end[nrow(out)])
      if(gene.dist > delta.tts){
        start = as.numeric(out$end[nrow(out)]) + delta.tss
        end = fr.new$ann.start[ind.up] - delta.tss
        out = rbind(out, c(chr, start, end, 0, ov.class, strand))
        end = fr.new$ann.end[ind.up]
        start = fr.new$ann.start[ind.up]
        id = fr.new$ann.gene[ind.up]
        out = rbind(out, c(chr, start, end, id, ov.class, strand))
        fr.new = fr.new[-ind.up,]
      } else {
        out$end[nrow(out)] = fr.new$ann.start[ind.up] - delta.tss
        end = fr.new$ann.end[ind.up]
        start = fr.new$ann.start[ind.up]
        id = fr.new$ann.gene[ind.up]
        out = rbind(out, c(chr, start, end, id, ov.class, strand))
        fr.new = fr.new[-ind.up,]
      }
    } # jj, intermediate subsets

    # assign the last subset
    last.dist = fr$infr.end[1] - as.numeric(out$end[nrow(out)])
    if(last.dist > delta.tts){
      start = as.numeric(out$end[nrow(out)]) + delta.tss
      end = fr$infr.end[1]
      out = rbind(out, c(chr, start, end, 0, ov.class, strand))
    } else {
      out$end[nrow(out)] = fr$infr.end[1]
    } # last subset

  } # plus

  # assing TU identifiers for the minus strand
  if(strand == "-"){

    # assign the first subset
    ind.upstream = which(fr$ann.end == max(fr$ann.end))
    tss.dist = fr$infr.end[1] - fr$ann.end[ind.upstream]
    if(tss.dist > tss.thresh){
      out$chr[1] = chr
      out$end[1] = fr$infr.end[1]
      out$start[1] = fr$ann.end[ind.upstream] + delta.tss
      out$class[1] = ov.class
      out$strand[1] = strand
      end = fr$ann.end[ind.upstream]
      start = fr$ann.start[ind.upstream]
      id = fr$ann.gene[ind.upstream]
      out = rbind(out, c(chr, start, end, id, ov.class, strand))
    } else {
      out$chr[1] = chr
      out$end[1] = fr$infr.end[1]
      out$start[1] = fr$ann.start[ind.upstream]
      out$id[1] = fr$ann.gene[ind.upstream]
      out$class[1] = ov.class
      out$strand[1] = strand
    } # first subset

    # assign the intermediate subsets
    fr.new = fr[-ind.upstream,]
    for(jj in 1:(nrow(fr)-1)){
      ind.up = which(fr.new$ann.end == max(fr.new$ann.end))
      gene.dist = as.numeric(out$start[nrow(out)]) - fr.new$ann.end[ind.up]
      if(gene.dist > delta.tts){
        end = as.numeric(out$start[nrow(out)]) - delta.tss
        start = fr.new$ann.end[ind.up] + delta.tss
        out = rbind(out, c(chr, start, end, 0, ov.class, strand))
        end = fr.new$ann.end[ind.up]
        start = fr.new$ann.start[ind.up]
        id = fr.new$ann.gene[ind.up]
        out = rbind(out, c(chr, start, end, id, ov.class, strand))
        fr.new = fr.new[-ind.up,]
      } else {
        out$start[nrow(out)] = fr.new$ann.end[ind.up] + delta.tss
        end = fr.new$ann.end[ind.up]
        start = fr.new$ann.start[ind.up]
        id = fr.new$ann.gene[ind.up]
        out = rbind(out, c(chr, start, end, id, ov.class, strand))
        fr.new = fr.new[-ind.up,]
      }
    } # jj, intermediate subsets

    # assign the last subset
    last.dist = as.numeric(out$start[nrow(out)]) - fr$infr.start[1]
    if(last.dist > delta.tts){
      end = as.numeric(out$start[nrow(out)]) - delta.tss
      start = fr$infr.start[1]
      out = rbind(out, c(chr, start, end, 0, ov.class, strand))
    } else {
      out$start[nrow(out)] = fr$infr.start[1]
    } # last subset

  } # minus

  # return data:
  return(out)

} # multi.overlap.assign

