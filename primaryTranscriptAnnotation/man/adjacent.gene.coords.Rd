% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_ann.R
\name{adjacent.gene.coords}
\alias{adjacent.gene.coords}
\title{Function to identify coordinates for adjacent gene pairs}
\usage{
adjacent.gene.coords(fix.genes = NULL, bed.long = NULL, exon1 = NULL,
  bw.plus = NULL, bw.minus = NULL, bp.bin = NULL, shift.up = NULL,
  dist.from.start = 50, delta.tss = NULL, knot.div = 4,
  knot.thresh = 5, diff.tss = 1000, pause.bp = 120,
  fname = "adjacentSplines.pdf")
}
\arguments{
\item{fix.genes}{A frame with upstream and downstream genes}

\item{bed.long}{Long gene annotations, see get.largest.interval()}

\item{exon1}{A bed6 frame with first exon gene annotations}

\item{bw.plus}{Plus strand bigWig data}

\item{bw.minus}{Minus strand bigWig data}

\item{bp.bin}{The interval will be separated into adjacent bins of this size}

\item{shift.up}{Look upstream of the long gene start by this amount to search for a viable TSS}

\item{dist.from.start}{Distance between the TTS of the upstream gene and the TSS of the downstream gene}

\item{delta.tss}{Amount by which to shift the identied TSS upstream of the identified interval}

\item{knot.div}{The number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number}

\item{knot.thresh}{Minimum number of knots}

\item{diff.tss}{Minimum distance between TSSs}

\item{pause.bp}{Number of bp by which a pause site should be considered downstrean of a TSS}

\item{fname}{File name (.pdf) for output plots}
}
\value{
A frame with the adjusted coordinates for the input gene set, along with a plot. 
The titles for the plots indicate the upstream gene, the downstream gene, and the strand.
For genes on the plus strand, the upstream gene is on the left. 
For the minus strand, the upstream gene is on the right.
The red line denotes the spline fit and the vertical lines indicate the pause peaks.
The plot is intended for diagnostic purposes.
}
\description{
Adjacent gene pairs can be identified based on manual analyses, 
despite overlaps in the largest interval coordinates.
We empirically define pause sites for these genes by binning the region spanned by both genes,
fitting smooth spline curves to the binned read counts, and identifying the two largest peaks separated
by a specified distance. We set a bin size and apply the constraint that the identified
pause peaks must be a given distance apart. 
For the spline fits, we set the
number of knots to the number of bins divided by specified setting.
We identify the TSSs as the exon 1 coordinates closest to the pause site peaks.
}
\examples{
# get the start sites for adjacent gene pairs
fix.genes = rbind(fix.genes.id, fix.genes.ov)
bp.bin = 5
knot.div = 40
shift.up = 100
delta.tss = 50
diff.tss = 1000
dist.from.start  = 50
adjacent.coords = adjacent.gene.coords(fix.genes=fix.genes, bed.long=TSS.gene,
                                       exon1=gencode.firstExon,
                                       bw.plus=bw.plus, bw.minus=bw.minus,
                                       knot.div=knot.div, bp.bin=bp.bin,
                                       shift.up=shift.up, delta.tss=delta.tss,
                                       dist.from.start=dist.from.start,
                                       diff.tss=diff.tss, fname="adjacentSplines.pdf")

}
