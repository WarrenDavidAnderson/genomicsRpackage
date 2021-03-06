% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_ann.R
\name{TTS.count.dist}
\alias{TTS.count.dist}
\title{Get read counts around transcription termination sites (TTSs)}
\usage{
TTS.count.dist(bed = NULL, bw.plus = NULL, bw.minus = NULL,
  window = NULL, bp.bin = NULL, frac.max = NULL,
  ratio.region = NULL)
}
\arguments{
\item{bed}{A bed frame for TSS evaluation}

\item{bw.plus}{Plus strand bigWig data}

\item{bw.minus}{Minus strand bigWig data}

\item{window}{Region size, centered on the TTS for analysis}

\item{bp.bin}{The interval will be separated into adjacent bins of this size}

\item{frac.max}{Fraction of cumulative distribution for sorting entries are sorted by indices min( cumsum(x)/sum(x) ) < frac.max}

\item{ratio.region}{Number of bp on either side of the TTS to evaluate the ratio entries}
}
\value{
A list with three elements: raw, scaled, and ratio.
All outputs are organized based on TTS position (left-upstream, see frac.max).
raw = binned raw counts.
scaled = binned scaled counts (0,1).
ratio = downstream(ratio.region) / upstream(ratio.region).
}
\description{
We set a window around the TTS and segment the window into bins. To sort the read data,
we compute cumulative counts of reads, and we sort based on the bin at which a specified percentage of
the reads are found. We also take a ratio of gene counts downstream / upstream of the
TTS with regions within an interval.
}
\examples{
# look at read distribution around identified TTSs
window = 1000
bp.bin = 10
frac.max = 0.8
ratio.region = 300
tts.dists = TTS.count.dist(bed=coords, bw.plus=bw.plus, bw.minus=bw.minus,
     window=window, bp.bin=bp.bin, frac.max=frac.max,
     ratio.region=ratio.region)
}
