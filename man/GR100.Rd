\name{GR100}
\alias{GR100}
\alias{GR}
\docType{data}
\title{
Peaks of a ChIP-seq experiment
}
\description{
The data-set \code{GR100} contains 100 enriched regions on chromosome 18 of the transcription 
factor c-Myc on murine cells, obtained by calling peaks with MACS [Zhang et al., 2008] on 
a ChIP-Seq experiment. The genomic coordinates of the peaks
are stored in the \linkS4class{GRanges} object \code{GR100}. 
}
\usage{data("GR100")}
\format{
  \linkS4class{GRanges} object with 100 ranges and 0 metadata columns.
  }
\references{
Zhang et al., 2008.  Model-based Analysis of ChIP-Seq (MACS). Genome Biology, vol. 9 (9) pp. R137.
}
\examples{
data(GR100)
GR
}
\keyword{datasets}
