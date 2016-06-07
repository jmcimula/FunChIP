\name{pileup_peak}
\docType{methods}
\alias{pileup_peak-method}
\alias{pileup_peak}
\alias{pileup_peak,GRanges-method}
\title{Computing read counts on a \linkS4class{GRanges} 
object.
}
\description{Given a \linkS4class{GRanges} object and the path of a .bam
file, it creates the corresponding pileup, containing the 
read counts on each nucleotide of the peaks
of the \linkS4class{GRanges} object. Reads can be extended to a length \code{d}, which is an estimate of the length of the sequencing fragment. See the function  \link[=compute_fragments_length]{compute_fragments_length} for details. For each peak this method creates a
vector containing these counts, i.e. the coverage function for the extended reads along the whole peak.
}
\usage{
\S4method{pileup_peak}{GRanges}(object, bamf = NULL, d = NULL)
}

\arguments{
  \item{object}{\linkS4class{GRanges} object containing the genomic coordinates of the peaks.
}
    \item{bamf}{Path to the .bam file
    used to compute the coverage function.
    The associated .bam.bai index file must also be present.}
    
    \item{d}{integer. Total length of the fragments. Positive and negative reads are extended in
    their 3' direction. Default is NULL; this value can be estimated by  \link[=compute_fragments_length]{compute_fragments_length}.}
}

\value{
the \linkS4class{GRanges} \code{object} with the new metadata column \code{counts} 
containing the coverage functions of the peaks.
}

\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
# GRanges object

data(GR)

# import the .bam file

bamf <- system.file("extdata", "test.bam", package="FunChIP",
                      mustWork=TRUE)

# extract the first 10 peaks of the GRange
# and compute the corresponding read counts
# with fragment length 160.

peaks <- pileup_peak(GR[1:10], bamf, d = 160)

}
