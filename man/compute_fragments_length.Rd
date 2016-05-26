\name{compute_fragments_length}
\alias{compute_fragments_length}
\title{
Computing the length of the fragments in the .bam file.
}
\description{
Given a .bam file and \linkS4class{GRanges} object, it computes
the positive and negative coverage for each  \linkS4class{GRanges} element, estimates the distance between positive and negative peaks, and finally the fragment length \code{d}, i.e. the sum of the length of the reads and the 
distance between positive and negative peaks. See Details and the package vignette for the description of the method.
}
\usage{
compute_fragments_length(object, bamf, min.d = 0, max.d = 200)
}

\arguments{
  \item{object}{ \linkS4class{GRanges} object of length \eqn{N}.
  }

    \item{bamf}{Path to the .bam file
    used to compute the coverage function.
    The associated .bam.bai index file must also be present.}
    
    \item{min.d}{integer. Minimum value for the distance between positive and 
    negative peaks. Default is 0.}
    
    
    \item{max.d}{integer. Maximum value for the distance between positive and 
    negative peaks. Default is 200.}

}
\details{
Given a set of \eqn{n= 1, \ldots N} regions, characterized by their positive and negative coverages, the function computes the distance between the positive peak \eqn{f_{n+}} and the negative peak, shifted by \eqn{\delta} \eqn{f_{n-}^{\delta}}:
\deqn{
D(f_{n+}, f_{n-}^{\delta}) = \frac{\| f_{n+} - f_{n-}^{\delta}\|^2_{L^2}}{\textrm{width}(\textrm{union}(f_{n+}, f_{n-}^{\delta}))}
}
The function computes the \eqn{d_{pn}} minimizing the distance between postive and negative peaks
\deqn{
d_{pn} = \textrm{argmin}_{\delta \in [\code{min.d}, \code{max.d}]} \sum_{n=1 }^N D(f_{n+}, f_{n-}^{\delta})
}
The function returns both the plot of the global distance vs the fragment length \eqn{d = d_{pn} + r}, where \eqn{r} is the length of the reads, and the optimum value for \eqn{d}. \eqn{r} is also estimated from the .bam file as the average of the read lengths.
}
\value{
optimum value of the parameter \eqn{d}, to be used in the \link[=pileup_peak-method]{pileup_peak} method.
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
# GRanges object

data(GR100)

# import the .bam file

bamf <- system.file("extdata", "test.bam", package="FunChIP",
                      mustWork=TRUE)
                      
# compute the estimated fragment length

d <- compute_fragments_length(GR[1:10], bamf, min.d = 0, max.d = 200)                      
}
