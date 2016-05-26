\name{summit_peak}
\docType{methods}
\alias{summit_peak-method}
\alias{summit_peak}
\alias{summit_peak,GRanges-method}
\title{
Finding the summits of the peaks of the \linkS4class{GRanges} object}
\description{
It identifies the summit of the peak and stores it in a new metadata column.
}
\usage{
\S4method{summit_peak}{GRanges}(object, summit = NULL)
}
\arguments{
  \item{object}{\linkS4class{GRanges} object  of length \eqn{N}. If \code{summit} is not provided, 
  \code{object} must contain the metadata column \code{spline}.
}
  \item{summit}{vector of length \eqn{N}. It contains the x coordinate of the summit of the peaks, i.e.
  the distance of the summit from the starting position of the spline approximation of peak (distance
  from \code{start_spline}). 
  If \code{summit} is \code{NULL} the
  summit of each peak is identified as the maximum point of the \code{spline}.}
}
\value{
the \linkS4class{GRanges} object with the new metadata column \code{summit_spline}. 
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
data(peaks)

# Computing the summits of the peaks from
# the spline-smoothed approximation.

peaks.spline.summit <- summit_peak(peaks.data.smooth)

}
