\name{smooth_peak}
\docType{methods}
\alias{smooth_peak-method}
\alias{smooth_peak}
\alias{smooth_peak,GRanges-method}
\title{
Spline smoothing of the peak}
\description{
It approximates the read counts associated to every peak with a suitable B-spline function, 
so that a smoothing representation of the peaks is obtained.
The first derivative of the spline is also computed.  To obtain a smooth representation, the peak is extended and  new initial and final points are identified. See the Vignette of the \link{FunChIP} package for a graphical representation of the spline approximation.
}

\usage{
\S4method{smooth_peak}{GRanges}(object, n.breaks = 100, subsample = TRUE, 
    subsample.data = 100, order = 4,  
    lambda = (10^(seq(-5,5, by = 0.5))),
    GCV.derivatives = TRUE , plot.GCV = FALSE)
}
\arguments{
  \item{object}{ \linkS4class{GRanges} object. It must contain the metadata column \code{counts}.
  }

  \item{n.breaks}{ integer. Number of breaks, or knots, for the B-spline basis domain definition. 
  Default is \code{100}.
  }

  \item{subsample}{logical. If \code{TRUE}, only a random subset (of size fixed by
  the parameter \code{subsample.data}) is used to identify the optimal value of lambda
  for the penalization via cross-validation. If \code{subsample=FALSE}, all the peaks of the
  \linkS4class{GRanges} data will be used. To contain running times,
  it is suggested to maintain the default value \code{subsample = TRUE}.
  }
  
  \item{subsample.data}{integer. Number of data used for the 
  cross-validation (if \code{subsample.data} is \code{TRUE}).
  Default value is 100. If \code{subsample = FALSE}, all data points will be used and
  \code{subsample.data} is ignored.
  }

  \item{order}{integer. Order of the B-spline basis used for the smoothing. The order is
  one higher than the degree of the spline. Default is 4 (cubic splines).
  }

  \item{lambda}{vector (or single value). Contains all the possible values of the smoothing
  parameter to be considered for the final choice. If a single value is provided, this will
  be automatically chosen for the smoothing. Default value is \eqn{10^{\textrm{\code{seq(-5,5,by=0.5)}}}}
  to analyze a sufficiently wide set of values. See details below. 
  }

  \item{GCV.derivatives}{logical. If \code{TRUE} the Generalized Cross
  Validation index (GCV) on the derivatives is considered
  as criteria to identify \eqn{\lambda}, 
  otherwise the GCV is computed on the data.
  Default is \code{TRUE}.
  }

  \item{plot.GCV}{logical. If \code{TRUE}, the plot of the GCV of the data and derivatives
  is shown as a function of \eqn{\lambda}. Default value is \code{FALSE}.}
}
\details{

 
It creates a piece-wise polynomial of fixed \code{order} \eqn{s}
approximating the data (B-spline expansion, Ramsay and Silverman, 2005). 
Given the point wise defined function \eqn{f: (x,f(x))}, the
\code{smooth_peak} method returns
the evaluation of \eqn{s} on the \eqn{x}
grid (\eqn{s(x)}) minimizing, for a fixed \eqn{\lambda},

\deqn{
ERR(\lambda) = \| f - s \|_{L^2}^2 + \lambda \|s''\|_{L_2}^2
}{ERR(\lambda) = || f - s ||^2_{L^2} + \lambda ||s''||^2_{L_2}},
with \eqn{s''} being the second derivative of the function \eqn{s} and \eqn{\| s \|_{L^2}}{||s||_{L^2}} 
the \eqn{L^2} norm of the function, i.e. the integral on the domain of \eqn{s} of \eqn{s^2}.

The choice of \eqn{\lambda} is crucial for the definition of the
spline, and it can be selected by 
minimizing the Generalized Cross-Validation index

\deqn{
GCV(\lambda) = \frac{n SSE}{(n-df(\lambda))^2}
}{
GCV(\lambda) = (n SSE)/(n-df(\lambda))^2
},
with SSE the error computed as
\deqn{SSE = \| f - s \|_{L^2}^2}{SSE = || f - s ||^2_{L^2}},
if \code{GCV.derivatives = FALSE}, or
\deqn{SSE = \| \nabla f - s' \|_{L^2}^2}{SSE = ||  grad(f) - s' ||^2_{L^2}},
if \code{GCV.derivatives = TRUE}, and \eqn{df(\lambda)} is the
number of the degrees of freedom of the basis expansion automatically
computed from \eqn{s}. For further details on the cross-validation procedure
and on the computation of the number of degrees of freedom see Ramsay and Silverman, 2005.

If \code{plot.GCV} is \code{TRUE}, the plot of the GCV index as a function of \eqn{\lambda}
is presented, which can be used to identify the optimal value of the parameter. If the plot
is decreasing in \eqn{\lambda}, one could consider to increase the allowed values of \eqn{\lambda}
to find the  minimum of the curve.

}
\value{
the \linkS4class{GRanges} \code{object} with new metadata columns:
\itemize{
\item{\code{width_spline} integer. Value containing the width of the smoothed peak, 
i.e. the number of non-zero values of the spline approximation.
This value is not necessarily equal to the original width of the peak, as the approximation
can stretch outside the original width of the peak: to ensure smoothness some 0 values can be
introduced at the edges of the region.}
\item{\code{spline} vector. Evaluation of the spline on the grid of size \code{width_spline}.}
\item{\code{spline_der} vector. Evaluation of the derivatives
of the spline on the grid of size \code{width_spline}.}
\item{\code{start_spline} integer. Genomic coordinate of the initial point of the
spline approximation.}
\item{\code{end_spline} integer. Genomic coordinate of the final point of the
spline approximation.}
}
}

\references{
Ramsay, J.O., Silverman, B.W., 2005. Functional Data Analysis, 2nd ed. Springer, New York, NY.
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
data(peaks)

# it computes the spline approximation
# of the pealks given the
# GRange with the metadata counts. 
# It is obtained by the pileup_peak method

# Default paramters are used: GCV is
# computed on the derivatives.

peaks.spline <- smooth_peak(peaks.data, lambda = 10^(-4:6), 
                            subsample.data = 50, GCV.derivatives = TRUE )

}
