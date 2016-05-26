\name{distance_peak}
\alias{distance_peak}
\title{
Computing distance matrices among spline-smoothed peaks and among their derivatives.
}
\description{
Given a metric (\eqn{L^1}, \eqn{L^2} or \eqn{L^\infty}), it computes the
pairwise distance among the spline approximation of the peaks and  among their
derivatives.
}
\usage{
distance_peak(object, p = 1)
}

\arguments{
  \item{object}{ \linkS4class{GRanges} object of length \eqn{N}. It has to contain the metadata columns
  \code{spline}, \code{spline_der}, \code{width_spline}.
  }

  \item{p}{integer. It must assume values in \code{{0, 1, 2}}. They correspond
  respectively to the \eqn{L^{\infty}}, \eqn{L^1} and \eqn{L^2} norm. See details
  for the definition of the norms.
  }

}
\details{
This function computes the pairwise distance of a set of \eqn{N} peaks. 
Given the spline-smoothed peaks \eqn{s_i} and their derivatives \eqn{s^{\prime}_i}
(\eqn{i = 1, \ldots, N}),
it returns two matrices \code{dist_matrix_d0} and \code{dist_matrix_d1} whose
elements \eqn{(i,j)}, with \eqn{i, j = 1, \ldots,  N}, are
\deqn{
  \textrm{\code{dist\_matrix\_d0}}  (i,j) =  \| s_i - s_j \|_{p} \,\,\,\,    \textrm{\code{dist\_matrix\_d1}} (i,j) =  \| s^{\prime}_i - s^{\prime}_j \|_{p}
  }{dist_matrix_d0 (i.j) = || s_i - s_j||_p    dist_matrix_d1 (i.j) = || s'_i - s'_j||_p  }.

In particular, in order to define the distance between two functions \eqn{f} and \eqn{g}:
\itemize{

  \item{define a common domain \eqn{U}, given by the union of the domains of \eqn{f}
  and \eqn{g}. If a function is not defined on the whole domain, it is
  extended with 0's on the missing parts. The value of 0 is chosen because 
  the background of the peaks has been removed during the
  definition of the splines, and
  hence they can now be continuously extended with 0's.}
  \item{choose the order of the norm \code{p}. Given the function \eqn{f} defined on \eqn{U}
  we have:
  \deqn{\| f \|_0 = \| f \|_{L^{\infty}} = \textrm{max}_{x \in U} | f(x) | }{|| f ||_0 = || f ||_{L^\infty} = max_{x in U} | f(x) | },

  \deqn{\| f \|_1 = \| f \|_{L^{1}} = \int_U |f(x)| dx  }{|| f ||_1 = || f ||_{L^{1}} = int (|f(x)| dx) },

  \deqn{\| f \|_2 = \| f \|_{L^{2}} = \sqrt{\int_U (f(x))^2 dx}}{|| f ||_2 = || f ||_{L^{2}} = sqrt(int (f(x))^2 dx)}.
      }

  }
}
\value{
list with two components
 \item{dist_matrix_d0 }{ \eqn{N} x \eqn{N} matrix of the pairwise distances
 between the splines.}

 \item{dist_matrix_d2 }{ \eqn{N} x \eqn{N} matrix of the pairwise distances
 between the derivatives of splines.}

}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
data(peaks)

# Compute the pairwise
# L2 distance between the peaks
dist_matrices <- distance_peak(peaks.data.summit, p = 2)

}
