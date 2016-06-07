\name{cluster_peak}
\docType{methods}
\alias{cluster_peak-method}
\alias{cluster_peak}
\alias{cluster_peak,GRanges-method}
\title{
Clustering the peaks with the k-mean alignment algorithm
}
\description{
It classifies and aligns the peaks
stored in the \linkS4class{GRanges} object. The method applies the k-mean alignment algorithm
with shift of the peaks and distance based on the convex combination of the
\eqn{L^p} distances between the spline-smoothed peaks and their derivatives. The order \eqn{p}
can be one of \eqn{1}, \eqn{2} and \eqn{\infty}.

}
\usage{
\S4method{cluster_peak}{GRanges}(object, parallel = FALSE, num.cores = NULL,
    n.clust = NULL,  seeds = NULL, shift.peak = NULL, weight = NULL,
    subsample.weight = 100, alpha = 1, p = 1, t.max = 0.5,
    plot.graph.k = TRUE, verbose = TRUE )
}
\arguments{
  \item{object}{ \linkS4class{GRanges} object of length \eqn{N}. It must contain the metadata columns
  \code{spline}, \code{spline_der}, \code{width_spline}, computed by \link[=smooth_peak-method]{smooth_peak}.
  }

  \item{parallel}{logical. If \code{TRUE}, the clustering for different values of the 
  parameter \eqn{k} in \code{n.clust} are run in parallel. Default is \code{FALSE}.
  }

  \item{num.cores}{integer. If \code{parallel} is \code{TRUE}, it indicates the number of
  cores used in the parallelization. If \code{NULL} (default), the number of cores is
  automatically identified.
  }

  \item{n.clust}{integer vector (or scalar). Number of clusters in which the data set 
  is divided (possibly one, if \code{n.clust} is a scalar). For each value of the vector,
  the cpp function \code{kmean_function} is called. 
  }

  \item{seeds}{vector. Indices of the initial centers of the clusters, needed to initialize the k-mean procedure. The k-mean alignment, like all the k-mean-like algorithms, 
  is dependent on the choice of the initial centers 
  of the clusters, and each initialization
  of the seeds can generate slightly different results. The 
  values must be included in \eqn{1, \ldots, N}. The length of the vector must be equal to
  the maximum number of clusters analyzed (\code{max(n.clust)}), otherwise
  it is truncated to this value, or the missing values are randomly generated.
  If \code{NULL} (default), the seeds are detected as the most central values (i.e. peaks with minimum distance from the others) of
  the set of peaks. If \code{seeds='random'}, the
  centers are randomly generated.  
  }

  \item{shift.peak}{logical. It indicates whether the alignment via a translation of the abscissae
  is performed (\code{shift.peak = TRUE}) or not (\code{shift.peak = FALSE}). If no value is
  provided (\code{shift.peak = NULL}, default), both analyses are performed.
  }

  \item{weight}{real. Weight \eqn{w} of the distance function (see Details for the 
  definitions of the distance function), needed to
  make the distance between splines and  derivatives comparable.
  If no value is provided (default is \code{NULL}), it is computed as the median of the
  ratio between the pairwise distances of the data (\eqn{d_0 (i,j)}) and of the derivatives (\eqn{d_1(i,j)})
  \deqn{
      w= \textrm{median} \frac{d_0(i,j)}{d_1(i,j)}
      }{ w = median d0(i,j)/d1(i,j)}
  with \eqn{i, j = 1: \ldots N}. 
  }

  \item{subsample.weight}{integer value. Number of data points used 
  to define the \code{weight}, if not assigned.
  Using all the peaks to define the weight can be computationally expensive and therefore
  a subsampling is suggested. If \code{subsample.weight=NULL} all the data will be used. 
  Default is 100, which is a reasonable trade off between running time 
	and reliability of the estimation.
  }

  \item{alpha}{real value between 0 and 1. Value of the convex weight \eqn{\alpha} of the distance to
  balance the distance between data and derivatives.
  See details for the definition. Default is 1.
  }

  \item{p}{integer value in \{0, 1 , 2\}. Order of the \eqn{L^p} distance
  used. In particular \code{p = 0} stands for the \eqn{L^{\infty}} distance,
  \code{p = 1} for \eqn{L^1} and \code{p = 2} for \eqn{L^2}}

  \item{t.max}{real value. It tunes the maximum shift allowed. In
  particular the maximum shift at each iteration is computed as
  \deqn{
    \textrm{max\_shift} = \textrm{\code{t.max}} * \textrm{range}({\textrm{\code{object}}})
  }{max_shift = t.max * range(object)}
  and the optimum registration coefficient will be identified between - max_shift and
  + max_shift. range({\code{object}}) is the maximum amplitude of the peaks. Default is 0.5.
  }

  \item{plot.graph.k}{logical. If \code{TRUE} the graph of the average distance between the
  data and corresponding center of the cluster, varying the number of clusters is plotted. 
  If \code{align=NULL}, both the analysis with and without alignment
  are performed, two lines are drown to show the decrease of the
  global distance introduced by the alignment procedure. Default is \code{TRUE}
  }

  \item{verbose}{logical. If \code{TRUE}, some parameters of the
  algorithm and the progress of the iterations are shown, if \code{FALSE} no
  information is provided. Default is \code{TRUE}, but consider to set 
  the parameter to \code{FALSE} in case of  parallel runs, to avoid the overlap of their outputs.
  }
}

\details{
See [Sangalli et al., 2010] and the package vignette for the complete 
description of the algorithm.
The algorithm is completely defined once we fix the family of the
warping function for the alignment and the distance function.
In this function we focus only on the specific case of
\itemize{
  \item{warping functions: shifts with integer coefficients
  \deqn{
  h(t) = t + c,
    }
  with \eqn{c} an integer value;
  }
  \item{distance: convex combination of the \eqn{L^p}
  distance between data and derivatives. The distance
  between \eqn{f} and \eqn{g} is
  \deqn{
  d (f, g) = (1 - \alpha) \| f - g \|_{p}
   + \alpha \,\, w \| f^{\prime}
   - g^{\prime} \|_{p}
    }{
  d(f, g) = (1 - \alpha) || f - g ||_p
   + \alpha   w || f'
   - g' ||_p
    }
  } The choice of \eqn{\| \cdot \|_{p}}{|| . ||_p} corresponds to
  the value of \code{p} in input. In particular \code{p = 0} stands for
  \eqn{\| \cdot \|_{L^{\infty}}}{||.||_L^\infty}, \code{p = 1} for \eqn{\| \cdot \|_{L^{1}}}{|| . ||_L^1}
  and \code{p = 2} for \eqn{\| \cdot \|_{L^{2}}}{|| . ||_L^2}
}

}

\value{
the \linkS4class{GRanges} \code{object}  with new metadata columns:
\itemize{
  \item{if \code{align} is \code{TRUE} or \code{NULL}, i.e. the
  clustering with alignment is performed the following metadata columns are added:
  \itemize{
      \item{\code{cluster_shift}: for each peak, a vector of length equal to the 
      maximum number of chosen clusters,
      containing at each position \eqn{k} the label of the cluster the peak is
      assigned to, when the total number of clusters is \eqn{k}
      and alignment is performed during the clustering.
      If \eqn{k} is not present in the \code{n.clust} vector, the corresponding
      value is \code{NA}. }
      \item{\code{coef_shift}: for each peak, a vector of length equal to the 
      maximum number of chosen clusters,
      containing at each position \eqn{k} the shift coefficient
      assigned to the peak, when the total number of clusters is \eqn{k}
      and alignment is performed during clustering.
      If \eqn{k} is not present in the vector \code{n.clust} the corresponding
      value is \code{NA}.}
      \item{\code{dist_shift}: for each peak, a vector of length equal to the 
      maximum number of chosen clusters,
      containing at each position \eqn{k} the distance of the specific peak
      from the corresponding center of the cluster, when
      the total number of clusters is \eqn{k}
      and alignment is performed during clustering.
      If \eqn{k} is not present in the vector \code{n.clust} the corresponding
      value is \code{NA}.}
  }}
   \item{if \code{shift.peak} is \code{FALSE} or \code{NULL}, i.e. clustering is performed
   without alignment, the following metadata columns are added:
  \itemize{
      \item{\code{cluster_NOshift}: for each peak, a vector of length equal to the 
      maximum number of chosen clusters, 
      containing at each position \eqn{k} the label of the cluster
      the peak is assigned to, when the total number of clusters is \eqn{k}
      and no alignment is performed during clustering.
      If \eqn{k} is not present in the vector \code{n.clust} the corresponding
      value is \code{NA}.}
      \item{\code{dist_NOshift}: for each peak, vector of length equal to the 
      maximum number of chosen cluster,
      containing at each position \eqn{k} the distance of the peak
      from the corresponding center of the cluster , when
      the total number of clusters is \eqn{k}
      and no alignment is performed during clustering.
      If \eqn{k} is not present in the vector \code{n.clust} the corresponding
      value is \code{NA}.}
  }}
}
}
\references{
Sangalli, L. M., Secchi, P., Vantini, S. and Vitelli, V., 2010.
  K-mean alignment for curve clustering. Computational Statistics and Data Analysis, 54 1219 - 1233.
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\seealso{
\link[=choose_k-method]{choose_k}
}

\examples{

# load the data
data(peaks)

# cluster and align the data as a 
# function of the
# number of cluster (from 1 to 5)
# with and without alignment.
# The automathically generated plot 
# can be used to detect the
# optimal number of clusters and the
# classification method to be used
# (with or without alignment)

clustered_peaks <- cluster_peak ( peaks.data.summit, parallel = FALSE ,
                                  n.clust = 1:5, shift.peak = NULL,
                                  weight = 1, alpha = 1, p = 2,
                                  plot.graph.k = TRUE, verbose = TRUE )


}
