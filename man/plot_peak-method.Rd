\name{plot_peak}
\docType{methods}
\alias{plot_peak-method}
\alias{plot_peak}
\alias{plot_peak,GRanges-method}

\title{
Plotting the peaks.}
\description{
It plots the peaks both as counts and splines (if provided).
Peaks are centered
around their summit. If the clustering has been performed, peaks can be aligned and 
divided in cluster, shown in different panels.
}
\usage{
\S4method{plot_peak}{GRanges}(object, index = 1:length(object),
line.plot = "spline", col = NULL, 
shift = NULL, k = NULL, cluster.peak = FALSE )
}

\arguments{
  \item{object}{\linkS4class{GRanges} object of length \eqn{N}. It must contain the metadata column \code{counts}.
  If it contains the metadata column \code{summit_spline}, peaks are plotted centered around their summits
  (i.e. the 0 of the abscissa is fixed as the summit of the peak). If summit is not
  provided peaks are not centered.
  }

  \item{index}{vector. Indices of the peaks to be plotted. Default is \code{1:length(object)}
  to plot all the peaks.}

  \item{line.plot}{string. Type of plot. If \code{'spline'} (default), the spline approximation of the 
  peaks is plotted. If \code{'counts'}, only the raw data are plotted. If \code{'both'}, both the 
  raw data and the approximation spline are plotted.}
  
  \item{col}{vector. Colors used to plot the peaks. If \code{NULL}, the
  rainbow color palette is used. If it has a single value, all the peaks are
  plotted with the same color. If it is a vector shorter than \eqn{N},
  only the first element of \code{col} is used. Default is \code{NULL}.}

  \item{shift}{logical. This parameter controls the abscissae of the
  plotted peaks, and plays two different roles, depending on the \code{cluster.peak} 
  parameter. \\
  If \code{cluster.peak} is \code{FALSE}, \code{shift = TRUE} means that
  peaks are plotted aligned 
  around the \code{summit_spline} point. The \linkS4class{GRanges} object must
  contain the results of the \link[=summit_peak-method]{summit_peak} method.
  If \code{cluster.peak} is \code{FALSE} and \code{shift = FALSE},
  peaks are plotted with no centering around the summit. If \code{cluster.peak} is \code{FALSE} and \code{shift = NULL} (default), peaks are
  centered, if the metadata column \code{summit_spline} is present in \code{object}, 
  otherwise they are plotted with the original abscissae, with no centering around the summit. \\
  If \code{cluster.peak = TRUE}, 
  the parameter \code{shift} sets the clustering result to be plotted. 
  If \code{TRUE} the \code{shift} results are plotted, otherwise the
  results associated to \code{NOshift} clustering are presented.}
  
  \item{k}{integer. If \code{cluster.peak = TRUE}, 
  results corresponding to \code{k} number of clusters are plotted.
  It must be set to a value included in the \code{n.clust} parameter 
  of the correspondent \code{object}.
  }
  
  \item{cluster.peak}{logical. It \code{FALSE}, \code{object} contains the set of
  peaks to be plotted without classification (centered or not around the
  summit). If \code{TRUE}, the \code{plot_peak} method plots the result
  of the classification associated to \code{shift} and \code{k} parameters;
  in this case, \code{object} must be 
  the output of the \link[=cluster_peak-method]{cluster_peak} 
  and must contain the correspondent classification. Default is \code{FALSE}.} 
}
\value{Graphical method to graphically represent data. No output returned,}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}
\examples{
# load the data
data(peaks)

# First example:
# plot of the spline approximation
# of the first 10 peaks
# centered around their summit

plot_peak(peaks.data.summit, index = 1:10,
          shift = TRUE)

# Second example:
# plot of the peaks 
# divided in the k=3 
# and shift = TRUE clusters
# obtained with the cluster_peak method
plot_peak(peaks.data.cluster,
shift = TRUE, k = 3, cluster.peak = TRUE)

}
