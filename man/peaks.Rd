\name{peaks}
\alias{peaks}
\alias{peaks.data}
\alias{peaks.data.classified}
\alias{peaks.data.cluster}
\alias{peaks.data.smooth}
\alias{peaks.data.summit}
\docType{data}
\title{
Coverage of a ChIP-Seq experiment
}
\description{
It contains a data-set used in all the examples of the \link{FunChIP}, together with all the metadata 
columns generated in the intermediate steps.
}
\usage{data("peaks")}
\format{
Several data-set are included
\itemize{
    \item{peaks.data: a \linkS4class{GRanges} object with 10 ranges and 1 metadata
     column: 
     \itemize{
        \item{\code{counts}. A vector for each range, with length equal to the
       width of the range, containing the coverage of the range, i.e. the base-level read counts. 
       It can be computed with the \link[=pileup_peak-method]{pileup_peak} method.}
        }
    }
    \item{peaks.data.smooth:  a \linkS4class{GRanges} object with 10 ranges and 6 metadata
     columns: 
     \itemize{
        \item{\code{counts}. As in \code{peaks.data}}.
        \item{\code{spline}. A vector for each range, containing the evaluation of the 
        spline approximation of the peak for each genomic base. It can be computed with the 
        \link[=smooth_peak-method]{smooth_peak} method. }
        \item{\code{spline_der}. A vector for each range, containing the evaluation of the derivatives of the
        spline approximation of the peak for each genomic base. It can be computed with the 
        \link[=smooth_peak-method]{smooth_peak} method. }
        \item{\code{width_spline}. Integer. The number of evaluated points of the 
        spline approximation, i.e. the number of non-zero points, for each range. It can be
        computed with the \link[=smooth_peak-method]{smooth_peak} method. }
        \item{\code{start_spline}. Integer. The starting point of the
        spline approximation. It could be smaller than \code{start(object)} since the
        approximation can increase the values at the border to make the curve smooth. It can be 
        computed with the \link[=smooth_peak-method]{smooth_peak} method. }
        \item{\code{end_spline}. Integer. The end point of the
        spline approximation. It could be larger than \code{end(object)} since the
        approximation can increase the values at the border to make the curve smooth. It can be 
        computed with the \link[=smooth_peak-method]{smooth_peak} method. }
     }
    }
    \item{peaks.data.summit:  a \linkS4class{GRanges} object with 10 ranges and 7 metadata
     columns: 
     \itemize{
        \item{\code{counts}. As \code{peaks.data}.}
        \item{\code{spline}. As \code{peaks.data.smooth}. }
        \item{\code{spline_der}. As \code{peaks.data.smooth}. }
        \item{\code{width_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{start_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{end_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{summit_spline}. The distance from the starting point of the spline
        of the maximum point (integer) of
        the spline (or the summit of the peak), for each range. It can be computed with the 
        \link[=summit_peak-method]{summit_peak} method}
     }
    }
     \item{peaks.data.cluster:  a \linkS4class{GRanges} object with 10 ranges and 12 metadata
     columns:
     \itemize{
        \item{\code{counts}. As \code{peaks.data}.}
        \item{\code{spline}. As \code{peaks.data.smooth}. }
        \item{\code{spline_der}. As \code{peaks.data.smooth}. }
        \item{\code{width_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{start_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{end_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{summit_spline}. As \code{peaks.data.summit}.}
        \item{\code{cluster_NOshift}. A vector of length 5 for each range,
        containing the label of the cluster assigned to the peak in case of clustering without alignment. 
        For example, the second element of the vector is the label of the corresponding peak 
        when the k-mean alignment algorithm is run with 2 clusters.  
        It can be computed with the \link[=cluster_peak-method]{cluster_peak} with \code{n.clust = 1:5}
        and \code{shift.peak=FALSE}.}
        \item{\code{dist_NOshift}. A vector of length 5 for each range,
        containing the distance from the center of the cluster assigned to the peak in case of  
        clustering without alignment. For example, the second element of the vector 
        is the distance of the corresponding peak from center of the corresponding cluster when the k-mean 
        alignment algorithm is run with 2 clusters. 
        It can be computed with the \link[=cluster_peak-method]{cluster_peak} with \code{n.clust = 1:5}
        and \code{shift.peak=FALSE}.}
        \item{\code{cluster_shift}. A vector of length 5 for each range,
        containing the label of the cluster assigned to the peak in case of clustering with alignment. 
        For example, the second element of the vector is the label of the corresponding peak 
        when the k-mean alignment algorithm is run with 2 clusters. 
        It can be computed with the \link[=cluster_peak-method]{cluster_peak} with \code{n.clust = 1:5}
        and \code{shift.peak=TRUE}.}
        \item{\code{coef_shift}. A vector of length 5 for each range,
        containing the optimal shift coefficient of the peak. 
        For example, the second element of the vector is the shift coefficient of the corresponding
        peak when the k-mean alignment algorithm is run with 2 clusters. 
        It can be computed with the \link[=cluster_peak-method]{cluster_peak} with \code{n.clust = 1:5}
        and \code{shift.peak=TRUE}.}
        \item{\code{dist_shift}. A vector of length 5 for each range,
        containing the distance from the center of the cluster assigned to the peak in case of 
        clustering with alignment. For example, the second element of the vector 
        is the distance of the corresponding peak from the center of the corresponding cluster when the k-mean 
        alignment algorithm is run with 2 clusters. 
        It can be computed with the \link[=cluster_peak-method]{cluster_peak} with \code{n.clust = 1:5}
        and \code{shift.peak=TRUE}.}
     }
    }
    \item{peaks.data.classified:  a \linkS4class{GRanges} object with 10 ranges and 13 metadata
     columns:
     \itemize{
        \item{\code{counts}. As \code{peaks.data}.}
        \item{\code{spline}. As \code{peaks.data.smooth}. }
        \item{\code{spline_der}. As \code{peaks.data.smooth}. }
        \item{\code{width_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{start_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{end_spline}. As \code{peaks.data.smooth}.  }
        \item{\code{summit_spline}. As \code{peaks.data.summit}.}
        \item{\code{cluster_NOshift}. As \code{peaks.data.cluster}. }
        \item{\code{dist_NOshift}. As \code{peaks.data.cluster}. }
        \item{\code{cluster_shift}.  As \code{peaks.data.cluster}.}
        \item{\code{coef_shift}.  As \code{peaks.data.cluster}..}
        \item{\code{dist_shift}.  As \code{peaks.data.cluster}.}
        \item{\code{cluster}. Integer. The index of the final label assigned, for each range.
        It can be computed with the \link[=choose_k-method]{choose_k} 
        with \code{k = 3} and \code{shift.peak = TRUE}.}
     }
    }
    

}
}

\examples{
data(peaks)
}
\keyword{datasets}
