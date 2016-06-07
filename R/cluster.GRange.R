
setGeneric("cluster_peak", function(object, ...) setGeneric("cluster_peak"))

cluster.GRange <- function(object, parallel = FALSE, num.cores = NULL,
                           n.clust = NULL,  seeds = NULL, shift.peak = NULL, weight = NULL, subsample.weight = 100,
                           alpha = 1, p = 1, t.max = 0.5,  plot.graph.k = TRUE, verbose = TRUE)
{
  # This function will call the .cpp funciton
  #   SEXP kmean_function (SEXP x, SEXP spline, SEXP spline_der,
  #                        SEXP lenght, SEXP seeds, SEXP align,
  #                        SEXP k, SEXP weight_input, SEXP alpha_input,
  #                        SEXP p_input, SEXP t_max_input, SEXP verbose)


  if (is.null(object$counts)) 
  {
        stop('No information on the peaks provided!')
  }
  
  if (is.null(n.clust))
  {
      stop('No information on the number of clusters provided')
  }
    
  if (is.null(object$spline) || is.null(object$spline_der) || is.null(object$width_spline) )
  {
    stop('spline approximation not provided. Run smooth_peak before the clustering')
  }

  if (is.null(object$summit_spline))
  {
      stop('summit is not a metadata of the GRange. Run summit_peak before the clustering')
  }
  
  if (max(n.clust)>=length(object))
  {
      stop ('Maximum number of cluster greater then the number of observations')
  }
  if (length(object)==1) 
  {
        x_centered_list <- vector('list', 1)
        x_centered_list[[1]] <- (-object$summit_spline +1) : (object$width_spline - object$summit_spline)
  }else
  {
        x_centered_list <- mapply(function(x,y){(-x+1):(y-x)}, object$summit_spline, object$width_spline)
  }
    
  x_matrix <- unlist.counts(x_centered_list, object$width_spline)
  spline_matrix <- unlist.counts(object$spline, object$width_spline)
  spline_der_matrix <- unlist.counts(object$spline_der, object$width_spline)

  # check the dimentions of x, spline, spline_der

  n.data <- dim(x_matrix)[1]
  n.points <- dim(x_matrix)[2]

  if ((dim(spline_matrix)[1] != n.data)  | (dim(spline_matrix)[2] != n.points))
  {
    stop('dimentions of x and spline don\'t coincide')
  }

  if ((dim(spline_der_matrix)[1] != n.data)  | (dim(spline_der_matrix)[2] != n.points))
  {
    stop('dimentions of x and spline_der don\'t coincide')
  }

  
  # check on the definition of the parameters of the distance
  
  if ( (p != 0) && (p != 1) && (p != 2) )
  {
      stop ('invalid value for p, It must be 0,1,2 ')
  }
  
  if ( (alpha < 0) || (alpha > 1) )
  {
      stop ('invalid value for alpha. It must be included in [0, 1]')
  }
  
  if (is.null(weight))
  {
      if ( (is.null(subsample.weight)) || (subsample.weight >= length(object)))
      {
          dist_matrix <- distance_peak(object, p)
          upper_triang_d0 <- dist_matrix$dist_matrix_d0[upper.tri(dist_matrix$dist_matrix_d0)]
          upper_triang_d1 <- dist_matrix$dist_matrix_d1[upper.tri(dist_matrix$dist_matrix_d1)]
          weight <- median(upper_triang_d0/upper_triang_d1)
      }else
      {
          # subsample the data to define the weights. Using all the
          # data can be computationally expensive
          
          object_here <- object[sort(sample(1:length(object), subsample.weight))]
          dist_matrix <- distance_peak(object_here, p)
          upper_triang_d0 <- dist_matrix$dist_matrix_d0[upper.tri(dist_matrix$dist_matrix_d0)]
          upper_triang_d1 <- dist_matrix$dist_matrix_d1[upper.tri(dist_matrix$dist_matrix_d1)]
          weight <- median(upper_triang_d0/upper_triang_d1)
          
      }
      
  }
  
  
  
  if (is.null(seeds))
  {
      dist_matrix <- distance_peak(object, p)
      dist_matrix_global <- (1 -alpha) * dist_matrix$dist_matrix_d0 + alpha * weight * dist_matrix$dist_matrix_d1
      dist_elem <- colSums(dist_matrix_global)
      seeds <- t(as.vector(sort(dist_elem, index.return = TRUE, decreasing = FALSE)$ix[1:max(n.clust)]))
      
  }
  
  if (!is.numeric(seeds))
  {
      if (seeds == 'random')
      {
          seeds <- sort(sample(1:n.data, max(n.clust)))
      }else
      {
          stop('invalid value for seeds')
      }
  }

  #check the uniqness of the seeds
  seeds <- unique(seeds)

  if ( ( min(seeds) < 0 ) | ( max(seeds) > n.data ) )
  {
    stop('invalid values in seed. Seeds are the indices of the data
         choosen as starting centers of the clusters and then must be inclued in 0, n.data')
  }

  if (length(seeds) < max(n.clust))
  {
    warning ('not enough values of seeds provided, random choice of the missing values')
    seeds_new <- sample((1:n.data)[-seeds], max(n.clust)-length(seeds))
    seeds <- c(seeds, seeds_new)
  }

  if (length(seeds) > max(n.clust))
  {
    warning ('too many seeds provided. Only the first used')
    seeds <- seeds[1:max(n.clust)]
  }


  if (is.null(shift.peak))
  {
    warning ('no value for shift.peak provided. Both results with no registration and
             results with registration provided')
  }else
  {
    if (shift.peak != TRUE & shift.peak != FALSE )
    {
      stop ('invalid value for shift.peak It must be TRUE or FALSE')
    }
  }


  if (verbose)
  {
    verb <- 'T'
    if (parallel)
    {
        warning('both parallel and verbose set to TRUE, the output will be unclear. ')
    }
  }else
  {
    verb <- 'F'
  }

  # to define the variable k as global.
  k = NULL
  
  if (parallel==TRUE)
  {
      
    if (is.null(num.cores))
    {
      warning ('the number of cores will be automatically identified.')
      num.cores = detectCores()
    }
      
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)
    
  
    
    if (is.null(shift.peak))
    {
      registration_NOalignment <- foreach(k=n.clust) %dopar%
      {
        .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
              object$width_spline, seeds[1:k]-1, 'F', k, weight, alpha, p, t.max, verb) 
          #NB the c++ codification of the position on vectors is from 0
      }

      registration_shift <- foreach(k=n.clust) %dopar%
      {
        .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
              object$width_spline, seeds[1:k]-1, 'T', k, weight, alpha, p, t.max, verb)
      }
    }else
    {
      if (shift.peak)
      {
        registration_NOalignment <- NULL

        registration_shift <- foreach(k=n.clust) %dopar%
        {
          .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
                object$width_spline, seeds[1:k]-1, 'T', k, weight, alpha, p, t.max, verb)
        }
      }
      if (!shift.peak)
      {
        registration_NOalignment <- foreach(k=n.clust) %dopar%
        {
          .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
                object$width_spline, seeds[1:k]-1, 'F', k, weight, alpha, p, t.max, verb)
        }

        registration_shift <- NULL
      }
    }

    stopCluster(cl)
    
  }else
  {
    if (is.null(shift.peak))
    {
      registration_NOalignment <- foreach(k=n.clust) %do%
      {
        .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
              object$width_spline, seeds[1:k]-1, 'F', k, weight, alpha, p, t.max, verb)
      }

      registration_shift <- foreach(k=n.clust) %do%
      {
        .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
              object$width_spline, seeds[1:k]-1, 'T', k, weight, alpha, p, t.max, verb)
      }
    }else
    {
      if (shift.peak)
      {
        registration_NOalignment <- NULL

        registration_shift <- foreach(k=n.clust) %do%
        {
          .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
                object$width_spline, seeds[1:k]-1, 'T', k, weight, alpha, p, t.max, verb)
        }
      }
      if (!shift.peak)
      {
        registration_NOalignment <- foreach(k=n.clust) %do%
        {
          .Call(kmean_function, x_matrix, spline_matrix , spline_der_matrix,
                object$width_spline, seeds[1:k]-1, 'F', k, weight, alpha, p, t.max, verb)
        }

        registration_shift <- NULL
      }
    }

  }

  # plot
  if (plot.graph.k)
  {
    ylim <- NULL
    no_al <- FALSE
    shi <- FALSE
    if (!is.null(registration_NOalignment))
    {
      mean_dist_NOal <- sapply(registration_NOalignment, function(x){mean(x$distances)})
      ylim <- c(min(ylim, mean_dist_NOal), max(ylim, mean_dist_NOal))
      no_al <- TRUE
    }
    if(!is.null(registration_shift))
    {
      mean_dist_shi <- sapply(registration_shift, function(x){mean(x$distances)})
      ylim <- c(min(ylim, mean_dist_shi), max(ylim, mean_dist_shi))
      shi <- TRUE
    }
    if (no_al & shi)
    {
      plot(n.clust, mean_dist_NOal, col=1, pch=19, type='b', xlab='n.clust', ylab ='average distance', ylim=ylim)
      points(n.clust, mean_dist_shi, col=2, pch=19, type='b')
      legend('topright', legend=c('NO alignment', 'shift'), col=c(1,2), pch=19, lty=1)
    }

    if (no_al & !shi)
    {
      plot(n.clust, mean_dist_NOal, col=1, pch=19, type='b', xlab='n.clust', ylab ='average distance', ylim=ylim)
      legend('topright', legend=c('NO alignment'), col=c(1), pch=19, lty=1)
    }

    if (!no_al & shi)
    {
      plot(n.clust, mean_dist_shi, col=2, pch=19, type='b', xlab='n.clust', ylab ='average distance', ylim=ylim)
      legend('topright', legend=c( 'shift'), col=c(2), pch=19, lty=1)
    }

  }

  # save the results
  if (!is.null(registration_NOalignment))
  {
    label <- lapply(registration_NOalignment, function(x){x$labels+1})
    labels_NOalignment <- convert_vectors_to_list(label, n.clust)

    dist <- lapply(registration_NOalignment, function(x){x$distances})
    dist_NOalignment <- convert_vectors_to_list(dist, n.clust)

    elementMetadata(object)[['cluster_NOshift']] <- labels_NOalignment
    elementMetadata(object)[['dist_NOshift']] <- dist_NOalignment

  }

  if (!is.null(registration_shift))
  {
    label <- lapply(registration_shift, function(x){x$labels+1})
    labels_shift <- convert_vectors_to_list(label, n.clust)

    coef <-  lapply(registration_shift, function(x){x$shift})
    coef_shift <- convert_vectors_to_list(coef, n.clust)

    dist <- lapply(registration_shift, function(x){x$distances})
    dist_shift <- convert_vectors_to_list(dist, n.clust)

    elementMetadata(object)[['cluster_shift']] <- labels_shift
    elementMetadata(object)[['coef_shift']] <- coef_shift
    elementMetadata(object)[['dist_shift']] <- dist_shift
  }

  return(object)
  }


setMethod("cluster_peak", signature=(object="GRanges") ,function(object, parallel = FALSE, num.cores = NULL,
                                                                 n.clust = NULL,  seeds = NULL, shift.peak = NULL,
                                                                 weight = NULL, subsample.weight = 100,
                                                                 alpha = 1, p = 1, t.max = 0.5, plot.graph.k = TRUE,
                                                                 verbose = TRUE)
                                                        cluster.GRange(object, parallel, num.cores,
                                                                       n.clust,  seeds, shift.peak,
                                                                       weight, subsample.weight, alpha,
                                                                       p, t.max,  plot.graph.k, verbose ) )






# function which defines from the list of the vectors of the counts
# to a matrix (n x max(length(counts)) ) with in each row the values 
# of the counts

unlist.counts <- function(counts, lenght.counts)
{
 matrix_peaks <- matrix(NA, nrow=length(lenght.counts), ncol=max(lenght.counts))
 for(i in 1:length(lenght.counts))
 {
   matrix_peaks[i,1:lenght.counts[i]] <- counts[[i]]
 }
 return(matrix_peaks)
}



# given a vector of positions (length p)
# and a list of p vectors (with fixed length n.elem)
# we aim to create a list of n.elem vectors each one of length
# max(positions) containing in the position i (i in positions)
# the value of vector[[n]][i] (n in 1:n.elem). The positions not 
# included in positions will be NA
convert_vectors_to_list <- function(vect, positions)
{
 matrix_data <- matrix(unlist(vect), length(vect), length(vect[[1]]), byrow=TRUE)
 matrix_NA_data <- matrix(NA, max(positions), length(vect[[1]]))

 matrix_NA_data[positions,] <- matrix_data

 list_temp <- apply(matrix_NA_data, 2, list)
 list <- lapply(list_temp, unlist)

 return(list)
}