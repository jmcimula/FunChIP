setGeneric("choose_k", function(object, ...) setGeneric("choose_k"))

choose.k.GRange <- function(object, k = NULL, shift.peak = NULL, cleaning = TRUE)
{

  if (class(object) != "GRanges")
  {
        stop('The first object is not a GRanges object.')
  }
    
    
  if (is.null(shift.peak))
  {
    stop('shift.peak must be provided. TRUE is for shift and FALSE for NO shift ' )
  }else
  {
    align <- ifelse(shift.peak, 'T', 'F')
  }
  if(is.null(k))
  {
    stop('Choose the number of cluster. ' )
  }
  
  if (align == 'T')
  {
    if  (is.null(object$cluster_shift))
    {
      stop ('shift.peak is TRUE but the GRange doesn\'t have the correspondent metadata. Run cluster_peak before!')
    }
    
    if (length(object$cluster_shift) < k)
    {
      stop ('k is higher than the maximum value of n.clust provided in cluster_peak')
    }
    
    if (is.na(object$cluster_shift[[1]][k]))
    {
      stop ('k is not considered as possible value of n.clust in cluster_peak')
    }
    
    labels_def <- sapply(object$cluster_shift, function(x){x[k]})

    
    elementMetadata(object)[['cluster']] <-  labels_def
    
  }
  
  if (align == 'F')
  {
    if  (is.null(object$cluster_NOshift))
    {
      stop ('shift.peak is FALSE but the GRange doesn\'t have the correspondent metadata. Run cluster_peak before!')
    }
    
    if (length(object$cluster_NOshift) < k)
    {
      stop ('k is higher than the maximum value of n.clust provided in cluster_peak')
    }
    
    if (is.na(object$cluster_NOshift[[1]][k]))
    {
      stop ('k is not a possible value of n.clust in cluster_peak')
    }
    
    labels_def <- sapply(object$cluster_NOshift, function(x){x[k]})
    dist_def <-  sapply(object$dist_NOshift,  function(x){x[k]})
    
    elementMetadata(object)[['cluster']] <-  labels_def
    
    
  }
  
    
  new_Grange <- object
    
  # remove all the metadata columns created in all the functions of the pacakge.
  # return if cleaning = TRUE
  
  elementMetadata(new_Grange)[["counts"]] <- NULL  
  elementMetadata(new_Grange)[["spline"]] <- NULL
  elementMetadata(new_Grange)[["spline_der"]] <- NULL
  elementMetadata(new_Grange)[["width_spline"]] <- NULL
  elementMetadata(new_Grange)[["summit_spline"]] <- NULL
  elementMetadata(new_Grange)[["cluster_NOshift"]] <- NULL
  elementMetadata(new_Grange)[["dist_NOshift"]] <- NULL
  elementMetadata(new_Grange)[["cluster_shift"]] <- NULL
  elementMetadata(new_Grange)[["coef_shift"]] <- NULL
  elementMetadata(new_Grange)[["dist_shift"]] <- NULL
  elementMetadata(new_Grange)[["start_spline"]] <- NULL
  elementMetadata(new_Grange)[["end_spline"]] <- NULL
  
  if (cleaning)
  {
      return(new_Grange)
  }else
  {
      return(object)
  }
}

setMethod("choose_k", signature = (object = "GRanges"), function(object, k = NULL, shift.peak = NULL, cleaning=TRUE)
    choose.k.GRange(object, k, shift.peak, cleaning) )

