
setGeneric("summit_peak", function(object, ... ) standardGeneric("summit_peak"))

summit.GRange <- function(object, summit = NULL)
{
    
  if (class(object) != "GRanges")
  {
        stop('The first object is not a GRanges object.')
  }
    
  if (is.null(summit))
  {
    # if the summit vector is not provided
    # the input must be a GRange object
    # with the metadata column spline. If it 
    # is not present, the summit will be computed 
    # from the counts, but it can be a noisy definition.
    # error if both spline and counts are NULL
    
    if (is.null(object$spline))
    {
        if (is.null(object$counts))
        {
            stop('Spine and counts are not metadata columns of the GRange object. 
                 Summit cannot be computed.')
            
        }else
        {
            warning('The metadata column spline is not present, the summit is 
                    computed from the row data counts.')
            
            width_peaks <- width(object)      
            
            matrix_peaks <- unlist.counts(object$counts, width_peaks)
            
            value_max <- apply(matrix_peaks, 1, max, na.rm=TRUE)
            point_max <- apply(matrix_peaks, 1, function(x){which(x == max(x, na.rm=TRUE))[1]})
            elementMetadata(object)[["summit_spline"]] <- point_max
        }

    }else
    {
        # compute the summit from the spline
          
        width_peaks <- object$width_spline    
        
        matrix_peaks <- unlist.counts(object$spline, width_peaks)
        
        value_max <- apply(matrix_peaks, 1, max, na.rm=TRUE)
        point_max <- apply(matrix_peaks, 1, function(x){which(x == max(x, na.rm=TRUE))[1]})
        elementMetadata(object)[["summit_spline"]] <- point_max
        
    }
  }else
  {
    # if summit is a vector, it will be used as the summit of the peak.
    if (length(summit)!= length(object))
    {
      stop ('summit and object must have the same length.')
    }
      
    # check that the summit is inside the ragne of the peak
    if ( ( length( which(summit > width(object)) ) != 0 ) || 
         ( length( which(summit < 0) ) != 0 ) )
    {
        stop ('summit must be inside the peak, i.e. greater than 0 and
              lower than the width of the peak')
    }
    elementMetadata(object)[["summit_spline"]] <- summit
  }
  
  return(object)
}


setMethod("summit_peak", signature=(object = "GRanges"), function(object, summit=NULL) summit.GRange(object, summit))

####################################
####### Auxiliary R funation #######
#################################### 

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

