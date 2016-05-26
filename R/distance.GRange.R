# This function creates a list made up by two distance matrices:
# one for the Lp distance of between the data stored in the
# GRange object and one for the Lp distance between derivatives

# This function calls the cpp funcion
#  SEXP distance_matrix( SEXP x, SEXP spline, SEXP spline_der,
#                     SEXP len, SEXP p_input)
# which has as imput:
# - a Nxgrid matrix x of the aligned abscissa
# - a Nxgrid matrix spline of the spline approximation of the peaks
# - a Nxgrid matrix spline_der of the derivatived of the splines
# - a N vector len of the lengths of the splines
# - a value in {0,1,2} of the order of the distance where
#       0 means the L\infty distance
#       1 the L1 distance
#       2 the L2 distance
# and which has as output
# - a list made up of 2 NxN matrices.


distance_peak <- function(object, p = 1)
{
  
  if (class(object) != "GRanges")
  {
      stop('The first object is not a GRanges object.')
  }
    
        
  if (is.null(object$spline) || is.null(object$spline_der) || is.null(object$width_spline) )
  {
    stop('spline is not provided. Run spline_peak before computing the distance')
    object <- smooth_peak(object)
  }

  if (length(object)==1)
  {
        warning('only one peak provided. Distance from itself is 0')
      return (list(dist_matrix_d0 = 0, dist_matrix_d1 = 0))
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

  if (p != 0 & p!= 1 & p!=2)
  {
    stop ('invalid value for p, It must be 0,1,2 ')
  }


  distances <- .Call("distance_matrix", x_matrix, spline_matrix, spline_der_matrix, object$width_spline, p)

  dist_matrix_d0 <- matrix(unlist(distances$dist_D0), n.data, n.data)
  dist_matrix_d1 <- matrix(unlist(distances$dist_D1), n.data, n.data)

  return (list(dist_matrix_d0 = dist_matrix_d0, dist_matrix_d1 = dist_matrix_d1))
}



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

