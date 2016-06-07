
setGeneric("smooth_peak", function(object, ...) standardGeneric("smooth_peak") )

smooth.GRange <- function(object, n.breaks=100, subsample=TRUE, subsample.data=100,
                          order = 4, lambda=(10^(seq(-5,5,by=0.5))),
                          GCV.derivatives = TRUE , plot.GCV = FALSE)
{

  # check the first element is a GRange object        
  if (class(object) != "GRanges")
  {
      stop('The first object is not a GRanges object.')
  }
    
  # check the presence of counts and width in the GRobject
  if (is.null(object$counts))
  {
    stop('couts not present in the GR object. Definiton of the spline impossible')
  }
 
  length.data <- width(object)
  
  matrix_peaks <-  unlist.counts(object$counts, length.data)
  
  data_extended <- t(apply(matrix_peaks, 1, extension_with_min))
  data_no_background <- t(apply(data_extended, 1, function(x){x-min(x,na.rm=TRUE)} ) )
  
  # extend with 0 also befor the starting point to have a spline approximation
  # smooth on the two sides of the peaks and not juns on the right
  
  Nstart_zeros <- 200
  Nend_zeros <- 200
  data <- cbind(matrix(0, dim(data_extended)[1], Nstart_zeros), data_no_background,
                matrix(0, dim(data_extended)[1], Nend_zeros))
  
  # defalut lambda is NULL, but it can be a vector or a single value.
  # if it is a single value this is the value of lamda to be used for the smoothing
  # if it is a vector the function choice_lambda_spline identifies the correct value.
  
  if (length(lambda)!=1)
    {
      lambda_selected <- choice_lambda_spline(data, length.data+Nend_zeros, Nstart_zeros,
                                              subsample = subsample, subsample.data = subsample.data,
                                              order = order, lambda = lambda,  n.breaks = n.breaks,
                                              GCV.derivatives = GCV.derivatives , plot.GCV = plot.GCV)
  }else
  {
      lambda_selected <- lambda
  }
  
  
  smoothing <- definition_spline(data, length.data+Nend_zeros,
                                 order = order, lambda = lambda_selected,  n.breaks = n.breaks)

  spline_list <- lapply(seq_len(nrow(smoothing$spline)), function(i) smoothing$spline[i,])
  spline_der_list <- lapply(seq_len(nrow(smoothing$spline_der)), function(i) smoothing$spline_der[i,])
  
  thres <- 0.1
  
  # function which has as input
  #             - the vector x of the spline approximation
  #             - the number of non zero points of the input y
  #             (global variable Nstart_zeros the number of zeros added in the first part of the spline)
  
  
  isolate_non_zeros <- function(x, y)
  {
      minimum <- min(which(x>thres))
      maximum <- min(which(x[(y+Nstart_zeros+1): length(x)]<=thres), (length(x)-y-Nstart_zeros)) + Nstart_zeros + y
      
      return(minimum:maximum)
  }
  
  
  if(length(spline_list)==1)
  {
    warning('just one peak analyzed')
    points_non_zero_spline <- isolate_non_zeros(spline_list[[1]], length.data+Nend_zeros)
    start_spline <- start(object) - Nstart_zeros + points_non_zero_spline[1] - 1 # new starting position of the GRanges
    final_spline <- start(object) - Nstart_zeros + points_non_zero_spline[length(points_non_zero_spline)] - 1# new end of the GRanges
  }else
  {
    points_non_zero_spline <- mapply(isolate_non_zeros,
                                     spline_list, length.data) 
    # the smoothing creates new non zero values. New starting and ending position
    # of the peak                             
    start_spline <- mapply(function(x,y){x - Nstart_zeros + y[1] -1}, start(object), points_non_zero_spline)
    final_spline <- mapply(function(x,y){x - Nstart_zeros + y[length(y)] -1}, start(object), points_non_zero_spline)
  }
  if (length(spline_list)==1)
  {
    spline_list_new <- vector("list",1)
    spline_list_new[[1]] <- spline_list[[1]][points_non_zero_spline]
    spline_der_list_new <- vector("list",1)
    spline_der_list_new[[1]] <- spline_der_list[[1]][points_non_zero_spline]
    lung_new <- length(points_non_zero_spline)
  }else{
    spline_list_new <- mapply(function(x, y){x[y]}, spline_list, points_non_zero_spline)
    spline_der_list_new <- mapply(function(x, y){x[y]}, spline_der_list, points_non_zero_spline)
    lung_new <- sapply(points_non_zero_spline, length)
  }
  

  elementMetadata(object)[["spline"]] <- spline_list_new
  elementMetadata(object)[["spline_der"]] <- spline_der_list_new
  elementMetadata(object)[["width_spline"]] <- lung_new
  elementMetadata(object)[["start_spline"]] <- start_spline
  elementMetadata(object)[["end_spline"]] <- final_spline
  
  
  return(object)
  
}

setMethod("smooth_peak", signature=(object="GRanges"), function(object, n.breaks = 100, subsample = TRUE, subsample.data = 100,
                                                                order = 4, lambda = (10^(seq(-5,5,by=0.5))),
                                                                GCV.derivatives = TRUE , plot.GCV = FALSE)
  smooth.GRange (object, n.breaks, subsample, subsample.data,
                 order, lambda,
                 GCV.derivatives , plot.GCV)
)


#####################################
####### Auxiliary R funations #######
##################################### 

# funtion to indentify the correct value of lambda, using the GCV on the data
# or on the derivatives

choice_lambda_spline <- function(data, length.data, Nstart_zeros,
                                 subsample=TRUE, subsample.data=100,
                                 order = 4, lambda=(10^(seq(-5,5,by=0.5))),  n.breaks=100,
                                 GCV.derivatives = TRUE , plot.GCV = TRUE)
{
  
  if (subsample && (subsample.data > dim(data)[1]))
  {
    warning('subsample is TRUE, but subsample.data is higher than the number of data, 
            subsample automatically set to FALSE and no subsample will be performed.')
    subsample <- FALSE
  }
  if (subsample)
  {
    choosen <- sort(sample(1:dim(data)[1], subsample.data))
    data <- data[choosen, 1:(max(length.data[choosen])+Nstart_zeros)]
  }
    
  
  
  NT=dim(data)[2]
  
  num_points_projection <- n.breaks
  
  
  GCV_d0 <- rep(NA, length(lambda))
  GCV_d1 <- rep(NA, length(lambda))
  
  m=order   # order of the  spline
  
  degree=m-1  # degree
  
  spline <- matrix(NA, nrow=dim(data)[1], ncol=NT)
  spline_der <- matrix(NA, nrow=dim(data)[1], ncol=NT)
  matrix_x <- matrix(rep(1:dim(data)[2],dim(data)[1]), nrow=dim(data)[1], ncol=dim(data)[2], byrow=TRUE)
  for (j in 1:length(lambda))
  {
    print(paste('iteration on lambda: ', j))
    breaks=seq(1,NT,length=num_points_projection) 
    
    basis=create.bspline.basis(breaks,norder=m)
    
    functionalPar = fdPar(fdobj=basis,Lfdobj=2,lambda=lambda[j])
    mod=smooth.basis(argvals=1:dim(data)[2], t(as.matrix(data)),functionalPar)
    
    GCV=mod$gcv
    GCV_d0[j] <- mean(GCV)
    spline = t(eval.fd(1:dim(data)[2], mod$fd, Lfdobj=0))
    spline_der = t(eval.fd(1:(dim(data)[2]), mod$fd, Lfdobj=1))
    SSE_d1 <- rep(0, dim(data)[1])
    GCV_d1_vett <- rep(0, dim(data)[1])
    for (t in 1:dim(data)[1])
    {
        # SSE is the L2 distance between the difference
      SSE_d1[t] <- sum((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)*1 - 
        ((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)[1]/2 - 
        ((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)[NT-1]/2 #1 in the distance between 
        # two points (delta t). and then the Trapezi rule is applied
      GCV_d1_vett[t] <- NT*SSE_d1[t]/(NT-mod$df)^2
    }
    GCV_d1[j] <- mean(GCV_d1_vett)

  }
  
  
  if (plot.GCV)
  {
    par(mfrow=c(1,2))
    # plot of the GCV on the data
    plot(log10(lambda), GCV_d0, pch=19, type='b')
    
    # plot of the GCV on the data
    plot(log10(lambda), GCV_d1, pch=19, type='b')
  }
  if (GCV.derivatives){
    
    lambda_selected <- lambda[which.min(GCV_d1)]
    
  }else{
        lambda_selected <- lambda[which.min(GCV_d0)]
  }
  
  return(lambda_selected)
}


# function to define th spline, with a fixed lambda (lambda)
# it will retur the splines (and the derivatives of these splines)
# evaulated in the points of the grid

definition_spline <- function(data, length.data,
                              order = 4, lambda=0,  n.breaks=100)
{
  
  NT=dim(data)[2]
  
  num_points_projection <- n.breaks
  
  m=order   # order of the  spline
  
  degree=m-1  # degree
  
  spline <- matrix(NA, nrow=dim(data)[1], ncol=NT)
  spline_der <- matrix(NA, nrow=dim(data)[1], ncol=NT)
  matrix_x <- matrix(rep(1:dim(data)[2],dim(data)[1]), nrow=dim(data)[1], ncol=dim(data)[2], byrow=TRUE)
  
  breaks=seq(1,NT,length=num_points_projection) 
  
  basis=create.bspline.basis(breaks,norder=m)
  
  functionalPar = fdPar(fdobj=basis, Lfdobj=2, lambda=lambda)
  mod=smooth.basis(argvals=1:dim(data)[2], t(as.matrix(data)), functionalPar)
  spline = t(eval.fd(1:dim(data)[2], mod$fd, Lfdobj=0))
  spline_der = t(eval.fd(1:(dim(data)[2]), mod$fd, Lfdobj=1))
  
  return(list(spline= spline, spline_der=spline_der))
}


# function to extend a peak (stored in the vector x)
# substituting the NAs in the first and last part
# of the peak with the closer value
# e.g. c(NA, 1, 2, 3, 4, NA, NA)
# would become
# c(1, 1, 2, 3, 4, 4, 4)

extension <- function(x){
  no_nas <- which(!is.na(x))
  min_nonas <- min(no_nas)
  max_nonas <- max(no_nas)
  y <- x
  if (min_nonas!=1) 
  {
    y[1:(min_nonas -1)] <- rep(x[min_nonas], (min_nonas-1))
  }
  if (max_nonas!= length(x)) 
  {
    y[(max_nonas+1):length(x)] <- rep(x[max_nonas], (length(x)-max_nonas))
  }
  y
}


# function to extend a peak (stored in the vector x)
# substituting the NAs in the first and last part
# of the peak with the minimun value of the peak
# e.g. c(NA, 1, 2, 3, 4, NA, NA)
# would become
# c(1, 1, 2, 3, 4, 1, 1)

extension_with_min <- function(x){
  no_nas <- which(!is.na(x))
  min_nonas <- min(no_nas)
  max_nonas <- max(no_nas)
  y <- x
  if (min_nonas!=1) 
  {
    y[1:(min_nonas -1)] <- rep(min(x, na.rm=TRUE), (min_nonas-1))
  }
  if (max_nonas!= length(x)) 
  {
    y[(max_nonas+1):length(x)] <- rep(min(x, na.rm=TRUE), (length(x)-max_nonas))
  }
  y
}

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