
setGeneric("plot_peak", function(object, ...) setGeneric("plot_peak") )

plot.GRange <-  function(object, index=1:length(object), 
                         line.plot = 'spline', col = NULL, shift = NULL, 
                         k = NULL,  cluster.peak = FALSE)
{
  
    
  if (is.null(object$counts))
  {
    stop('No information on the peaks provided!')
  }
  
  if(max(index)>length(object) || min(index)<=0)
  {
    stop('invalid values of index. indicies must be positive and lower than the number of peaks')
  }
  
  if ( ( line.plot != 'spline' ) && ( line.plot != 'count' ) && ( line.plot != 'both' ))
  {
      stop('invalid value for line.plot It must be \'spline\', \'count\' or \'both\'')
  }
    
    
 
    # set the color.
    # if cluster.peak = FALSE
    # rainbow palette if col = NULL
    # if a vector of the same lenght as index, these are the colors used
    # if the length is differnt only the first element used.
    # 
    # if cluster.peak = TRUE
    # rainbow palette with size the numeber of clusters if col = NULL
    # if a vector of the same length as 
  if (!cluster.peak)
  {
      if (is.null(col)) 
      {
          colour <- rainbow(length(index))
      }else
      {
          if (length(col)!= length(index))
          {
              warning('number of colours different from number of peaks to be plotted.
                      All the peaks plot with the first colour.')
              colour <- rep(col[1], length(index))
          }else
          {
              colour <- col  
          }
          
      }
  }else
  {
      
      if (is.null(k))
      {
          stop ('set the number of clusters k.')
      }
      
      if (is.null(col))
      {
          colour <- rainbow(k)
      }else
      {
          if (length(col)==k)
          {
              colour <- col
          }else
          {
              warning('lenght of color not equal to the number of clusters. 
                      Value is ignored. Rainbow palette is used')
              colour <- rainbow(k)
          }
      }
  }
 
  
 if(!cluster.peak)
 {
     if (is.null(shift))
     {
         if (is.null(object$summit_spline))
         {
             shift <- FALSE
         }else
         {
             shift <- TRUE
         }
     }
     
     if (shift)
     {
         if (is.null(object$summit_spline))
         {
             stop('shift is TRUE, but summit_spline not provided')
         }else
         {
             summit_here <- object$summit_spline
         }
         plot_row_data(object, summit_here, line.plot, 
                       index, shift = TRUE, main='peaks', colour)
     }else
     {
         plot_row_data(object, summit_here=NULL, line.plot, 
                       index, shift = FALSE, main='peaks', colour)
     }
     
 }else
 {
     ## definition of the two objects matching index and labels
     if (is.null(shift))
     {
         if (!is.null(object$cluster_NOshift) & !is.null(object$cluster_shift))
         {
             stop('To plot the clustering results choose whether to plot the 
                  shifted or not shifted classification')
         }else
         {
             if (is.null(object$cluster_NOshift))
             {
                 warning('shift is NULL, shift result are plotted since no other results
                         are provided')
                 shift <- TRUE
             }else
             {
                 warning('shift is NULL, NOshift result are plotted since no other results
                             are provided')
                 shift <- FALSE
                
             }
         }
     }
     
     if (shift)
     {
         summit_here <- object$summit_spline - sapply(object$coef_shift, function(x){x[k]})
         labels_list <- object$cluster_shift
         labels <- sapply(labels_list, function(x){x[k]})
     }else
     {
         summit_here <- object$summit_spline
         labels_list <- object$cluster_NOshift
         labels <- sapply(labels_list, function(x){x[k]})
     }
  
     number_elements_cluster <- sapply(as.vector(table(labels)),
                                       function(x){min(x, max(index))})
     elements_selected <- which(labels==1)[1:number_elements_cluster[1]]
     if (k > 1)
     {
         for (i in 2:k)
         {
             elements_selected <- c(elements_selected, 
                                which(labels==i)[1:number_elements_cluster[i]])
         }
     }
     
     object <- object[elements_selected]
     summit_here <- summit_here[elements_selected]
     labels <- labels[elements_selected]
     
     if (k!=1)
     {
         if ( (k%%3) == 0)
         {
             par (mfrow = c(ceiling(k/3), 3), oma = c(5,0,0,0) + 0.5,
                  mar = c(1,1,1,1) + 0.5)
         }else
         {
             par (mfrow = c(ceiling(k/2), 2), oma = c(5,0,0,0) + 0.5,
                  mar = c(1,1,1,1) + 0.5)
         }
     }
     
    
    
     for (i in 1:k)
     {

         order <- c(which(labels != i), which(labels == i))
         col_cluster <- c(rep('grey', length(which(labels != i))),
                               rep(colour[i], length(which(labels == i))))
         if (k != 1)
         {
             lwd_cluster <- c(rep(2, length(which(labels != i))),
                                                     rep(3, length(which(labels == i))))
          }else
          {
              lwd_cluster <- rep(2, length(labels))
          }                                           
         
         plot_row_data(object, summit_here= summit_here, line.plot, 
             index=order, shift = TRUE, main=paste('cluster ', i), colour = col_cluster,
             lwd_here = lwd_cluster)
         
     }
     
 }

}
    
setMethod("plot_peak", signature=(object="GRanges") , function(object, index=1:length(object), 
                                                               line.plot = 'spline', col = NULL, shift = NULL, 
                                                               k = NULL,  cluster.peak = FALSE) 
                                                            plot.GRange(object, index, 
                                                                        line.plot, col, shift, 
                                                                        k,  cluster.peak ))



#### AUXILIARY: PLOT
#### 
#### 

plot_row_data <- function(object, summit_here=NULL, line.plot, index, shift, main='peaks', colour, lwd_here=NULL)
{
    
    if (is.null(lwd_here))
    {
        lwd_here <- rep(2, length(object))
    }
    if (is.null(object$spline)) # no information on the smoothing -> plot the original counts.
    {
        
        if (( line.plot == 'spline' ) || ( line.plot == 'both' ) )
        {
            stop('Metadata spline not provided, but line.plot is \'spline\' or \'both\'')
        }
        
        
        
        
        to_plot <- lapply(object$counts[index], function(x){x-min(x)}) # plot the data once the background has been removed
        
        ylim <- c(min(sapply(to_plot, min)), max(sapply(to_plot, max)))
        
        warning('No smoothing information. Row data are plot')
        
        
        for (i in 1:length(index))
        {
            if (i==1) plot(to_plot[[i]], type='l', col=colour[i], lty=i, ylim=ylim, xlim=c(0, max(sapply(to_plot, length))), xlab='bp', ylab='counts', main='original peaks', lwd = lwd_here[i])
            else lines(to_plot[[i]], col=colour[i], lty=i, lwd=lwd_here[i])
        }
        
    }else
    { # information on the smoothing are provided -> the spline and/or the counts are plot (see line.plot)
        
        to_plot_counts <- lapply(object$counts[index], function(x){x-min(x)}) 
        to_plot_spline <- lapply(object$spline[index], function(x){x})
        
        # plot the data once the background has been removed
        
        # definition of the aligned abscissa (0 is in the summit, if present)
        
        if (length(index)==1)
        {
            x_centered_list_spline <- vector('list', 1)
            x_centered_list_counts <- vector('list', 1)
            
            if(shift)
            {
                # convert the summit of the spline into the summit of the counts
                summit_counts <- summit_here[index] + object$start_spline[index] - start(object)[index]
                
                x_centered_list_spline[[1]] <- (-summit_here[index] +1) : 
                    (object$width_spline[index] - summit_here[index])
                x_centered_list_counts[[1]] <- (-summit_counts +1) : 
                    (width(object)[index] - summit_counts)
            }else
            {
                
                x_centered_list_spline[[1]] <- (1:object$width_spline)[index]
                x_centered_list_counts[[1]] <- (1:width(object))[index]
            }
        }else
        {
            if(shift)
            {
                summit_counts <- summit_here[index] + object$start_spline[index] - start(object)[index]
                
                x_centered_list_spline <- mapply(function(x,y){(-x+1):(y-x)}, summit_here[index], 
                                                 object$width_spline[index])
                x_centered_list_counts <- mapply(function(x,y){(-x+1):(y-x)}, summit_counts, 
                                                 width(object)[index])
                
                
            }else
            {
                x_centered_list_spline <- lapply(object$width_spline[index], function(x){1:x})
                x_centered_list_counts <- lapply(width(object)[index], function(x){1:x})
            }
            
        }
        
        if (length(index)==1)
        {
            
            if (line.plot=='both')
            {
                ylim <- c(min(to_plot_counts[[1]], to_plot_spline[[1]]), 
                          max(to_plot_counts[[1]], to_plot_spline[[1]]))
            }
            if (line.plot=='spline')
            {
                ylim <- c(min(to_plot_spline[[1]]), 
                          max(to_plot_spline[[1]]))
            }
            if (line.plot=='count')
            {
                ylim <- c(min(to_plot_counts[[1]]), 
                          max(to_plot_counts[[1]]))
            }
            
            xlim <- c(min(x_centered_list_counts[[1]], x_centered_list_spline[[1]] ),
                      max(x_centered_list_counts[[1]], x_centered_list_spline[[1]]  ))
        }else{
            
            
            if (line.plot=='both')
            {
                ylim <- c(min(sapply(to_plot_counts, min), sapply(to_plot_spline, min )), 
                          max(sapply(to_plot_counts, max), sapply(to_plot_spline, max ) )) 
            }
            if (line.plot=='spline')
            {
                ylim <- c(min(sapply(to_plot_spline, min )), 
                          max(sapply(to_plot_spline, max ) )) 
            }
            if (line.plot=='count')
            {
                ylim <- c(min(sapply(to_plot_counts, min)), 
                          max(sapply(to_plot_counts, max))) 
            }
            
            
            xlim <- c(min(sapply(x_centered_list_counts, min), sapply(x_centered_list_spline, min )), 
                      max(sapply(x_centered_list_counts, max), sapply(x_centered_list_spline, max ) )) 
            
            
        }
        
        for (i in 1:length(index))
        {
            
            if (i==1)
            {
                if (line.plot=='both')
                {
                    plot(x_centered_list_counts[[i]], to_plot_counts[[i]], type='l', col=colour[i], 
                          lty=1, ylim=ylim, xlim=xlim,  xlab='bp', ylab='counts', main=main, lwd=lwd_here[i])
                    lines(x_centered_list_spline[[i]], to_plot_spline[[i]], col=colour[i], lty=1, lwd=lwd_here[i])
                }
                if (line.plot=='spline')
                {
                    plot(x_centered_list_spline[[i]], to_plot_spline[[i]], type='l', col=colour[i], 
                          lty=1, ylim=ylim, xlim=xlim,  xlab='bp', ylab='smoothed counts', main=main, lwd=lwd_here[i])
                }
                if (line.plot=='count')
                {
                    
                    plot(x_centered_list_counts[[i]], to_plot_counts[[i]], type='l', col=colour[i], 
                          lty=1, ylim=ylim, xlim=xlim,  xlab='bp', ylab='counts', main=main, lwd=lwd_here[i])
                }
            }else
            {
                if ( (line.plot == 'spline') || (line.plot == 'both'))
                {
                    lines(x_centered_list_spline[[i]], to_plot_spline[[i]], col=colour[i], 
                         lty=1, lwd=lwd_here[i])
                }
                if ((line.plot == 'count') || (line.plot == 'both'))
                {
                    lines(x_centered_list_counts[[i]], to_plot_counts[[i]], col=colour[i], 
                          lty=1, lwd=lwd_here[i])
                }
            }
            
            
        }
    }
    
}