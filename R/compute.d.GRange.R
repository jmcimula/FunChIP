compute_fragments_length <- function( object, bamf , min.d = 0, max.d = 200)
{
    
    if (class(object) != "GRanges")
    {
        stop('The first object is not a GRanges object.')
    }
    
    if (is.null(bamf))
    {
        stop ('provide the name of the .bam file. In the same folder there must be also the index file.')
    }
    
    bf=BamFile(bamf)
    # SET SCANNING PARAMTERES (THE GRANGE)
    param = ScanBamParam(which=object)
    
    # SET OUTPUT PARAMETERS (IGNORE DIFFERENT STRANDS AND NUCLEOTIDES
    p_param = PileupParam(distinguish_nucleotides=FALSE, distinguish_strands=TRUE)
    
    # OBTAIN PILEUP
    p = pileup(bf, scanBamParam=param, pileupParam=p_param)
    
    
    list_counts <-  by(p, p$which_label, function(x){
        return(list(counts = x$count, pos = x$pos, strand = x$strand))})
    
    data_grange <- mapply(function(x,y){list(start=x, width=y)}, 
                          start(object), width(object), SIMPLIFY = FALSE)
    
    st <- sapply(list_counts, function(x){x$strand})
    list_counts_correct_length_strand <- mapply(function(x, y){
        start_object <- y$start
        positions <- x$pos
        vect_pos <- rep(0, y$width)
        vect_neg <- rep(0, y$width)
        
        # isolate the positive and negative components
        position_positive <- which(x$strand == '+')
        vect_pos[positions[position_positive]-start_object+1] <- x$counts[position_positive]
        
        position_negative <- which(x$strand == '-')
        vect_neg[positions[position_negative]-start_object+1] <- x$counts[position_negative]
        
        vector_list <- vector('list', 2)
        vector_list[['vect_pos']] <- vect_pos
        vector_list[['vect_neg']] <- vect_neg
        
        return(vector_list)
        
    }, 
            list_counts, data_grange, 
            SIMPLIFY = FALSE)
    
    d_opt <- choose_d(list_counts_correct_length_strand, c(min.d, max.d), bamf)
    
    
    return(d_opt)
    
}

#################################
###### AUXILIARY FUNCTIONS ######
#################################

compute_distance <- function(x_1, y_1, x_2, y_2)
{
    x_common <- union(x_1, x_2)
    y_1_common <- rep(0, length(x_common))
    y_1_common[x_common %in% x_1] <- y_1
    
    y_2_common <- rep(0, length(x_common))
    y_2_common[x_common %in% x_2] <- y_2
    
    vect_diff_sq <- (y_2_common - y_1_common)^2
    # Trapez. rule with equispaced points (distance = 1)
    # squared of the distance normalized by the length of the domain
    L2_dist <- (2*sum(vect_diff_sq) - vect_diff_sq[1] - vect_diff_sq[length(vect_diff_sq)])/(2*length(x_common))
    
    return(L2_dist)
}


compute_distance_with_shift <- function(data_1, data_2, shift)
{
    if (length(data_1) != length(data_2))
    {
        stop('different length for the two vectors. no distance computed.')
    }
    x_1 <- 1:length(data_1)
    x_2 <- 1:length(data_1)
    
    x_2 <- x_2 - shift # if data_2 is the negative strand
    # it has to be moved decreasing the grid
    
    d <- compute_distance(x_1 , data_1, x_2, data_2)
    return(d)
}

choose_d <- function(x, range_d, bamf)
{
    l2_dist <- sapply(range_d[1]:range_d[2], function(shift)
    {
        l2_vec <- sapply(x, function(x){
            compute_distance_with_shift(x$vect_pos, x$vect_neg, shift)
        })
        return(sum(sqrt(l2_vec)))                 
        
    }
    )
    # if you need the plot of d vs distance
    plot(range_d[1] : range_d[2], l2_dist, xlab = 'dist positive - negative', ylab = 'global distance',
         main = 'L2 distance squared normalized by the length of the domain')
    
    minimum_l2 <- which.min(l2_dist)
    optimum_d <- (range_d[1] : range_d[2])[minimum_l2]
    points(optimum_d, l2_dist[minimum_l2], pch = 19, col= 2)
    
    print(paste('estimated distance positive - negative read ', optimum_d, sep=''))
    
    read_length <- round(mean(scanBam(bamf)[[1]]$qwidth, na.rm=TRUE))
    print(paste('estimated read length ', read_length, sep=''))
    
    d_opt <- optimum_d+read_length
    return(d_opt)
}


