setGeneric("pileup_peak", function(object, ...) setGeneric("pileup_peak"))

pileup.GRange <- function( object, bamf , d)
{

  if (is.null(d))
  {
      stop('d has not been defined. Use compute_fragments_length or provide a value.')
  }
  if (class(object) != "GRanges")
  {
    stop('The first object is not a GRanges object.')
  }
    
  if (is.null(bamf))
  {
    stop ('provide the name of the .bam file. In the same folder there must be also the index file.')
  }
  

  # 1) EXTEND PEAKS BY d


  GR2 <- object
  start(GR2) = start(object)-d
  end(GR2) = end(object)+d

  # 2) SUBSET BAM FILE TO EXTENDED PEAKS
  bv = BamViews(bamf, bamRanges=GR2)
  GA = readGAlignments(bv)[[1]]

  # 3) EXTEND READS IN SUBSET BAM FILE
  bam_GR = granges(GA)
  bam_GR2 = resize(bam_GR, width=d)

  # 4) PILEUP ON EXTENDED READS

  cov = coverage(bam_GR2)

  chrs = unique(seqnames(object))
  counts = list()
  for(chrom in chrs){
  	p = as.vector(cov[[chrom]])
  	ind = which(as.vector(seqnames(object)) == chrom)	
  	for(j in 1:length(ind)){
  		sj = start(object)[ind[j]]
  		ej = end(object)[ind[j]]
  		counts[[ind[j]]] = p[sj:ej]		
  	}	
  }

  elementMetadata(object)[["counts"]] <- counts
  
  
  # check the lenght of the counts and the width of the peaks. They must be the same
  
  length_1 <- as.vector(sapply(object$counts, length)) # lenght of the metadata count
  length_2 <- width(object) # amplitude of the peak
  

  if ( length ( which(length_2 == length_1) ) != length(length_1) )
  {
      stop('Lengths of count and widths of peak don\'t match!')
  }
  
  return(object)
}


setMethod("pileup_peak", signature = (object = "GRanges"), function(object, bamf = NULL, d = NULL)
    pileup.GRange(object, bamf, d) )

