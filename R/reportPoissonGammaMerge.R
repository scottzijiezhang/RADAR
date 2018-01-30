#' @title reportPoissonGammaMerge
#' @param x The RNADMethyl data list
#' @param cutoff The p_value cutoff to merge the bin
#' @param est By default "auto", the function takes estimates from given The RNADMethyl data list. One can also pass a test estimates by assigning it to stats.
reportPoissonGammaMerge <- function(x,
                                    cutoff,
                                    est = "auto"
                                    ){
  
  if( !is.matrix(est) & !is.data.frame(est) ){
    if(est == "auto"){
      stats <- x$all.est
    }else{ stop("Cannot recognize parameter est...") }
  }else if(all(c("p_value","beta","padj") %in% colnames(est)) | all(c("p_value3","beta1","padj") %in% colnames(est)) ){
    stats <- est
  }else{
    stop("The est must have column named p_value beta and padj ")
  }
  
  
  cat("Getting significant bins....\n")

  if("p_value" %in% colnames(stats)){
    sig.bins <- rownames(stats[stats[,"padj"] < cutoff ,])
  }else if("p_value3" %in% colnames(stats) ){
    sig.bins <- rownames(stats[stats[,"padj"] < cutoff ,])
    colnames(stats)[which(colnames(stats) == "p_value3")] = "p_value"
    colnames(stats)[which(colnames(stats) == "beta1")] = "beta"
  }else{
    stop("Cannot find p value column in the given Stats...\n")
  }

  ## Get the gene names of significant bins
  aa <- strsplit(sig.bins, ",")
  gene.name <- unique(unlist(lapply(aa, function(x){
    return(x[1])
  })))

  cat(paste("Reporting ",length(sig.bins)," bins in ",length(gene.name)," genes...\n"))
  ## Get the gene and continuous bins for significant bins
  geneBins <- .getGeneBins(x$geneModel,gene.name,x$binSize )
  rownames(geneBins) <- paste(geneBins$geneName,geneBins$slidingStart,sep = ",")
  geneBins$diff <- FALSE
  geneBins[sig.bins,"diff"] <- TRUE
  
  ID <- geneBins$diff
  num_lines <- length(ID)
  
  # start ids of checkpoints
  ## find peak-starting checkpoints (either from nonpeak to peak, or peak in a different batch)
  start_id <- which((ID[2:num_lines]-ID[1:num_lines-1]==1) |
                      ((geneBins$geneName[2:num_lines]!=geneBins$geneName[1:num_lines-1]) & (ID[2:num_lines] == TRUE)) )
  start_id <- start_id + 1 # add 1 since ID was counted from 2 to num_lines
  if ( ID[1]==TRUE ) { start_id <- c(1,start_id) } # if the first checkpoint bin is peak

  # end ids of checkpoints
  ## find peak-ending checkpoints (either from peak to nonpeak, or peak in a different batch)
  end_id <- which((ID[1:num_lines-1]-ID[2:num_lines]==1) |
                    ((geneBins$geneName[1:num_lines-1]!=geneBins$geneName[2:num_lines]) & (ID[1:num_lines-1] == TRUE)) )
  if (ID[num_lines]==TRUE) {end_id <- c(end_id,num_lines)} # if the last checkpoint bin is peak

  peak_id_pairs <- cbind(start_id, end_id)

  num_peaks <- nrow(peak_id_pairs)
  cat(paste("Reporting ",num_peaks," peaks...\n"))
  if (num_peaks == 0){return(data.frame())
  }else {
    start_time <- Sys.time()
    registerDoParallel(cores = 6)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
    merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
      peak_row_id <- peak_id_pairs[p,]
      tmp=GRanges(seqnames = geneBins$chr[peak_row_id[1]] , ranges = IRanges(geneBins$start[peak_row_id[1]],geneBins$end[peak_row_id[2]]),strand = geneBins$strand[peak_row_id[1]])
      tmp=GenomicRanges::intersect(tmp,x$geneModel[geneBins$geneName[peak_row_id[1]]][[1]])
      data.frame(chr=geneBins$chr[peak_row_id[1]],
                 start = geneBins$start[peak_row_id[1]],
                 end = geneBins$end[peak_row_id[2]],
                 name = as.character(geneBins$geneName[peak_row_id[1]]),
                 score = 0,
                 strand = geneBins$strand[peak_row_id[1]],
                 thickStart = geneBins$start[peak_row_id[1]],
                 thickEnd = geneBins$end[peak_row_id[2]],
                 itemRgb=0,
                 blockCount = length(tmp),
                 blockSizes = paste(as.data.frame(tmp)[,4],collapse=","),
                 blockStarts = paste(as.data.frame(tmp)[,2] - replicate(length(tmp),geneBins$start[peak_row_id[1]]),collapse=","),
                 logFC = min( stats[rownames(geneBins[peak_row_id[1]:peak_row_id[2],]), "beta"] ),
                 p_value = .fishersMethod(  stats[rownames(geneBins[peak_row_id[1]:peak_row_id[2],]), "p_value"] )
                 )
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  }
  return( merged.report )
}

## a helper function to calculate merged p_values
.fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)


