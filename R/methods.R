#################################################################################################################################################3

################################################################################################################################################################
#' @export
#' @rdname IP.files
setMethod("IP.files", signature("MeRIP"), function(object){
  object@bamPath.ip
})

#' @export
#' @rdname Input.files
setMethod("Input.files", signature("MeRIP"), function(object){
  object@bamPath.input
})

#' @export
#' @rdname counts
setMethod("counts", signature("MeRIP.RADAR"), function(object){
  object@reads
})

#' @export
#' @rdname sizeFactors
setMethod("sizeFactors", signature("MeRIP.RADAR"), function(object){
  object@sizeFactor
})


#######################################################################################################################################################################################################################


## helper function
.peakExons <- function(peak,y){
  exonID <- peak$start <= y$end & peak$end >= y$start
  if(sum(exonID) == 1){
    return(data.frame(start = peak$start, end = peak$end, width = peak$end - peak$start + 1))
  }else if(sum(exonID) > 1){
    peakexon <- y[exonID,]
    peakexon[1,"start"] <- peak$start
    peakexon[sum(exonID),"end"] <- peak$end
    return(data.frame(start = peakexon$start, end = peakexon$end, width = peakexon$end - peakexon$start + 1))
  }
}

## helper function
.getPeakBins <- function(geneGRList,geneName,slidingStarts,binSize){
  
  geneModel =GenomicRanges::reduce( geneGRList[geneName][[1]] )## merge overlapping exons
  
  # DNA location to gene location conversion
  df.geneModel= as.data.frame(geneModel) ##data frame of gene model
  dna.range = as.data.frame(range(geneModel) )
  df.geneModel$end = df.geneModel$end - dna.range$start + 1
  df.geneModel$start = df.geneModel$start - dna.range$start + 1
  DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
  no.exon = dim(df.geneModel)[1]
  for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
  exon.length = sum(DNA2RNA)
  #creat a corresponding map from RNA to DNA
  RNA2DNA = 1:exon.length
  pointer = 1
  for (j in 1:no.exon){
    RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
    pointer = pointer + df.geneModel$width[j]
  }
  RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates
  #create center points of continuous window
  if(exon.length <= binSize | (dna.range$strand == "+" & slidingStarts[2] == ( exon.length - binSize - exon.length %% binSize + 1) ) ){ # for genes that is shorter than bin size || if it is the rightmost bin (buffer bin) to be retrieved (positive strand only)
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[exon.length]  )
  }else if( dna.range$strand == "-" & slidingStarts[2] == 1 ){ # when retrieve the leftmost (buffer) bin on reverse strand
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize + exon.length %% binSize - 1 ]  )
  }else{ # retrieve regular bins
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize - 1 ] )
  }
  
  mapping$chr = as.character(dna.range$seqnames)
  return(mapping[,c("chr","start","end")])
  
}

## a helper function to calculate merged p_values
.fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)

#' @export
#' @title reportResult
#' @description merge significant bins and report result as BED12 format
#' @param 
setMethod("reportResult", signature("MeRIP.RADAR"), function(object, cutoff = 0.1, Beta_cutoff = 0.5, threads = 1){
  if( nrow(object@test.est) == 0 ){
    stop("Need to run diffIP/diffIP_parallel to test for differential methylation before report result! Currently no test statistics available...")
  }
  
  stats <- object@test.est
  num_bins <- length( which(stats[,"fdr"] < cutoff & abs(stats[,"beta"] )> Beta_cutoff) )
  if(num_bins < 1 ){
    stop("There is no bin passing the threshold...\n No differential peaks can be reported at current cutoff...")
  }else {
    sig.bins <- rownames(stats)[which(stats[,"fdr"] < cutoff & abs(stats[,"beta"] )> Beta_cutoff)]
  }
  
  geneBins <- geneBins(object)
  ID <- (rownames(geneBins) %in% sig.bins)
  
  num_lines <- length(ID)
  
  # start ids of checkpoints
  ## find peak-starting checkpoints (either from nonpeak to peak, or peak in a different batch)
  start_id <- which((ID[2:num_lines]-ID[1:num_lines-1]==1) |
                      ((geneBins$gene[2:num_lines]!=geneBins$gene[1:num_lines-1]) & (ID[2:num_lines] == TRUE)) )
  start_id <- start_id + 1 # add 1 since ID was counted from 2 to num_lines
  if ( ID[1]==TRUE ) { start_id <- c(1,start_id) } # if the first checkpoint bin is peak
  
  # end ids of checkpoints
  ## find peak-ending checkpoints (either from peak to nonpeak, or peak in a different batch)
  end_id <- which((ID[1:num_lines-1]-ID[2:num_lines]==1) |
                    ((geneBins$gene[1:num_lines-1]!=geneBins$gene[2:num_lines]) & (ID[1:num_lines-1] == TRUE)) )
  if (ID[num_lines]==TRUE) {end_id <- c(end_id,num_lines)} # if the last checkpoint bin is peak
  
  peak_id_pairs <- cbind(start_id, end_id)
  num_peaks <- nrow(peak_id_pairs)
  geneGRList <- object@geneModel
  peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])
  
  if (num_peaks == 0){return(data.frame())
  }else{
    start_time <- Sys.time()
    registerDoParallel(cores = threads)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
    merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
      peak_row_id <- peak_id_pairs[p,]
      geneExons <- GenomicRanges::reduce ( geneGRList[peakGenes[p]][[1]] )
      
      peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),object@binSize )
      
      peakE <- .peakExons(peak,as.data.frame(geneExons))
      data.frame(chr=peak$chr,
                 start = peak$start,
                 end = peak$end,
                 name = peakGenes[p],
                 score = 0,
                 strand = as.character(strand(geneExons))[1],
                 thickStart = peak$start,
                 thickEnd = peak$end,
                 itemRgb=0,
                 blockCount = nrow(peakE),
                 blockSizes = paste(peakE$width,collapse=","),
                 blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=","),
                 logFC = stats[rownames(geneBins[peak_row_id[1]:peak_row_id[2],]), "beta"][which.max(abs( stats[rownames(geneBins[peak_row_id[1]:peak_row_id[2],]), "beta"] ))],
                 p_value = .fishersMethod(  stats[rownames(geneBins[peak_row_id[1]:peak_row_id[2],]), "p_value"] )
      )
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  }
  
  cat("When merging neighboring significant bins, logFC was reported as the max logFC among these bins.\np-value of these bins were combined by Fisher's method.\n")

  object@mergedBins <- merged.report
  object@reportBin.fdr <- cutoff
  object@reportBin.logFoldChange <- Beta_cutoff

  return( object )

  
})

########################################################################################################################################################


#' @export
#' @title extract geneBins
#' @description geneBins extractor
setMethod("geneBins", signature("MeRIP"), function(object){
  if( nrow(object@geneBins) > 0 ){
    return(object@geneBins)
  }else{
    ## split gene and bin names
    aa <- strsplit(rownames(object@reads), ",")
    gene.name <- unlist(lapply(aa, function(x){
      return(x[1])
    }))
    bin.name <- unlist(lapply(aa, function(x){
      return(x[2])
    }))
    geneBins <- data.frame(gene=gene.name,bin=as.integer(bin.name))
    rownames(geneBins) <- rownames(object@reads)
    return(geneBins)
  }
})
#' @export
setMethod("geneBins", signature("MeRIP.RADAR"), function(object){
  callNextMethod()
})

##########################################################################################################################################################

#' @title normalizeLibrary
#' @description Normalized the input as RNA-seq data and normalize IP by enrichment. Specifically, we normalize ip libraries sizes so that the geometry mean of enrichment are the same.
#' @param object MeRIP.RADAR object.
#' @import DESeq2
#' @export
#' @return returns a MeRIP.RADAR object
setMethod("normalizeLibrary", signature("MeRIP"), function(object, boxPlot = TRUE){
  
  ## load data from input
  m6A <- object@reads[,(1+length(object@samplenames)):(2*length(object@samplenames))]
  input <- object@reads[,1:length(object@samplenames)]
  colnames(input) <- colnames(m6A) <-  object@samplenames
  object@geneBins <- geneBins<- geneBins(object)
  
  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,geneBins$gene,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- object@samplenames
  
  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
  
  ## compute normalized input data
  norm.input <-t( t(input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)
  
  
  ## estimate enrichment using top IP count bins
  ave.ip <- rowMeans(m6A)
  ave.top <- order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]
  
  ## Get the gene level input count for corresponding bins
  geneCounts.window <- geneSum.norm[geneBins[rownames(m6A),"gene"],]
  
  enrich <- as.data.frame(m6A[ave.top,]/geneCounts.window[ave.top,])
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]
  
  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(object@samplenames)])
  
  ## calculate normzlied ip read count
  norm.ip <-t( t(m6A)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)
  
  object@geneSum <- geneSum.norm
  outObj <- as(object, "MeRIP.RADAR" )
  outObj@norm.ip <- norm.ip
  outObj@norm.input <- norm.input
  outObj@sizeFactor <- sizeFactor
  
  if(boxPlot){
    plot.new()
    par(mfrow=c(2,2))
    boxplot(log(geneSum[rowSums(geneSum)!=0,]+1),main = "INPUT")
    boxplot(log(geneSum.norm[rowSums(geneSum.norm)!=0,]+1),main = "Normalized INPUT")
    boxplot(log(enrich[rowSums(enrich)!=0,]+0.1), main = "IP (Estimated enrichment)")
    enrich.norm <- as.data.frame(norm.ip[ave.top,]/geneCounts.window[ave.top,])
    boxplot(log(enrich.norm[rowSums(enrich.norm)!=0,]+0.1), main = "Normalized IP (estimated enrichment)")
    par(mfrow=c(1,1))
  }
  
  return(outObj)
  
})


#######################################################################################################################################################

## a helper function to remove 0 from vector
.noZero <- function(x){sapply(x,max,1)}

#' @title adjustExprLevel
#' @param object The MeRIP.RADAR object that has been normalized for library size.
#' @param adjustBy By default, adjust post-IP count by INPUT geneSum. Can also choose "pos" to use current position count to adjust for expression level.
#' @return  object The MeRIP.RADAR object now with IP-count adjusted for expression level.
#' @export
setMethod("adjustExprLevel", signature("MeRIP.RADAR"), function(object, adjustBy = "geneSum" ){
  if( nrow(object@norm.ip)<0 ){
    stop("Please normalize library size before running expression level adjustment!")
  }else if(adjustBy == "geneSum"){
    cat("Adjusting expression level using Input geneSum read count...\n")
    geneSum <- geneExpression(object)
    geneSum <- t(apply(geneSum,1,.noZero))
    gene.size <- t( apply(geneSum,1,function(x){x/mean(x)}) )
    gene.size.factor <- gene.size[as.character( geneBins(object)[rownames(object@norm.ip),"gene"] ) ,]
    ip_norm_geneSum <- object@norm.ip/gene.size.factor
    ip_norm_geneSum <- round(ip_norm_geneSum)
    object@ip_adjExpr <- ip_norm_geneSum
    return(object)
  }else if(adjustBy == "pos"){
    cat("Adjusting expression level using Input bin-level read count...\n")
    norm.input <- t(apply(object@norm.input,1,.noZero))
    pos.size <-  t( apply(norm.input,1,function(x){x/mean(x)}) )
    ip_norm_pos <- object@norm.ip/pos.size
    ip_norm_pos <- round(ip_norm_pos)
    object@ip_adjExpr <- ip_norm_pos
    return(object)
  }else{
    stop("Must specify adjustBy = \"geneSum\" or by \"pos\"...")
  }
})

############################################################################################################################################
## A helper function to filter group mean
.groupMeanFilter <- function(
  x, # the matrix to be filtered
  X, ## The vector of grouping (study design)
  cuttoff ## the cutoff of min window read counts for furthur analysis
){
  tmp <-  t( apply(x,1,tapply,X,mean) )
  
  return( x[which(apply(tmp,1,max) > cuttoff),] )
}

#' @title filterBins
#' @param x The MeRIP.RADAR object
#' @param minCountsCutOff The minimal read count required in each bin, default is 10. This cutoff
#' @export
setMethod("filterBins", signature("MeRIP.RADAR"), function(object, minCountsCutOff = 10 ){
  ## check if predictor variable has been set
  if(! nrow(variable(object) ) > 0 ){stop( "Please set the predictor variable/covariates by variable(MeRIP.RADAR_object) <- ..." )}
  
  # organized group info
  X <- variable(object)[,1]
  if( length(levels( as.factor(X) ) ) !=2 ){stop("The predictor variable must have two groups")}
  group1 <- which(X==levels( as.factor(X) )[1])
  group2 <- which(X==levels( as.factor(X) )[2])
  ngroup1 <- length(group1)
  ngroup2 <- length(group2)
  
  ## filter the bin with low counts
  cat("Filtering bins with low read counts...\n")
  keep.rowname <- rownames( .groupMeanFilter( x = object@norm.ip[rownames(object@ip_adjExpr),] , X = X, cuttoff = minCountsCutOff ) )
  filtered <- object@ip_adjExpr[keep.rowname,]
  cat(paste("Bins with average counts lower than ",minCountsCutOff," in both groups have been removed...\n"))
  ###############################
  
  input <- object@norm.input[keep.rowname,]
  ip <- object@norm.ip[keep.rowname,]
  input <- t(apply(input,1, .noZero) )
  FoldEnrichs <- ip/input
  na.flag <- apply(FoldEnrichs, 1, function(x){any(is.na(x)) | any(is.infinite(x))})
  FoldEnrichs <- FoldEnrichs[!na.flag, ]
  zero.flag <- apply(FoldEnrichs, 1, function(x){sum(x[group1] == 0) >= ngroup1 | sum(x[group2] == 0) >= ngroup2})
  FoldEnrichs <- FoldEnrichs[!zero.flag, ]
  ## generate flags for enriched windows
  enrichFlag <- vector(length=dim(FoldEnrichs)[1])
  cat("Filtering bins that is enriched in IP experiment...")
  for( i in 1:dim(FoldEnrichs)[1]){
    enrichFlag[i] <- max( tapply(FoldEnrichs[i,],X,median) ) > 1.2 # test if at least one group has median enrichment > 1.2
  }
  names(enrichFlag) <- rownames(FoldEnrichs)
  enrichedBins <- rownames(FoldEnrichs[which(enrichFlag),])
  enrich.filtered <-  filtered[enrichedBins, ] 
  
  object@ip_adjExpr_filtered <- enrich.filtered
  return(object)
  
})

########################################################################################################################################################

#' @export
#' @title Prepare coverage plot
#' @description import GTF into the MeRIP object for plot
#' @param object The MeRIP object
#' @param gtf optional gtf file if the stored path to gtf file has changed.
setMethod("PrepCoveragePlot", signature("MeRIP"), function(object , gtf = NULL){
  
  if( length(object@GTF) > 0 ){
    cat("GTF has already been imported as Granges, importing it again...\n")
  }
  ## check gtf slot
  if( file.exists(object@gtf) ){
    object@GTF <- rtracklayer::import(object@gtf, format = "gtf")
    return(object)
  }else if(! is.null(gtf) ){
    object@GTF <- rtracklayer::import( gtf, format = "gtf")
    object@gtf <- gtf
    cat(paste0("assigning new path to gtf file: ",gtf," \n"))
    return(object)
  }else{
    stop("The gtf file doesn't exist! Please suply path to gtf file in by PrepCoveragePlot(MeRIP, gtf = \"path/to/gtf\" )")
  }
})

###################################################################################################################################
#' @export
#' @title extractIP
#' @param object The MeRIP.RADAR object
#' @param normalized logical option, whether to return normalized IP read counts. Default is TRUE.
#' @param adjusted logical option, whether to return normalized and expression level adjusted IP read counts. Default is FALSE.
#' @param filtered logical option, whether to return normalized, adjusted and low-read-count-filtered count matrix. Default is FALSE.
#' @description The extractor to return the IP read count matrix of the experiment.
setMethod("extractIP", signature("MeRIP.RADAR"), function(object, normalized = TRUE, adjusted = FALSE , filtered = FALSE){
  if(nrow(object@norm.ip)>0 & normalized){
    if( adjusted ){
      if(filtered & nrow(object@ip_adjExpr_filtered)> 0 ){
        cat("Returning normalized, expression level adjusted and low-counts-filtered IP read counts.\n")
        return( object@ip_adjExpr_filtered )
      }else if(nrow(object@ip_adjExpr)>0){
        cat("Returning normalized and expression level adjusted IP read counts.\n")
        return( object@ip_adjExpr )
      }else{stop("IP read count has not yet been adjusted for expression variation. Call adjustExprLevel() first!")}
    }else{
      cat("Returning normalized IP read counts.\n")
      return( object@norm.ip )
    }
  }else{
    if(normalized){stop("IP read count has not yet been adjusted for expression variation. Call adjustExprLevel() first!")}
    cat("Returning raw IP read counts.\n")
    return( object@reads[,(1+length(object@samplenames)):(2*length(object@samplenames))] )
  }
})

#' @title extractINPUT
#' @param object The MeRIP.RADAR object
#' @param normalized logical option, whether to return normalized IP read counts. Default is TRUE
#' @description The extractor to return the INPUT read count matrix of the experiment.
#' @export
setMethod("extractInput", signature("MeRIP.RADAR"), function(object, normalized = TRUE){
  if(nrow(object@norm.input)>0 & normalized){
    cat("Returning normalized INPUT read counts.\n")
    return(object@norm.input)
  }else{
    if(normalized){stop("Trying to return normalized input read count, but not available. Please run normalizeLibrary() first...")}
    cat("Returning raw INPUT read counts.\n")
    return( object@reads[,1:length(object@samplenames) ] )
  }
})

#######################################################################################################################################
#' @title plotGeneCov
#' @param object The MeRIP. object
#' @param geneName The gene symbol to be ploted.
#' @param GTF The GRanges object containing gtf annotation. Can obtain by rtracklayer::import("file.gtf", format= "gtf").
#' @param libraryType Specify whether the library is the same or opposite strand of the original RNA molecule. Default is "opposite".
#' @param center Specify the method to calculate average coverage of each group. Could be mean or median.
#' @param ZoomIn c(start,end) The coordinate to zoom in at the gene to be ploted.
#' @param adjustExprLevel logical parameter. Specify whether normalize the two group so that they have similar expression level.
#' @param split Logical option. Whether to split plots for each individual (TRUE), or plot each group by group mean/median coverage (FALSE, default).
#' @export
setMethod("plotGeneCov", signature("MeRIP.RADAR"), function(object, geneName, libraryType = "opposite", center = mean,ZoomIn = NULL, adjustExprLevel = FALSE , split = FALSE){
  
  if(adjustExprLevel){
    adj<- object@geneSum[geneName,]/mean(object@geneSum[geneName,])
  }else{
    adj <- "none"
  }
  
  if( nrow(variable(object)) > 0 ){
  X <- factor(variable(object)[,1])
  if(split){
    plotGeneCoverageSplit(IP_BAMs = IP.files(object),
                   INPUT_BAMs =Input.files(object),
                   size.IP = object@sizeFactor$ip,
                   size.INPUT = object@sizeFactor$input,
                   X = X,
                   geneName = geneName,
                   geneModel = object@geneModel,
                   libraryType = libraryType, center = center  ,GTF = object@GTF, ZoomIn = ZoomIn, adjustExprLevel = adj, Names = samplenames(object)  )+
    theme(plot.title = element_text(hjust = 0.5,size = 18,color = "black"),legend.title =  element_text(hjust = 0.5,size = 15,color = "black"),legend.text =  element_text(size = 13.5,color = "black"))
  }else{
    plotGeneCoverage(IP_BAMs = IP.files(object),
                   INPUT_BAMs =Input.files(object),
                   size.IP = object@sizeFactor$ip,
                   size.INPUT = object@sizeFactor$input,
                   X = X,
                   geneName = geneName,
                   geneModel = object@geneModel,
                   libraryType = libraryType, center =  center  ,GTF = object@GTF,ZoomIn =  ZoomIn,adjustExprLevel =  adj )+
    theme(plot.title = element_text(hjust = 0.5,size = 18,color = "black"),legend.title =  element_text(hjust = 0.5,size = 15,color = "black"),legend.text =  element_text(size = 13.5,color = "black"))
  }
  
}else{
  
  if(split){
    plotGeneCoverageSplit(IP_BAMs = IP.files(object), 
                          INPUT_BAMs = Input.files(object),
                          size.IP = object@sizeFactor$ip,
                          size.INPUT = object@sizeFactor$input,
                          X = rep("All samples",length(object@samplenames)),
                          geneName = geneName,
                          geneModel = object$geneModel,
                          libraryType = libraryType, center = center, GTF =  object@GTF ,ZoomIn = ZoomIn, adjustExprLevel = adj,  Names = samplenames(object)  )+
      theme(plot.title = element_text(hjust = 0.5,size = 18,color = "black"),legend.position="none" )
    
  }else{
    plotGeneCoverage(IP_BAMs = IP.files(object),
                   INPUT_BAMs = Input.files(object),
                   size.IP = object@sizeFactor$ip,
                   size.INPUT = object@sizeFactor$input,
                   X = rep("All samples",length(object@samplenames)),
                   geneName = geneName,
                   geneModel = object$geneModel,
                   libraryType = libraryType, center =  center, GTF = object@GTF ,ZoomIn = ZoomIn, adjustExprLevel = adj  )+
    theme(plot.title = element_text(hjust = 0.5,size = 18,color = "black"),legend.position="none" )
  }
 
  
}
})



#######################################################################################################################################################

#' @title getTPM from geneSum
#' @name geneExpressionTMP
#' @param object The MeRIP.RADAR object
#' @param meanFragmentLength The mean length of RNA fragment (insert of RNA library). Default is 150bp.
#' @param normalize Logical indicating whether normalized TPM or raw TPM should be returned.
#' @return A data.frame of gene quantification in Transcript Per Million (reads).
#' @export
setMethod("geneExpressionTMP", signature("MeRIP.RADAR") , function(object, meanFragmentLength = 150, normalize = T){
  input <- object@reads[,1:length(object@samplenames)]
  gene.name <- geneBins(object)$gene
  geneSum <- geneExpression(object)
  colnames(geneSum) <- object@samplenames
  
  genes <- rownames(geneSum)
  
  cat("calculating gene length...\n")
  geneLength <- sapply(genes,function(yy){
    sum( as.data.frame( object@geneModel[[yy]] )$width )
  })
  
  
  effLength <- geneLength - meanFragmentLength
  effLength <- sapply(effLength,max,1) ## remove effective length smaller than 0.
  
  cat(paste0("computing TPM from read counts using mean fragment length = ",meanFragmentLength,".\n"))
  
  rate <- apply(geneSum,2,function(aa){aa/effLength})
  totalCounts <-colSums(rate)
  
  tpm <- t( t(rate)/totalCounts ) *1e6
  
  size.tpm <- DESeq2::estimateSizeFactorsForMatrix(tpm)
  tpm_norm <- t(t(tpm)/ size.tpm)
  
  if(normalize){
    return(tpm_norm)
  }else{
    return(tpm)
  }
})

###########################################################################################################################################################

#' @title plotTPM
#' @param TPM Dataframe of gene TPM
#' @param geneName The name of genes to be ploted.
#' @param group Categorical info for each sample.
#' @param logCount where to plot count at log scale
#' @export
plotTPM <- function(TPM,geneName,group,logCount = FALSE, facet_grid = FALSE){
  if( length(geneName) == 1){
    temp <- as.data.frame(t(TPM[geneName,] ) )
  }else{
    temp <- as.data.frame(TPM[geneName,] )
  }
  
  if(logCount){
    temp <- log(temp)
    colnames(temp) <- paste0(group,1:length(group))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(group)[1]
    for(i in 2:length(group)){
      temp_melt$Group[grep(unique(group)[i],temp_melt$variable)] <- unique(group)[i]
    }
    
    axis.font <- element_text( color = "black")
    if(facet_grid){
      ggplot(temp_melt, aes(x= Group,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log TPM")+facet_grid(.~ name)+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_blank(),
                           axis.title.y=element_text(size=20, color = "black", vjust=0.5, angle=90 ),
                           legend.title=element_text(size = 15,color = "black"),legend.text = element_text(size = 18, color = "black" ),
                           axis.text.x =element_blank() ,axis.text.y = element_text(size = 15,color = "black" ),
                           plot.title = element_text(size=22, color = "black", hjust=0.5,vjust=0.5 ),
                           axis.ticks.x = element_blank(),
                           strip.text.x = element_text(size = 15,color = "black") )+
        ggtitle("Gene expression level")
    }else{
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log TPM")+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_text(size=20, color = "black", hjust=0.5 ),
                           axis.title.y=element_text(size=20, color = "black", vjust=0.5, angle=90 ),
                           legend.title=element_text(size = 15,color = "black"),legend.text = element_text(size = 18, color = "black" ),
                           axis.text.x = element_text(size = 15,color = "black" ) ,axis.text.y = element_text(size = 15,color = "black" ),
                           plot.title = element_text(size=22, color = "black", hjust=0.5,vjust=0.5 ))+
        ggtitle("Gene expression level")
    }
    
  }else{
    colnames(temp) <- paste0(group,1:length(group))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(group)[1]
    for(i in 2:length(group)){
      temp_melt$Group[grep(unique(group)[i],temp_melt$variable)] <- unique(group)[i]
    }
    axis.font <- element_text( color = "black")
    if(facet_grid){
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="TPM")+facet_grid(.~ name)+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_blank(),
                           axis.title.y=element_text(size=20, color = "black", vjust=0.5, angle=90 ),
                           legend.title=element_text(size = 15,color = "black"),legend.text = element_text(size = 18, color = "black" ),
                           axis.text.x = element_blank() ,axis.text.y = element_text(size = 15,color = "black" ),
                           plot.title = element_text(size=22, color = "black", hjust=0.5,vjust=0.4 ),
                           axis.ticks.x = element_blank(),
                           strip.text.x = element_text(size = 15,color = "black") )+
        ggtitle("Gene expression level")
    }else{
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="TPM")+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_text(size=20, color = "black", hjust=0.5 ),
                           axis.title.y=element_text(size=20, color = "black", vjust=0.5, angle=90 ),
                           legend.title=element_text(size = 15,color = "black"),legend.text = element_text(size = 18, color = "black" ),
                           axis.text.x = element_text(size = 15,color = "black" ,colour = "black") ,axis.text.y = element_text(size = 15,color = "black" ),
                           plot.title = element_text(size=22, color = "black", hjust=0.5, vjust=0.4 ) )+
        ggtitle("Gene expression level")
    }
    
  }
}

##########################################################################################################################################

#' @title Perform inferential test using Poisson Random effect model in RADAR 
#' @description 
#' @name diffIP
#' @param object The MeRIP.RADAR object
#' @param fdrBy The method to control for false discovery rate. The default is "qvalue", can also be "fdr".
#' @export
setMethod( "diffIP", signature("MeRIP.RADAR"), function(object, exclude = NULL, maxPsi = 100, fdrBy = "fdr", steps = 10, gamma = 0.5, down = 0.1){
  
  if( nrow(object@variate) != length(object@samplenames)  ){
    stop(" Predictor variable lengthen needs to match the sample size! If you haven't set the predictor variable, please set it by variable(object) <- data.frmae(group = c(...)) ")
  }else if(length(unique(variable(object)[,1])) != 2 & !is.numeric(variable(object)[,1]) ){
    stop("The levels of predictor variable needs to be two!")
  }
  
  if(!is.null(exclude) & all( exclude %in% samplenames(object)) ){
    object <- select(object, setdiff(samplenames(object), exclude ))
  }
  
  allY <- object@ip_adjExpr_filtered
  psi <- 10 # start point
  
  ## convert predictor variable if it is not numeric
  if( is.numeric(variable(object)[,1]) ){
    X <- variable(object)[,1]
  }else{
    tmp <- as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] )
    names(tmp) <- as.character(variable(object)[,1])
    cat("The predictor variable has been converted:\n")
    print(tmp)
    X <- as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] ) # convert categorical variable into numerical variable.
  }
  
  if( ncol(variable(object)) == 1 ){
    
    cat("running PoissonGamma test at single beta mode\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissionGamma(Y, X, beta, psi, mu2, gamma = gamma, steps = steps, down = down,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    rownames(all.est) <- rownames(allY)[all.id]
    
    
  }else if(ncol(variable(object)) > 1){
    if(! any(sapply(2:ncol(variable(object)),function(x) is.numeric(variable(object)[,x])) ) ){stop("Please convert all covariates into numerical variables. Discrete variable should be converted to binary variable.")}
    X.all <- as.matrix( cbind(X,variable(object)[,2:ncol(variable(object))])  )
    colnames(X.all) <- colnames(variable(object))
    design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste(colnames(X.all), sep = "",collapse = " + ")) )
    cat("running PoissonGamma test at multi-beta mode...\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ] )
      aa <- unlist(summary( lm( design.multiBeta, data = as.data.frame(cbind(Y, X.all)) ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = gamma, steps = steps, down = down,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    
    rownames(all.est) <- rownames(allY)[all.id] ## assign window names to test statistics
    colnames(all.est) <- gsub("p_value3","p_value",colnames(all.est))
    colnames(all.est) <- gsub("beta1","beta",colnames(all.est))
  }
  
  if(fdrBy == "qvalue"){
    fdr <- qvalue::qvalue(all.est[,"p_value"] )$qvalue
    object@fdr.method = "qvalue"
  }else{
    fdr <- p.adjust(all.est[,"p_value"],method = fdrBy )
    object@fdr.method = "Benjamini & Hochberg"
  }
  
  
  object@test.est <- cbind(all.est,fdr)
  object@test.method <- "PoissonGamma test (RADAR)"
  cat("\n")
  return(object)
})

#' @title Perform inferential test using Poisson Random effect model in RADAR 
#' @description 
#' @name diffIP_parallel
#' @param object The MeRIP.RADAR object
#' @param exclude A vector of characters. The samples to be excluded from the differential test.
#' @param maxPsi The max random effect parameter Psi allowed in the model fitting.
#' @param fdrBy The method to control for false discovery rate. The default is "qvalue", can also be "fdr".
#' @param thread The number of thread to use for computing.
#' @export
setMethod( "diffIP_parallel", signature("MeRIP.RADAR"), function(object, exclude = NULL, maxPsi = 100, fdrBy = "fdr", thread = 8 ){
  
  if( nrow(object@variate) != length(object@samplenames)  ){
    stop(" Predictor variable lengthen needs to match the sample size! If you haven't set the predictor variable, please set it by variable(object) <- data.frmae(group = c(...)) ")
  }else if(length(unique(variable(object)[,1])) != 2 & !is.numeric(variable(object)[,1]) ){
    stop("The levels of predictor variable needs to be two!")
  }
  
  if(!is.null(exclude) & all( exclude %in% samplenames(object)) ){
    object <- select(object, setdiff(samplenames(object), exclude ))
  }
  
  allY <- object@ip_adjExpr_filtered
  psi <- 10 # start point
  
  ## convert predictor variable if it is not numeric
  if( is.numeric(variable(object)[,1]) ){
    X <- variable(object)[,1]
  }else{
    tmp <-  as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] )
    names(tmp) <- as.character(variable(object)[,1])
    cat("The predictor variable has been converted:\n")
    print(tmp)
    X <-  as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] ) # convert categorical variable into numerical variable.
  }
  
  if( ncol(variable(object)) == 1 ){
    
    cat("running PoissonGamma test at single beta mode\n")
    
    start_time <- Sys.time()  # track run time
    ## register cluster for hyperthread computing
    doParallel::registerDoParallel(cores=thread)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to run PoissonGamma test...\n"))
    
    error.id <- NULL
    all.est <- foreach(kk = 1:nrow(allY), .combine = rbind, .errorhandling = "remove") %dopar% {
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) == "try-error"){
        error.id <- c(error.id, kk)
      }
      est
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to run PoissonGamma test:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    
    all.id <- which(! 1:nrow(allY) %in% error.id )
    rownames(all.est) <- rownames(allY)[all.id]
    
    
  }else if(ncol(variable(object)) > 1){
    if(! any(sapply(2:ncol(variable(object)),function(x) is.numeric(variable(object)[,x])) ) ){stop("Please convert all covariates into numerical variables. Discrete variable should be converted to binary variable.")}
    X.all <- as.matrix( cbind(X,variable(object)[,2:ncol(variable(object))])  )
    colnames(X.all) <- colnames(variable(object))
    design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste(colnames(X.all), sep = "",collapse = " + ")) )
    cat("running PoissonGamma test at multi-beta mode...\n")
    
    error.id <- NULL
    start_time <- Sys.time()
    ## register cluster for hyperthread computing
    doParallel::registerDoParallel(cores=thread)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to run multi-beta PoissonGamma test...\n"))
    
    all1 <- foreach(kk = 1:nrow(allY),.combine = rbind, .errorhandling = "remove") %dopar% {
      Y <- unlist(allY[kk, ] )
      aa <- unlist(summary( lm( design.multiBeta, data = as.data.frame( cbind(Y, X.all) ) ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 25, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) == "try-error"){
        error.id <- c(error.id, kk)
      }
      est
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to run multi-beta PoissonGamma test:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    all.id <- which(! 1:nrow(all1) %in% error.id )
    rownames(all1) <- rownames(allY)[all.id] ## assign window names to test statistics
    
    all.est <- all1
    colnames(all.est) <- gsub("p_value3","p_value",colnames(all.est))
    colnames(all.est) <- gsub("beta1","beta",colnames(all.est))
  }
  
  if(fdrBy == "qvalue"){
    fdr <- qvalue::qvalue(all.est[,"p_value"] )$qvalue
    object@fdr.method = "qvalue"
  }else{
    fdr <- p.adjust(all.est[,"p_value"],method = fdrBy )
    object@fdr.method = "Benjamini & Hochberg"
  }
  
  object@test.est <- cbind(all.est,fdr)
  object@test.method <- "PoissonGamma test (RADAR)"
  cat("\n")
  return(object)
})

##########################################################################################################################################

#' @export
setMethod("samplenames", signature("MeRIP.RADAR"), function(object ){
  object@samplenames
})
#' @export
setMethod("samplenames", signature("MeRIP"), function(object ){
  object@samplenames
}) 
#' @export
setMethod("samplenames<-", signature("MeRIP.RADAR"), function(object ,value){
  object@samplenames <- value
  object
})
#' @export
setMethod("samplenames<-", signature("MeRIP"), function(object ,value){
  object@samplenames <- value
  object
})

###########################################################################################################################
#' @export
setMethod("variable", signature("MeRIP.RADAR"), function(object ){
  object@variate
})
#' @export
setMethod("variable<-", signature("MeRIP.RADAR"), function(object ,value){
  if(! is.data.frame(value) ){ stop("Please set the variable as data.frame where rows are samples and columns are variables.")}
  object@variate <- value
  object
})


##################################################################################################################################
#' @title extractor for RNAseq data
#' @name geneExpression
#' @description Extract gene level expression (RNAseq) data in normalized read counts
#' @param object The MeRIP object
#' @import DESeq2
#' @export
setMethod("geneExpression", signature("MeRIP.RADAR"), function( object ){
  if(nrow(object@geneSum)> 0 ){
    return(object@geneSum)
  }else{
    input <- object@reads[,1:length(object@samplenames)]
    colnames(input) <- object@samplenames
    geneBins <- geneBins(object)
    ## Get input geneSum (gene level quantification)
    geneSum <- NULL
    for(i in 1:ncol(input) ){
      y <- input[,i]
      gene.sum <- tapply(y,geneBins$gene,sum)
      geneSum <- cbind(geneSum,gene.sum)
    }
    colnames(geneSum) <- object@samplenames
    
    size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
    
    geneSum.norm <- t ( t(geneSum)/size.input)
    return(geneSum.norm)
  }
})

#########################################################################################################################

#########################################################################################################################################

#' @export
#' @title subset MeRIPdata
#' @name select
#' @description subset dataset by samples.
#' @param object The MeRIP object
#' @param samples The samplenames to be subset or the index number of samples to be subset.
#' @return an MeRIP object of selected samples.
setMethod("select", signature("MeRIP"),function(object , samples ){
  
  id <- if( is.integer(samples) & all(samples > 0) & all(samples < length(samplenames(object))) ){
    samples
  }else if( all(samples %in% samplenames(object)) ){
    match(samples , samplenames(object))
  }
  
  newOb <- new(Class = "MeRIP",
               reads = counts(object)[,c(id,id + length(samplenames(object)))],
               binSize = object@binSize,
               geneModel = object@geneModel,
               gtf = object@gtf,
               bamPath.input = Input.files(object)[id],
               bamPath.ip = IP.files(object)[id],
               samplenames = samplenames(object)[id],
               geneBins = geneBins(object),
               GTF = object@GTF)
  
  newOb@geneSum <- if( nrow(object@geneSum) > 1 ){object@geneSum[,id]}else{object@geneSum}
  return(newOb)
})


#' @export
#' @name select
#' @description subset dataset by samples.
#' @param object The MeRIP.RADAR object
#' @param samples The samplenames to be subset or the index number of samples to be subset.
#' @return an MeRIP.RADAR object of selected samples.
setMethod("select", signature("MeRIP.RADAR"),function(object , samples ){
  
  id <- if( is.integer(samples) & all(samples > 0) & all(samples < length(samplenames(object))) ){
    samples
  }else if( all(samples %in% samplenames(object)) ){
    match(samples , samplenames(object))
  }else{
    stop("Please specify valide samplenames or index to subset the dataset.")
  }
  
  newOb <- new(Class = "MeRIP.RADAR",
               reads = counts(object)[,c(id,id + length(samplenames(object)))],
               binSize = object@binSize,
               geneModel = object@geneModel,
               gtf = object@gtf,
               bamPath.input = Input.files(object)[id],
               bamPath.ip = IP.files(object)[id],
               samplenames = samplenames(object)[id],
               geneBins = geneBins(object),
               GTF = object@GTF
  )
  
  newOb@geneSum <- if( nrow(object@geneSum) > 1 ){object@geneSum[,id]}else{object@geneSum}
  newOb@norm.input <- if(nrow(object@norm.input) > 1){object@norm.input[,id]}else{object@norm.input}
  newOb@norm.ip <- if(nrow(object@norm.ip) > 1){object@norm.ip[,id]}else{object@norm.ip}
  newOb@ip_adjExpr <- if(nrow(object@ip_adjExpr) > 1 ){object@ip_adjExpr[,id]}else{object@ip_adjExpr}
  newOb@ip_adjExpr_filtered <- if(nrow(object@ip_adjExpr_filtered) > 1 ){object@ip_adjExpr_filtered[,id]}else{object@ip_adjExpr_filtered}
  newOb@sizeFactor <- if(nrow(object@sizeFactor) > 1 ){object@sizeFactor[id,]}else{object@sizeFactor}
  newOb@variate <- if(nrow(object@variate) > 1 ){ if(ncol(object@variate)>1){
    object@variate[id,]
  }else{
    newVar <- data.frame(object@variate[id,])
    colnames(newVar) <- colnames(object@variate)
    newVar
  } }else{
    object@variate
    }
  if(nrow(object@test.est) > 1 ){cat("Inferential test is not inherited because test result changes when samples are subsetted!\nPlease re-do test.\n")}
  return(newOb)
})

######################################################################################################################################

#' @title plot distribution of peaks on gene annotation
#' @name peakDistribution
#' @description plot distribution of differential peaks on gene annotation
#' @param object The MeRIP.RADAR object
#' @export
setMethod("peakDistribution", signature("MeRIP.RADAR"), function(object){
  ## collapes the peak to the center
  x <- results(object)
  for(i in 1:nrow(x)){
    start = x[i,2]
    end = x[i,3]
    blocks = x[i,10]
    blockStart = as.numeric( unlist(strsplit(as.character(x[i,12]),",")) )
    blockSize = as.numeric( unlist(strsplit(as.character(x[i,11]),",")) )
    if(blocks <2){
      x[i,2] =  x[i,3] = round(start + sum(blockSize)/2 )
    } else{
      pointer = 2
      while(sum(blockSize[1:pointer]) < sum(blockSize)/2 ){pointer = pointer+1}
      x[i,2] =  x[i,3] = round(start + blockStart[pointer] + sum(blockSize)/2 - sum(blockSize[1:(pointer-1)]) )
    }
    x[i,10] = 1
  }
  
  ################
  gr.peak <- GenomicRanges::reduce( makeGRangesFromDataFrame(x) )
  txdb <- makeTxDbFromGFF(file = object@gtf,format = "gtf")
  
  cds <-  cdsBy(txdb,by = "tx")
  n.cds <- length( which(countOverlaps(gr.peak,cds)>0) )
  
  fiveUTR <- fiveUTRsByTranscript(txdb)
  n.fiveUTR <- length( which(countOverlaps(gr.peak,fiveUTR)>0) )
  
  threeUTR <- threeUTRsByTranscript(txdb)
  n.threeUTR <-  length( which(countOverlaps(gr.peak,threeUTR)>0) )
  
  exon <- exonsBy(txdb,by = "tx")
  all_mRNA <- unique(c(names(fiveUTR),names(threeUTR),names(cds)))
  name_ncRNA <- setdiff(names(exon),all_mRNA)
  ncRNA <- exon[name_ncRNA]
  n.ncRNA <-  length( which(countOverlaps(gr.peak,ncRNA)>0) )
  
  slices <- c(n.fiveUTR,n.cds,n.threeUTR,n.ncRNA)
  pct <- round(slices/sum(slices)*100)
  lbls <- c("5'UTR","CDS","3'UTR", "ncRNA")
  lbls <- paste(lbls, pct,sep = "\n") # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  distr <- data.frame(Annotation = factor(c("5'UTR","CDS","3'UTR", "ncRNA"),levels = c("5'UTR","CDS","3'UTR", "ncRNA")),
                      fraction = pct)
  ggplot(distr,aes( x = factor(1), y = fraction, fill = Annotation)) + geom_bar(width = 1,stat = "identity") +coord_polar( theta = "y")+theme_minimal()+
    theme(axis.title = element_blank(), axis.text = element_blank(),axis.title.y = element_blank(),panel.grid=element_blank(), legend.title = element_text(color = "black"),legend.text = element_text(color = "black")) + 
    geom_text(aes(y = rev( fraction/2 + c(0, cumsum(fraction)[-length(fraction)])),label = lbls), size=4 )
})

#############################################################################################################################

#' @export
#' @name results
#' @title export results
#' @description The extractor for final test result.
#' @param object The MeRIP.RADAR object.
#' @return joint peaks (with test result) in a data.frame.
setMethod("results", signature("MeRIP.RADAR"), function(object){
  if(object@test.method != "none" & nrow(object@mergedBins)> 1){
    cat(paste0("There are ",nrow(object@mergedBins)," reported differential loci at FDR < ",object@reportBin.fdr, " and logFoldChange > ",object@reportBin.logFoldChange,".\n"))
    return( object@mergedBins )
  }else if(nrow(object@mergedBins)> 1){
    stop("Please use reportResult(MeRIP.Peak) function to report significant bins at a cutoff you set (or by default FDR<0.1 & Log fold change > 0.5)!\n")
  }else{
    stop("No inferential test has been performed!\n")
  }
})

#############################################################################################################################

#' @export
#' @name plotHeatMap
#' @title plotHeatMap
#' @description Plot heat map of methylation variations across samples been analyzed.
#' @param object The MeRIP.RADAR object.
#' @param covariates Logic parameter. Whether to regress out covariates (if variable has been set to have more than one column). Default is TRUE.
setMethod("plotHeatMap", signature("MeRIP.RADAR"), function(object, covariates=TRUE){
  if(object@test.method != "none" & nrow(object@mergedBins)> 1){
    cat(paste0("Plot heat map for differential loci at FDR < ",object@reportBin.fdr, " and logFoldChange > ",object@reportBin.logFoldChange,".\n"))
    
    stats <- object@test.est
    topBins.id <-  rownames(stats)[which(stats[,"fdr"] < object@reportBin.fdr & abs(stats[,"beta"] )> object@reportBin.logFoldChange)]
    topBins <- extractIP(object, normalized = TRUE, adjusted = T )[topBins.id,]
    log_topBins <- log(topBins+1)

    if(!covariates | ncol(variable(object) )<2 ){
      
      log_topBins_center <- t(apply(log_topBins,1,function(x){x-mean(x)})  )
      colnames(log_topBins_center) <- samplenames(object) 
      dist.pear <- function(x) as.dist(1-cor(t(x)))
      hclust.ave <- function(x) hclust(x, method="average")
      gplots::heatmap.2(log_topBins_center,scale="row",trace="none",labRow=NA,main = "Methylation level of significant bins",
                        distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
    }else{
      
      covariates <- variable(object)[,-c(1)]
      
      registerDoParallel(cores = 4)
      cov.out <- foreach(i = 1:nrow(log_topBins), .combine = rbind) %dopar% {
        Y = log_topBins[i,]
        tmp_data <- as.data.frame(cbind(Y,covariates))
        resi <- residuals( lm(Y~.,data=tmp_data  ) )
        resi
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      rownames(cov.out) <- rownames(log_topBins)
      colnames(cov.out) <- colnames(log_topBins)
      
      dist.pear <- function(x) as.dist(1-cor(t(x)))
      hclust.ave <- function(x) hclust(x, method="average")
      gplots::heatmap.2(cov.out,scale="row",trace="none",labRow=NA,main = "Methylation level of significant bins\n(Covariates regressed out)",
                        distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
      
    }
    
  }else if(nrow(object@mergedBins)> 1){
    stop("Please use reportResult(MeRIP.RADAR) function to report significant bins at a cutoff you set (or by default FDR<0.1 & Log fold change > 0.5)!\n")
  }else{
    stop("No inferential test has been performed!\n")
  }
})

###########################################################

.reportPeaks <- function(radar, est , cutoff = 0.1, Beta_cutoff = 0.5, threads = 1){
  stats <- est
  
  colnames(stats) <- gsub("fdr","padj",colnames(stats))
  colnames(stats) <- gsub("log2.OR","beta",colnames(stats))
  colnames(stats) <- gsub("logFC","beta",colnames(stats))
  colnames(stats) <- gsub("log2.OR","beta",colnames(stats))
  
  num_bins <- length( which(stats[,"fdr"] < cutoff & abs(stats[,"beta"] )> Beta_cutoff) )
  if(num_bins < 1 ){
    stop("There is no bin passing the threshold...\n No differential peaks can be reported at current cutoff...")
  }else {
    sig.bins <- rownames(stats)[which(stats[,"fdr"] < cutoff & abs(stats[,"beta"] )> Beta_cutoff)]
  }
  
}