#' @title filterBins
#' @param x The RNADMethyl data list
#' @param minCountsCutOff The minimal read count required in each bin, default is 10. This cutoff
#' @export x The RNAMethyl data list with filtered ip counts
filterBins <- function(
  x,
  minCountsCutOff = 10  ## the cutoff of min window read counts for furthur analysis
){
  
  # organized group info
  X <- as.factor(x$X)
  if(length(levels(X)) !=2 ){stop("The study design must have two groups")}
  group1 <- which(X==levels(X)[1])
  group2 <- which(X==levels(X)[2])
  ngroup1 <- length(group1)
  ngroup2 <- length(group2)
  
  ## filter the bin with low counts
  filtered <- .groupMeanFilter( x = x$ip_adjExpr , X = X, cuttoff = minCountsCutOff )
  cat(paste("Bins with average counts lower than ",minCountsCutOff," in both groups have been removed...\n"))
  ###############################
  
  keep.rowname <- rownames(filtered)
  input <- x$norm.input[keep.rowname,]
  ip <- x$norm.ip[keep.rowname,]
  input <- t(apply(input,1,noZero))
  FoldEnrichs <- ip/input
  na.flag <- apply(FoldEnrichs, 1, function(x){any(is.na(x)) | any(is.infinite(x))})
  FoldEnrichs <- FoldEnrichs[!na.flag, ]
  zero.flag <- apply(FoldEnrichs, 1, function(x){sum(x[group1] == 0) >= ngroup1 | sum(x[group2] == 0) >= ngroup2})
  FoldEnrichs <- FoldEnrichs[!zero.flag, ]
  ## generate flags for enriched windows
  enrichFlag <-vector(length=dim(FoldEnrichs)[1])
  for( i in 1:dim(FoldEnrichs)[1]){
    enrichFlag[i] <- max( tapply(FoldEnrichs[i,],X,median) ) > 1.2 # test if at least one group has median enrichment > 1.2
  }
  names(enrichFlag) <- rownames(FoldEnrichs)
  enrichedBins <- rownames(FoldEnrichs[which(enrichFlag),])
  enrich.filtered <-  filtered[enrichedBins, ] 
  cat("Filtering bins that is enriched in IP experiment....")
  
  x <- c(x,list('ip_adjExpr_filtered'=enrich.filtered))
  return(x)
}

#######################################################
## A helper function to filter out low counts window ##
#######################################################
.groupMeanFilter <- function(
  x, # the matrix to be filtered
  X, ## The vector of grouping (study design)
  cuttoff ## the cutoff of min window read counts for furthur analysis
){
  tmp <- matrix(nrow = dim(x)[1], ncol = 2)
  for( i in 1:dim(x)[1]){
    tmp[i,] <- tapply(as.numeric(x[i,]),X,mean)
  }
  return( x[which(apply(tmp,1,max) > cuttoff),] )
}

## a helper function to remove 0 from vector
noZero <- function(x){sapply(x,max,1)}
