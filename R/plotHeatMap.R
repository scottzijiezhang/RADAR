#' @title plotHeatMap
#' @param x The RADAR object
#' @param covariates The covariates matrix, must be numerical matrix.
#' @param topBins The number of top bins to plot the heat map. Default is 5000. 
#' @export
plotHeatMap <- function(x,covariates=NULL,topBins = "significant", label = NULL ){
  
  if(! "all.est" %in% names(x)){stop("Need to run diffIP before plotHeatMap...")
    }else if(!is.numeric(topBins)){
    topBins <- length(which(x$all.est[,"padj"]<0.1))
  }
  
  if(nrow(x$ip_adjExpr_filtered) < topBins){
    topBins <- nrow(x$ip_adjExpr_filtered)
    cat( paste0("Only ",topBins," bins passed filter...\n") )
    cat( paste0("Plot heat map for top ",topBins," bins ranked by differential test p value...\n") )
  }else{
    cat( paste0("Plot heat map for top ",topBins," bins ranked by differential test p value...\n") )
  }
  
  if(is.null(covariates) ){
    
    topBins <-  x$ip_adjExpr_filtered[order(x$all.est[,"padj"],decreasing = F)[1:topBins],]
    log_topBins <- log(topBins+1)
    log_topBins_center <- t(apply(log_topBins,1,function(x){x-mean(x)})  )
    if(!is.null(label) & length(label) == ncol(log_topBins_center) ){colnames(log_topBins_center) <- label }
    dist.pear <- function(x) as.dist(1-cor(t(x)))
    hclust.ave <- function(x) hclust(x, method="average")
    gplots::heatmap.2(log_topBins_center,scale="row",trace="none",labRow=NA,main = "Methylation level of top-lowest p value bins",
                      distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
  }else{
    
    topBins <-  x$ip_adjExpr_filtered[order(x$all.est[,"padj"],decreasing = F)[1:topBins],]
    log_topBins <- log(topBins+1)
    
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
    if(!is.null(label) & length(label) == ncol(cov.out) ){colnames(cov.out) <- label }
    gplots::heatmap.2(cov.out,scale="row",trace="none",labRow=NA,main = "Methylation level of top-lowest p value bins\n(Covariates regressed out)",
                      distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
    
  }
  
}