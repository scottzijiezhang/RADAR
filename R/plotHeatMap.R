#' @title plotHeatMap
#' @param x The RADAR object
#' @param covariates The covariates matrix, must be numerical matrix.
#' @export
plotHeatMap <- function(x,covariates=NULL){
  
  if(! "all.est" %in% names(x)){stop("Need to run diffIP before plotHeatMap...")}
  
  if(is.null(covariates) ){
    topBins <-  x$ip_adjExpr_filtered[order(x$all.est[,"padj"],decreasing = F)[1:5000],]
    log_topBins <- log(topBins+1)
    log_topBins_center <- t(apply(log_topBins,1,function(x){x-mean(x)})  )
    
    dist.pear <- function(x) as.dist(1-cor(t(x)))
    hclust.ave <- function(x) hclust(x, method="average")
    gplots::heatmap.2(log_topBins_center,scale="row",trace="none",labRow=NA,main = "Top Bins ranked by p value for IP counts",
                      distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
  }else{
    topBins <-  x$ip_adjExpr_filtered[order(x$all.est[,"padj"],decreasing = F)[1:5000],]
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
    gplots::heatmap.2(cov.out,scale="row",trace="none",labRow=NA,main = "Top Bins ranked by p value for IP counts",
                      distfun=dist.pear, hclustfun=hclust.ave,col=rev(RColorBrewer::brewer.pal(9,"RdBu")))
    
  }
  
}