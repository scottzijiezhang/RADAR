plotPCA.DESeq2 <- function(data,group = NULL, returnPC = FALSE){
  if(is.null(group)){
    colData <- data.frame(group=colnames(data))
  }else{
    colData <- data.frame(group=group)
  }

  rownames(colData) <- colnames(data)
  dds <- DESeq2::DESeqDataSetFromMatrix(data,colData,design = ~group)
  cat("Using regularized log transformation of DESeq2 to tranform data...\n")
  rld <- DESeq2::rlog(dds)

  if(returnPC){
    PCs <- DESeq2::plotPCA(rld,intgroup = "group",returnData=TRUE)
    return(PCs)
  }else{
    cat("Plot PCA using the rlog transformed data...\n")
    DESeq2::plotPCA(rld,intgroup = "group")
  }

}

#' @title plotPCAfromMatrix
#' @param m The matrix of count data
#' @param group The factor levels to color the samples. Should be the save number as the # of matrix columns
#' @param loglink Logic parameter determine whether to take log of the metrix data. Default is TRUE. If your input matrix is at log scale, use FALSE.
#' @export
plotPCAfromMatrix <- function(m,group,loglink = TRUE){
  if(loglink){
    mm <- log(m + 1)
  }else{
    mm <- m 
  }
  pc <- prcomp(t(mm))
  pca.df <- as.data.frame(pc$x)
  vars <- apply(pca.df ,2, var)
  props <- 100*(vars / sum(vars) )
  makeLab <- function(x,pc) paste0("PC",pc,": ",round(x,digits = 2),"% variance")
  ggplot(data = pca.df,aes(PC1,PC2,label = rownames(pca.df),colour=group) )+geom_text()+ xlab(makeLab(props[1],1)) + ylab(makeLab(props[2],2))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 1),
                       axis.title.x=element_text(size=20, hjust=0.5),
                       axis.title.y=element_text(size=20, vjust=0.4, angle=90),
                       legend.title=element_text(size = 15),legend.text = element_text(size = 18),
                       axis.text.x = element_text(size = 15,color = "black") ,axis.text.y = element_text(size = 15)  )

}
