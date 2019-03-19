#' @title DESeq2 
#' @description A wrapper function to call DESeq2
#' @param countdata Matrix of read count.
#' @param pheno Matrix defining the grouping of the samples. 
#' @param covariates Optional matrix defining the covariates.
#' @param normalizeLib Logic option to choose whether perform normalization. Default if "FALSE" for normalized count matrix.
#' @export 
DESeq2 <- function(countdata, pheno, covariates=NULL, Dispplot=FALSE,normalizeLib = F){
  library(DESeq2)
  ## check they have the numbers of samples
  if (ncol(countdata) != nrow(pheno)) {
    stop("ERROR - samples size doesn't agree with the number of phenotype")
  }
  
  # round data if there are non-integer columns
  intcol<-vector("logical")
  for(i in 1:ncol(countdata)){
    intcol<-c(intcol,is.integer(countdata[,i]))
  }
  
  if (!all(TRUE==intcol)) {
    warning("WARNING! Non-integer expression levels. Data rounded")
    countdata<-round(countdata)
  }
  
  # remove out zero count genes in all samples
  countdata.sel = countdata[rowSums(countdata)>0,]
  
  ## start analysis
  conditions<-factor(as.matrix(pheno))
  
  # create a coldata frame and instantiate the DESeqDataSet
  
  
  if(is.null(covariates)){
    print("without pc ...")
    coldata <- data.frame(row.names=colnames(countdata.sel), conditions)
    dds <- DESeqDataSetFromMatrix(countData=countdata.sel, colData=coldata, design=~conditions)
  }else{
    print("with pc ...")
    print(dim(covariates))
    coldata <- data.frame(row.names=colnames(countdata.sel), conditions, covariates)
    #coldata <- data.frame(conditions, covariates)
    #rownames(coldata) <- as.character(colnames(countdata.sel))
    #print(head(coldata))
    names <- colnames(covariates)
    
    Formula <- formula(paste("~ conditions+", paste(names, collapse="+") ) )
    dds <- DESeqDataSetFromMatrix(countData=countdata.sel, colData=coldata, design=Formula )
  }
  
  if(!normalizeLib){
    cat("Using input matrix without additional normalization...\n")
    sizeFactors(dds) <- rep(1, nrow(pheno))
  }else{
    cat("Normalize the input matrix by geometry mean implemented in DESeq2...\n")
  }
  
  # run the DESeq2 pipeline
  dds <- DESeq(dds)
  
  if(Dispplot){
    # Plot dispersions
    png("qc-dispersions.png", 1000, 1000, pointsize=20)
    plotDispEsts(dds, main="Dispersion plot")
    dev.off()
  }
  
  # get differential expression results
  res <- DESeq2::results(dds)
  
  # merge data and output
  res.DESeq2 = data.frame(GeneID=rownames(res),log2FC=round(res$log2FoldChange,digits=2),p.adjust=res$padj, pvalue=signif(res$pvalue,digits=3))
  
  return(res.DESeq2)
}

#' @title edgeR 
#' @description A wrapper function to call edgeR
#' @param countdata Matrix of read count.
#' @param pheno Matrix defining the grouping of the samples. 
#' @param covariates Optional matrix defining the covariates.
#' @param normalizeLib Logic option to choose whether perform normalization. Default if "FALSE" for normalized count matrix.
#' @export 
edgeR <- function(countdata, pheno, covariates=NULL,normalizeLib = FALSE){
  library(edgeR)
  ## check they have the numbers of samples
  if (ncol(countdata) != nrow(pheno)) {
    stop("ERROR - samples size doesn't agree with the number of phenotype")
  }
  
  # round data if there are non-integer columns
  intcol<-vector("logical")
  for(i in 1:ncol(countdata)){
    intcol<-c(intcol,is.integer(countdata[,i]))
  }
  
  if (!all(TRUE==intcol)) {
    warning("WARNING! Non-integer expression levels. Data rounded")
    countdata<-round(countdata)
  }
  
  # remove out zero count genes in all samples
  countdata.sel = countdata[rowSums(countdata)>0,]
  
  ## start analysis
  conditions<-factor(as.matrix(pheno))
  
  # creat main edgeR object
  d = DGEList(countdata.sel, group=conditions )
  
  if(is.null(covariates)){
    print("without covariates (PCs)...")
    
    if(normalizeLib){
      # do default normalisation for edgeR
      cat("Normalize library using TMM implemented in edgeR...\n")
      d <- calcNormFactors(d)
    }else{
      cat("Using input matrix without additional normalization...\n")
    }
    
    
    # estimate the common and tagwise dispersion
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    
    # determine differentially expressed genes (using exact test)
    dest <- exactTest(d)
    
  }else{
    print("with covariates (PCs) ... ")
    design <- model.matrix(~conditions + as.matrix(covariates))
    d <- calcNormFactors(d)
    d <- estimateGLMCommonDisp(d,design)
    d <- estimateGLMTagwiseDisp(d,design)
    # using exact test to calculate the p-value
    # edgeR.et <- exactTest(edgeR.res)
    # using Liklihood ratio test to calcluate the p-value
    d.fit <- glmFit(d, design)
    dest <- glmLRT(d.fit,coef=2)
  }
  # need to extract table data from de.com object,
  # then select only required columns, plus calculate adjust p-value
  res.edgeR <- data.frame(rownames(dest$table),round(dest$table$logFC,digits=2),signif(p.adjust(dest$table$PValue,method='BH'),digits=3),signif(dest$table$PValue,digits=3))
  
  # name
  names(res.edgeR) <- c('GeneID','log2FC','p.adjust','pvalue')
  
  return(res.edgeR)
}
