#' @title normalizeLibrary
#' @param readsOut Read counts data list from countReads() function.
#' @param X Study design or Grouping for the samples, should have 2 levels
#' @export norm_lib returns the list of normalized reads, study design and gene-bin names
normalizeLibrary <- function(readsOut,X){
  
  input <- readsOut$reads[,1:length(readsOut$samplenames)]
  m6A <- readsOut$reads[,(1+length(readsOut$samplenames)):(2*length(readsOut$samplenames))]
  colnames(input) <- colnames(m6A) <- readsOut$samplenames
  ## split gene and bin names
  aa <- strsplit(rownames(input), ",")
  gene.name <- unlist(lapply(aa, function(x){
    return(x[1])
  }))
  bin.name <- unlist(lapply(aa, function(x){
    return(x[2])
  }))
  geneBins <- data.frame(gene=gene.name,bin=bin.name)
  rownames(geneBins) <- rownames(input)
  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,gene.name,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- readsOut$samplenames
  
  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum[rowSums(geneSum)>0,])

  norm.input <-t( t(input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)
  
  ## estimate enrichment using top IP count bins
  ave.ip <- rowMeans(m6A)
  ave.top <- order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]

  geneCounts.window <- geneSum.norm[gene.name,]

  enrich <- as.data.frame(m6A[ave.top,]/geneCounts.window[ave.top,])
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(X)])

  norm.ip <-t( t(m6A)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)
  norm_lib <- c(readsOut, list('geneSum'=round(geneSum.norm),
                               'norm.input'=norm.input,
                               'norm.ip'=norm.ip,
                               'sizeFactor'=sizeFactor,
                               'geneBins'=geneBins,
                               'X'=X) 
                )
  
  par(mfrow=c(2,2))
  boxplot(log(geneSum[rowSums(geneSum)!=0,]+1),main = "INPUT")
  boxplot(log(geneSum.norm[rowSums(geneSum.norm)!=0,]+1),main = "Normalized INPUT")
  boxplot(log(enrich[rowSums(enrich)!=0,]+0.1), main = "Estimated enrichment")
  enrich.norm <- as.data.frame(norm.ip[ave.top,]/geneCounts.window[ave.top,])
  boxplot(log(enrich.norm[rowSums(enrich.norm)!=0,]+0.1), main = "Normalized estimated enrichment")
  par(mfrow=c(1,1))
  return(norm_lib)
}
