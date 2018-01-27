#' @title normalizeLibrary
#' @param reads Read counts matrix from countReads() function.
#' @param samplenames A vector of Samples names. Needs to be the same as the one input to countReads()
#' @param X Study design or Grouping for the samples, should have 2 levels
#' @export norm_lib returns the list of normalized reads, study design and gene-bin names
normalizeLibrary <- function(reads,samplenames,X){
  input <- reads[,1:length(samplenames)]
  m6A <- reads[,(1+length(samplenames)):(2*length(filenames))]

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
  colnames(geneSum) <- filenames

  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)

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
  norm_lib <- list('geneSum'=geneSum.norm,
                  'norm.input'=norm.input,
                  'norm.ip'=norm.ip,
                  'sizeFactor'=sizeFactor,
                  'geneBins'=geneBins,
                  'X'=X)

  return(norm_lib)
}
