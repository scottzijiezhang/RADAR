#' @title gtfToGeneModel
#' @description to remove ambiguous gene model and return gene model as genomic ranges object
#' @param gtf gtf file to build gene model
#' @export
##helper function to remove ambiguous gene model and return gene model as genomic ranges object
gtfToGeneModel = function(gtf){

  #load gene models for  genes
  ## make the tx database object
  if(is.character(gtf)){
    txdb=makeTxDbFromGFF(gtf,format="gtf")
  }else{
    txdb = gtf
  }

  cols <- c("tx_chrom", "tx_strand")
  single_strand_genes <- genes(txdb, columns=cols)

  single_strand_genes_names = names(single_strand_genes)
  all_gene_names = names(genes(txdb,single.strand.genes.only=F))
  ambiguious_genes = setdiff(all_gene_names,single_strand_genes_names)
  ambi_genes = genes(txdb,single.strand.genes.only=F)[ambiguious_genes]
  resol_ambi_genes = lapply(ambi_genes,function(x){ return(x[1])  })

  #deal with ambiguous chromosome
  wantedChrom=unique(sapply(single_strand_genes$tx_chrom,unlist) )[nchar(unique(sapply(single_strand_genes$tx_chrom,unlist) )) <=5]
  id = which(sapply(single_strand_genes$tx_chrom,unlist) %in% wantedChrom )
  wantedGeneNames = names(single_strand_genes[id])

  geneModel = exonsBy(txdb,by="gene")

  single_strand_geneModel = geneModel[wantedGeneNames]
  ambi_geneModel = geneModel[ambiguious_genes]
  resol_ambi_geneModel = GRangesList()
  if(length(ambiguious_genes) > 0 ){
    for(i in 1: length(ambi_geneModel)){
    id = as.data.frame(findOverlaps(ambi_geneModel[[i]],resol_ambi_genes[names(ambi_geneModel)[i]][[1]]) )$queryHits
    newModel = GRangesList(ambi_geneModel[[i]][id])
    names(newModel) = names(ambi_geneModel)[i]
    if( all(match (as.character( seqnames(newModel[[1]]) ), wantedChrom,nomatch=0 )>0 ) ) {
      resol_ambi_geneModel = c(resol_ambi_geneModel,newModel)
      }
    }
  }


  final_geneModel = c(single_strand_geneModel,resol_ambi_geneModel)
  seqlevels(final_geneModel) = wantedChrom

  return(final_geneModel)
}
