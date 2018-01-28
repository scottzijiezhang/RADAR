\name{gtfToGeneModel}
\alias{gtfToGeneModel}
\title{gtfToGeneModel}
\usage{
Function to resolve ambigious gene annotations
}
\description{
Take gtf file and resolve ambigious gene models. Then return a genomic range list with gene models.
}
\examples{
gene_model = gtfToGeneModel("~/work/genome/hg38/hg38_UCSC.gtf")
}
