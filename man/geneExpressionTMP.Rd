% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{geneExpressionTMP}
\alias{geneExpressionTMP}
\title{getTPM from geneSum}
\usage{
\S4method{geneExpressionTMP}{MeRIP.RADAR}(object,
  meanFragmentLength = 150, normalize = T)
}
\arguments{
\item{object}{The MeRIP.RADAR object}

\item{meanFragmentLength}{The mean length of RNA fragment (insert of RNA library). Default is 150bp.}

\item{normalize}{Logical indicating whether normalized TPM or raw TPM should be returned.}
}
\value{
A data.frame of gene quantification in Transcript Per Million (reads).
}
\description{
getTPM from geneSum
}
