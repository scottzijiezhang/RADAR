% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{normalizeLibrary,MeRIP-method}
\alias{normalizeLibrary,MeRIP-method}
\title{normalizeLibrary}
\usage{
\S4method{normalizeLibrary}{MeRIP}(object, boxPlot = TRUE)
}
\arguments{
\item{object}{MeRIP.RADAR object.}
}
\value{
returns a MeRIP.RADAR object
}
\description{
Normalized the input as RNA-seq data and normalize IP by enrichment. Specifically, we normalize ip libraries sizes so that the geometry mean of enrichment are the same.
}
