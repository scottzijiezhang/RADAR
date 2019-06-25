#' @import DESeq2
#' @export
setGeneric("counts", getGeneric("counts", package = "DESeq2"))

#' @export
#' @rdname counts
setMethod("counts", signature("MeRIP"), function(object){
  object@reads
})

#' @export
#' @rdname Input.files
setGeneric("Input.files", function(object) {
  standardGeneric("Input.files")
})

#' @export
#' @rdname IP.files
setGeneric("IP.files", function(object) {
  standardGeneric("IP.files")
})

#' @export
setGeneric("geneBins", function(object) {
  standardGeneric("geneBins")
})


#' @export
setGeneric("filterBins", function(object,  minCountsCutOff = 10 ) {
  standardGeneric("filterBins")
})

#' @export
setGeneric("extractInput", function(object, ... ) {
  standardGeneric("extractInput")
})

#' @export
setGeneric("extractIP", function(object, ...) {
  standardGeneric("extractIP")
})

#' @export
setGeneric("PrepCoveragePlot", function(object, ...) {
  standardGeneric("PrepCoveragePlot")
})

#' @export
setGeneric("normalizeLibrary",function(object,  boxPlot = TRUE){
  standardGeneric("normalizeLibrary")
})

#' @export
setGeneric("adjustExprLevel",function(object, adjustBy = "geneSum" ){
  standardGeneric("adjustExprLevel")
})

#' @export
setGeneric("geneExpression",function(object, ...){
  standardGeneric("geneExpression")
})


#' @export
setGeneric("variable",function(object){
  standardGeneric("variable")
})

#' @export
setGeneric("variable<-",function(object, value){
  standardGeneric("variable<-")
})

#' @export
setGeneric("samplenames",function(object){
  standardGeneric("samplenames")
})

#' @export
setGeneric("samplenames<-",function(object, value){
  standardGeneric("samplenames<-")
})


#' @export
setGeneric("peakDistribution",function(object){
  standardGeneric("peakDistribution")
})

#' @export
setGeneric("plotGeneCov", function(object, geneName, libraryType = "opposite", center = mean,ZoomIn = NULL, adjustExprLevel = FALSE , split = FALSE){ standardGeneric("plotGeneCov")})


#' @export
setGeneric("geneExpressionTMP",function(object, meanFragmentLength = 150, normalize = T){
  standardGeneric("geneExpressionTMP")
})

#' @export
setGeneric("diffIP",function(object,  exclude = NULL, maxPsi = 100, fdrBy = "fdr" ){
  standardGeneric("diffIP")
})

#' @export
setGeneric("diffIP_parallel",function(object,  exclude = NULL, maxPsi = 100, fdrBy = "fdr", thread){
  standardGeneric("diffIP_parallel")
})

#' @export
setGeneric("reportResult", function(object,  cutoff = 0.1, Beta_cutoff = 0.5, threads = 1) {
  standardGeneric("reportResult")
})



#' @export
setGeneric("select",  function(object, samples){standardGeneric("select")})


#' @export
setGeneric("results", function(object){standardGeneric("results")})

#' @export
setGeneric("plotHeatMap", function(object, covariates=TRUE ){standardGeneric("plotHeatMap")})





 