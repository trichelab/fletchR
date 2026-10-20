#' add an LSI result to an SCE (NOTE: this assumes prefiltering was done!) 
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF") 
#' @param name          name of reducedDim to store the LSI embedding ("LSI") 
#' @param scaleDims     Z-scale the dimensions? (TRUE)
#' @param corCutOff     drop dimensions where cor(LSIdim, nFrag) > corCutOff
#' @param excludeChr    Chromosomes to exclude
#' @param nDimensions   number of dimensions for the SVD (30) 
#' @param depth         colData column holding depth for correlations ('nFrags')
#' @param subsetLSI     subset to mcols(x)$usedForLSI? (TRUE) 
#'
#' @return              SingleCellExperiment with reducedDim(x, name)
#'
#' @details             This function assumes that ArchR selected the features!
#'                      Here we do a simple, noniterative one-pass SVD on TFIDF
#'
#' @seealso             addTfIdf
#' @seealso             ArchR::addIterativeLSI
#'
#' @import              SingleCellExperiment
#' @import              GenomeInfoDb
#' @import              irlba
#'
#' @export
#'
addLSI <- function(x,
  useMatrix="TfIdf",
  name="LSI",
  scaleDims=FALSE,
  corCutOff=0.75,
  excludeChr=c("chrM","chrX","chrY"),
  nDimensions=30,
  depth="nFrags",
  seed = 1,
  subsetLSI = FALSE,...) {
  set.seed(seed)
  stopifnot(depth %in% names(colData(x)))
 
  mat <- filterAndGetMat(x=x,useMatrix=useMatrix,excludeChr=excludeChr,subsetLSI=subsetLSI,prune=c(1,1),replaceZeros = FALSE)
  # message("Running SVD...")

  outliers = getOutliersIdx(x,useMatrix=useMatrix)
  lem <- calcLSI(mat=mat, outliers=outliers, nDimensions=nDimensions)
  
  reducedDim(x, name) <- lem
  message("Done.")
  return(x) 
}

#' Get the outlier cells from a SingleCellExperiment object
#'
#' @param x A SingleCellExperiment object.
#' @param useMatrix The name of the assay to use for outlier detection (default is "TfIdf").
#' @return A logical vector indicating which cells are outliers.
getOutliersIdx <- function(x, useMatrix="TfIdf") {
  outliers = metadata(x)[[useMatrix]][["outliers"]] # returns NULL if doesn't exist
  if (is.null(outliers)) {
    message("Either No outliers were detected or metadata(x)[[useMatrix]][['outliers'] does not exist. Confirm prior to exceeding")
  }
  return(colnames(x) %in% outliers)
}

