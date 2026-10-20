#' add iterative LSI results to an SCE
#'
#' @param x               a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix       name of assay to use ("counts")
#' @param iter            number of LSI iterations (2)
#' @param preTfIdf        name of precomputed TF-IDF assay, or NULL
#' @param binarize        binarize counts before TF-IDF? (TRUE)
#' @param corCutOff       exclude depth-correlated dimensions during clustering (0.75)
#' @param excludeChr      chromosomes to exclude
#' @param nDimensions     number of LSI dimensions (30)
#' @param outlierQuantile quantiles for excluding depth outliers during LSI
#' @param varFeatures     number of variable features to select (25000)
#' @param totalFeatures   number of features considered for feature selection (500000)
#' @param depthCol        colData column holding sequencing depth ("nFrags")
#' @param scaleTo         scale factor for normalization (10000)
#' @param name            name of reducedDim to store the LSI embedding ("iterativeLSI")
#' @param seed            random seed (1)
#' @param ...             additional arguments
#'
#' @return                SingleCellExperiment with reducedDim(x, name)
#'
#' @details               Runs featureSelectionLSI and stores the final LSI
#'                        embedding in reducedDim(x, name).
#'
#' @seealso               featureSelectionLSI
#' @seealso               ArchR::addIterativeLSI
#'
#' @import                SingleCellExperiment
#' @import                GenomeInfoDb
#' @import                irlba
#'
#' @export
addIterativeLSI <- function(x,
                            useMatrix = "counts",
                            iter = 2,
                            preTfIdf = NULL,
                            binarize = TRUE,
                            corCutOff = 0.75,
                            excludeChr = c("chrM", "chrX", "chrY"),
                            nDimensions = 30,
                            outlierQuantile = c(0.02, 0.98),
                            varFeatures = 25000,
                            totalFeatures = 500000,
                            depthCol = "nFrags",
                            scaleTo = 10000,
                            name = "iterativeLSI",
                            seed = 1,
                            ...) {
  stopifnot(depthCol %in% names(SummarizedExperiment::colData(x)))
  
  LSIres = featureSelectionLSI(
    x = x,
    useMatrix = useMatrix,
    iter = iter,
    preTfIdf = preTfIdf,
    binarize = binarize,
    corCutOff = 0.75,
    excludeChr = excludeChr,
    nDimensions = nDimensions,
    outlierQuantile = outlierQuantile,
    varFeatures = varFeatures,
    totalFeatures = totalFeatures,
    depthCol = depthCol,
    scaleTo = scaleTo,
    seed = seed
  )
  
  reducedDim(x, name) <- LSIres
  message("Done.")
  return(x)
}