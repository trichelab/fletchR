
#' Add iterative LSI results to an SCE
#'
#' @param x               a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix       name of assay to use ("counts")
#' @param iter            number of LSI iterations (2)
#' @param preTfIdf        name of precomputed TF-IDF assay, or NULL
#' @param binarize        binarize counts before TF-IDF? (TRUE)
#' @param corCutOff       exclude depth-correlated dimensions during clustering (0.75)
#' @param excludeChr      chromosomes to exclude
#' @param nDimensions     number of LSI dimensions (30)
#' @param filterQuantile  quantile for filtering highly accessible features (0.995)
#' @param outlierQuantile quantiles for excluding depth outliers during LSI
#' @param varFeatures     number of variable features to select (25000)
#' @param totalFeatures   number of features considered for feature selection (500000)
#' @param depthCol        colData column holding sequencing depth ("nFrags")
#' @param outlierCellsCol name of logical colData column indicating outlier cells;
#'                        NULL uses quantile-based outlier detection
#' @param scaleTo         scale factor for normalization (10000)
#' @param method          matrix decomposition method ("SVD" or "NMF")
#' @param multimodal      logical indicating multimodal data (FALSE)
#' @param scaleMethod     modality scaling method ("Frobenius" or "minMax")
#' @param balanceFirstSelection logical to balance features in first-pass selection
#' @param balanceIterativeSelection logical to balance features in iterative selection
#' @param name            name of reducedDim to store the LSI embedding ("iterativeLSI")
#' @param seed            random seed (1)
#' @param ...             additional arguments passed to featureSelectionLSI
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
addIterativeLSI = function(x,
                            useMatrix = "counts",
                            iter = 2,
                            preTfIdf = NULL,
                            binarize = TRUE,
                            corCutOff = 0.75,
                            excludeChr = c("chrM", "chrX", "chrY"),
                            nDimensions = 30,
                            filterQuantile = 0.995,
                            outlierQuantile = c(0.02, 0.98),
                            varFeatures = 25000,
                            totalFeatures = 500000,
                            depthCol = "nFrags",
                            outlierCellsCol = NULL,
                            scaleTo = 10000,
                            method = "SVD",
                            multimodal = FALSE,
                            scaleMethod = "Frobenius",
                            balanceFirstSelection = FALSE,
                            balanceIterativeSelection = FALSE,
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
    corCutOff = corCutOff,
    excludeChr = excludeChr,
    nDimensions = nDimensions,
    filterQuantile = filterQuantile,
    outlierQuantile = outlierQuantile,
    varFeatures = varFeatures,
    totalFeatures = totalFeatures,
    depthCol = depthCol,
    outlierCellsCol = outlierCellsCol,
    scaleTo = scaleTo,
    method = method,
    multimodal = multimodal,
    scaleMethod = scaleMethod,
    balanceFirstSelection = balanceFirstSelection,
    balanceIterativeSelection = balanceIterativeSelection,
    seed = seed,
    ...
  )
  
  reducedDim(x, name) = LSIres
  
  message("Done.")
  return(x)
}
