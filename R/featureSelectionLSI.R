#' perform iterative LSI feature selection
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
#' @param scaleTo         scale factor for normalization (10000)
#' @param method          matrix decomposition method (SVD,NMF)
#' @param seed            random seed (1)
#'
#' @return                list containing LSI, selected features, and TF-IDF result
#'
#' @details               Performs iterative LSI following the ArchR approach:
#'                        LSI, clustering, variable feature selection, and repeated LSI.
#'
#' @seealso               ArchR::addIterativeLSI
#'
#' @import                SingleCellExperiment
#' @import                GenomeInfoDb
#' @import                irlba
#'
#' @export
featureSelectionLSI = function(x,
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
                               scaleTo = 10000,
                               method = "SVD",
                               seed = 1) {
  stopifnot(depthCol %in% names(SummarizedExperiment::colData(x)))
  # set.seed(seed)
  LSIres = .firstPass(
    x = x,
    useMatrix = useMatrix,
    preTfIdf = preTfIdf,
    binarize = binarize,
    excludeChr = excludeChr,
    nDimensions = nDimensions,
    filterQuantile = filterQuantile,
    outlierQuantile = outlierQuantile,
    varFeatures = varFeatures,
    totalFeatures = totalFeatures,
    scaleTo = scaleTo,
    method = method,
    seed = seed
  )
  
  if (iter <= 1) {
    return(LSIres)
  }
  
  clusters = clusterSeurat(
    mat = removeDepthCorrelatedDims(
      projectedMat = ArchRRowZscores(SingleCellExperiment::sampleFactors(LSIres)),
      depth = colData(x)[[depthCol]],
      corCutOff = corCutOff)
    )
  
  prevLSI = LSIres
  
  for (i in seq_len(iter - 1)) {
    currentLSI = .iterativePass(
      x = x,
      clusters = clusters,
      useMatrix = useMatrix,
      preTfIdf = preTfIdf,
      binarize = binarize,
      excludeChr = excludeChr,
      nDimensions = nDimensions,
      varFeatures = varFeatures,
      totalFeatures = totalFeatures,
      scaleTo = scaleTo,
      outlierQuantile = outlierQuantile,
      method = method,
      seed = seed
    )
    
    clusters = clusterSeurat(
      mat = removeDepthCorrelatedDims(
        projectedMat = ArchRRowZscores(SingleCellExperiment::sampleFactors(currentLSI)),
        depth = colData(x)[[depthCol]],
        corCutOff = corCutOff)
    )
    
    if (length(unique(clusters)) == 1) {
      currentLSI = prevLSI
      break
    }
    prevLSI = currentLSI
  }
  finalLSI = currentLSI
  
  return(finalLSI)
}


#' perform the first pass of iterative LSI
#'
#' @param x               a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix       name of assay to use ("counts")
#' @param preTfIdf        name of precomputed TF-IDF assay, or NULL
#' @param binarize        binarize counts before TF-IDF? (TRUE)
#' @param excludeChr      chromosomes to exclude
#' @param nDimensions     number of LSI dimensions (30)
#' @param filterQuantile  quantile for filtering highly accessible features (0.995)
#' @param outlierQuantile quantiles for excluding depth outliers during LSI
#' @param varFeatures     number of variable features to select (25000)
#' @param totalFeatures   number of features considered for feature selection (500000)
#' @param depthCol        colData column holding sequencing depth ("nFrags")
#' @param scaleTo         scale factor for TF-IDF normalization (10000)
#' @param method          matrix decomposition method (SVD,NMF)
#' @param seed            random seed (1)
#'
#' @return                list containing LSI, selected features, and TF-IDF result
#'
#' @keywords internal
.firstPass = function(x,
                      useMatrix = "counts",
                      preTfIdf = NULL,
                      binarize = TRUE,
                      excludeChr = c("chrM", "chrX", "chrY"),
                      nDimensions = 30,
                      filterQuantile = 0.995,
                      outlierQuantile = c(0.02, 0.98),
                      varFeatures = 25000,
                      totalFeatures = 500000,
                      depthCol = "nFrags",
                      scaleTo = 10000,
                      method = "SVD",
                      seed = 1) {
  set.seed(seed)
  mat = filterAndGetMat(
    x = x,
    useMatrix = useMatrix,
    excludeChr = excludeChr,
    binarize = binarize
  )
  
  selectedFeatures = findTopFeatures(
    mat,
    filterQuantile = filterQuantile,
    varFeatures = varFeatures,
    totalFeatures = totalFeatures
  )
  if (!is.null(preTfIdf)) {
    rm(mat)
    gc()
    #errors could be handled more elegantly here for missing TfIdf slot
    tfidfRes = .getTfIdf(
      x = x,
      selectedFeatures = selectedFeatures,
      preTfIdf = preTfIdf,
      excludeChr = excludeChr
    )
  } else{
    mat = mat[rownames(mat) %in% selectedFeatures, ,drop = FALSE]
    tfidfRes <- calcTfIdf(
      mat,
      outlierQuantiles = outlierQuantile,
      excludeZeros = TRUE,
      scaleTo = scaleTo)
    rm(mat)
    gc()
  }
  
  outlierIdx = colnames(tfidfRes$tfidf) %in% tfidfRes$outliers
  if(method =="SVD"){
    lem = calcLSI(mat=tfidfRes$tfidf,nDimensions=nDimensions,outliers=outlierIdx,seed = seed)
  }else if(method == "NMF"){
    lem = calcNMF(mat=tfidfRes$tfidf,nDimensions=nDimensions,outliers=outlierIdx,seed = seed)
  }else{
    stop("Invalid Method. Choose SVD or NMF")
  }
  lem@metadata$selectedFeatures = selectedFeatures
  lem@metadata$tfidf = tfidfRes
  
  
  return(lem)
}

#' perform an iterative pass of LSI
#'
#' @param x               a SingleCellExperiment, usually from archRtoSCE
#' @param clusters        cluster assignments used for variable feature selection
#' @param useMatrix       name of assay to use ("counts")
#' @param preTfIdf        name of precomputed TF-IDF assay, or NULL
#' @param binarize        binarize counts before TF-IDF? (TRUE)
#' @param excludeChr      chromosomes to exclude
#' @param nDimensions     number of LSI dimensions (30)
#' @param outlierQuantile quantiles for excluding depth outliers during LSI
#' @param varFeatures     number of variable features to select (25000)
#' @param totalFeatures   number of features considered for feature selection (500000)
#' @param scaleTo         scale factor for normalization (10000)
#' @param method          matrix decomposition method (SVD,NMF)
#' @param seed            random seed (1)
#'
#' @return                list containing LSI, selected features, and TF-IDF result
#'
#' @keywords internal
.iterativePass = function(x,
                          clusters,
                          useMatrix = "counts",
                          preTfIdf = NULL,
                          binarize = TRUE,
                          excludeChr = c("chrM", "chrX", "chrY"),
                          nDimensions = 30,
                          outlierQuantile = c(0.02, 0.98),
                          varFeatures = 25000,
                          totalFeatures = 500000,
                          # depthCol = "nFrags",
                          scaleTo = 10000,
                          method = "SVD",
                          seed = 1) {
  set.seed(seed)
  mat = filterAndGetMat(
    x = x,
    useMatrix = useMatrix,
    excludeChr = excludeChr,
    binarize = binarize
  )
  
  selectedFeatures = findVariableFeaturesByCluster(
    mat = mat,
    clusters = clusters,
    nFeatures = varFeatures,
    # totalFeatures = totalFeatures,
    scaleTo = scaleTo
  )
  
  if (!is.null(preTfIdf)) {
    rm(mat)
    gc()
    #errors could be handled more elegantly here for missing TfIdf slot
    tfidfRes = .getTfIdf(
      x = x,
      selectedFeatures = selectedFeatures,
      preTfIdf = preTfIdf,
      excludeChr = excludeChr
    )
  } else{
    mat = mat[rownames(mat) %in% selectedFeatures, ,drop = FALSE]
    tfidfRes <- calcTfIdf(
      mat,
      outlierQuantiles = outlierQuantile,
      excludeZeros = TRUE,
      scaleTo = scaleTo
    )
    rm(mat)
    gc()
  }
  
  outlierIdx = colnames(tfidfRes$tfidf) %in% tfidfRes$outliers
  if(method =="SVD"){
    lem = calcLSI(mat=tfidfRes$tfidf,nDimensions=nDimensions,outliers=outlierIdx,seed = seed)
  }else if(method == "NMF"){
    lem = calcNMF(mat=tfidfRes$tfidf,nDimensions=nDimensions,outliers=outlierIdx,seed = seed)
  }else{
    stop("Invalid Method. Choose SVD or NMF")
  }
  
  lem@metadata$selectedFeatures = selectedFeatures
  lem@metadata$tfidf = tfidfRes
  
  
  return(lem)
}

#' retrieve a precomputed TF-IDF matrix
#'
#' @param x                a SingleCellExperiment containing precomputed TF-IDF values
#' @param selectedFeatures features to retain
#' @param preTfIdf         name of assay containing precomputed TF-IDF values ("TfIdf")
#' @param excludeChr       chromosomes to exclude
#'
#' @return                 list containing the TF-IDF matrix, IDF values, and outlier cells
#'
#' @keywords internal
.getTfIdf = function(x,
                     selectedFeatures,
                     preTfIdf = "TfIdf",
                     excludeChr = c("chrX", "chrY", "chrM")) {
  tfidfRes <- list(
    tfidf = filterAndGetMat(
      x = x,
      useMatrix = preTfIdf,
      features = selectedFeatures,
      excludeChr = excludeChr,
      binarize = FALSE,
      replaceZeros = FALSE
    ),
    idf = metadata(x)[[preTfIdf]][["idf"]][rownames(x) %in% selectedFeatures],
    outliers = metadata(x)[[preTfIdf]][["outliers"]]
  )
  return(tfidfRes)
}
# testCorrelation <- function(pCheck, pCheck2) {
#   lapply(seq_len(ncol(pCheck)), function(x) {
#     cor(pCheck[, x], pCheck2[, x])
#   })
# }

# if (!is.null(preTfIdf)) {
#   tfidfRes <- list(
#     tfidf = filterAndGetMat(
#       x = x,
#       useMatrix = preTfIdf,
#       features = selectedFeatures,
#       excludeChr = excludeChr,
#       binarize = FALSE,
#       replaceZeros = FALSE
#     ),
#     idf = metadata(x)[[preTfIdf]][["idf"]],
#     outliers = metadata(x)[[preTfIdf]][["outliers"]]
#   )
# } else {
#   tfidfRes <- calcTfIdf(
#     mat[rownames(mat) %in% selectedFeatures, ],
#     outlierQuantiles = outlierQuantile,
#     excludeZeros = TRUE,
#     scaleTo = scaleTo
#   )
# }