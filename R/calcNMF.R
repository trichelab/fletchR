#' Core computational engine: Compute NMF on a matrix
#'
#' @param mat Feature-by-cell matrix (e.g. TF-IDF).
#' @param nDimensions Number of NMF factors.
#' @param outliers Logical vector indicating outlier cells
#'   to ignore during NMF fitting.
#' @param seed set seed
#'
#' @return A LinearEmbeddingMatrix containing NMF cell factors.
#'
calcNMF = function(mat,
                   nDimensions = 30,
                   outliers = NULL,
                   seed = 1) {
  if (is.null(outliers)) {
    outliers = rep(FALSE, ncol(mat))
  }
  
  message("Running NMF...")
  mat = mat[Matrix::rowSums(mat[, !outliers, drop = FALSE]) > 0, , drop = FALSE]
  nmf = RcppML::nmf(data = mat[, !outliers, drop = FALSE], k = nDimensions, seed = seed)
  
  matNMF = projectNMF(mat = mat, nmf = nmf)

  lem = SingleCellExperiment::LinearEmbeddingMatrix(
    sampleFactors = matNMF,
    featureLoadings = Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, ncol(matNMF)),
      giveCsparse = TRUE
    ),
    metadata = list(nmf = nmf)
  )
  
  dimnames(lem) = dimnames(matNMF)
  
  stopifnot(identical(
    dimnames(SingleCellExperiment::sampleFactors(lem)),
    dimnames(matNMF)
  ))
  
  stopifnot(identical(SingleCellExperiment::sampleFactors(lem), matNMF))
  
  return(lem)
}