#' Core computational engine: Compute LSI on a matrix
#'
#' @param mat Feature-by-cell matrix (e.g. TF-IDF).
#' @param nDimensions Number of singular values to calculate.
#' @param outliers Logical vector indicating outlier cells to ignore during SVD fitting.
#' @param seed Set seed
#' @return A list containing cell-by-dimension projected matrix and svd result
#' @import SingleCellExperiment, irlba
calcLSI <- function(mat,
                    nDimensions = 30,
                    outliers = NULL,
                    seed = 1) {
  if (is.null(outliers)) {
    outliers <- rep(FALSE, ncol(mat))
  }
  message("Running SVD...")
  set.seed(seed)
  mat = mat[rowSums(mat[, !outliers]) > 0, ]
  svd <- irlba::irlba(mat[, !outliers, drop = FALSE], nDimensions, nDimensions)
  
  rownames(svd$u) = rownames(mat)
  rownames(svd$v) = colnames(mat)[!outliers]
  
  matSVD = projectSVD(mat, svd$u, svd$d, nDimensions)
  
  lem = SingleCellExperiment::LinearEmbeddingMatrix(
    sampleFactors = matSVD,
    featureLoadings = Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, ncol(matSVD)),
      giveCsparse = TRUE
    ),
    metadata = list(svd = svd)
  )
  
  dimnames(lem) = dimnames(matSVD)
  stopifnot(identical(
    dimnames(SingleCellExperiment::sampleFactors(lem)),
    dimnames(matSVD)
  ))
  stopifnot(identical(SingleCellExperiment::sampleFactors(lem), matSVD))
  return(lem)
}
