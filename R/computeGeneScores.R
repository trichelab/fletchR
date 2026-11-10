#' compute gene scores from fragment files given a model 
#' 
#' @param parq            parquet filename holding fragments
#' @param genes           gene representation to feed the model  
#' @param geneModel       "exp(-abs(x)/5000) + exp(-1)" (same as ArchR default) 
#' @param upstream        c(1000, 100000) by default
#' @param downstream      c(1000, 100000) by default
#' @param geneBoundaries  use gene start/end rather than TSS (TRUE) 
#' @param tileSize        size of bin for computing gene scores (500) 
#' @param ceiling         4, by default
#' @param scaleFactor     5, by default
#' @param scaleTo         10000, by default
#' @param excludeChr      c("chrY", "chrM"), by default
#' @param blacklist       exclude these regions by default (encodeV2) 
#' @param con             optional connection; will create a new one if needed
#'
#' @return                a dense matrix of gene scores for each cell
#'
#' @import DuckDBArray
#' @import Rduckhts
#' @import Matrix
#' @import DBI 
#'
#' @export
#'
computeGeneScores <- function(parq, genes, geneModel, upstream, downstream, geneBoundaries, tileSize, ceiling, scaleFactor, scaleTo, excludeChr, blacklist, con=NULL) {

  stop("Not finished yet") 

  newcon <- is.null(con)
  if (newcon) con <- rduckhts_connect(...)

  # need to scope the matrix in terms of columns (cells) and rows (regions)
  # this should probably be farmed out to a separate table for handling 
  #
  # It may be a better idea to just use DuckDBArray for this step
  # but most likely we will need to populate the sparse array first, 
  # THEN save to a parquet file for DuckDBArray to handle. 
  #
  res <- dbGetQuery(con, 
                    "SELECT row_idx, col_idx, value FROM matrix_table")

  # Construct the sparse matrix
  sparse_mat <- sparseMatrix(
    i = res$row_idx,
    j = res$col_idx,
    x = res$value
  )


}

