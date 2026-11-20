#' detect columns in a parquet file corresponding to chrom, start, end[, name]
#' 
#' @param parquet     parquet filename
#'
#' @return            a list with at least three non-NULLs
#'
#' @import nanoparquet
#'
#' @export
detectChromStartEnd <- function(parquet) { 

  stopifnot(file.exists(parquet))
  pq_names <- nanoparquet::read_parquet_schema(parquet)$name
  pq_cols <- c(chrom = intersect(c("chrom", "seqnames"), pq_names)[1],
               chromStart = intersect(c("chromStart", "start"), pq_names)[1],
               chromEnd = intersect(c("chromEnd", "end"), pq_names)[1],
               name = intersect(c("name", "names"), pq_names)[1])
  if (any(is.na(pq_cols[1:3]))) {
    notfound <- paste(names(which(is.null(pq_cols))), collapse=", ")
    stop("Could not detect ", notfound, " in ", parquet)
  }
  return(as.list(pq_cols))

}
