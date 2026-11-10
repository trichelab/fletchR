#' convert a fragment file to a parquet file, optionally restricting its scope
#' 
#' @param bed           the fragment (tabixed TSV or BED) file (must be indexed)
#' @param con           opt. duckhts-enabled duckdb connection (default: new)
#' @param parquet       opt. parquet file name (default: s/.bed.gz/.parquet/)
#' @param region        opt. region to load fragments (default: no restriction)
#' @param ...           opt. args to pass to rduckhts_connect() for new con
#' 
#' @return              TRUE on success, FALSE on failure. 
#'
#' @details Can be used as a DuckDBGRanges backend via ddbGRangesFromParquet.
#' Only the first four columns (chrom, chromStart, chromEnd, name) are imported.
#' region is expected to be a single string, e.g. "chr8:11808264-11808510"
#' although it sure would be handy to have multiple strings, e.g. for black-
#' listed standard chromosomes in the default case 
#' see https://duckdb.org/community_extensions/extensions/duckhts
#' and https://cran.r-project.org/web/packages/Rduckhts/refman/Rduckhts.html
#' We might start using DuckDBGRanges for region filtering at some point. 
#' 
#' @import Rduckhts
#' @import DBI
#'
#' @export
#'
fragmentsToParquet <- function(bed, con=NULL, parquet=NULL, region=NULL, ...) { 

  tabixed <- (file.exists(paste0(bed, ".tbi"))|file.exists(paste0(bed, ".csi")))
  if (!tabixed) stop("Your fragment (BED) file must be tabix-indexed.") 
  if (is.null(parquet)) parquet <- sub("\\.(bed|tsv)(.gz)?", ".parquet", bed)

  if (file.exists(parquet)) {

    message(parquet, " already exists! Skipping.")
    return(TRUE)

  } else { 

    newcon <- is.null(con)
    if (newcon) con <- rduckhts_connect(...)
    rduckhts_tabix_convert_parquet(
      con = con,
      path = bed,
      output = parquet, 
      header_names = c("chrom",
                       "chromStart",
                       "chromEnd",
                       "name"),
      column_types = c("VARCHAR",
                       "BIGINT",
                       "BIGINT",
                       "VARCHAR"),
      compression = "zstd"
    )
    if (newcon) DBI::dbDisconnect(con, shutdown = TRUE)

    # check
    if (file.exists(parquet)) {
      message(parquet, " written to disk.")
      return(TRUE)
    } else {
      message("ERROR: ", parquet, " was not written to disk.")
      return(FALSE)
    }

  }

}
