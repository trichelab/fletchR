#' simple wrapper to yield a DuckDBGRanges from a fletchR fragments parquet file
#'
#' @param parquet   the parquet filename 
#' @param ...       optional additional arguments to DuckDBGRanges
#' 
#' @details this function will attempt to figure out which columns in the 
#'          (presumably BED-derived) parquet file match chrom, start, end,
#'          and the 'name' mcol for cell barcodes or other naming information.
#'          If there is no 'name' column in the parquet file, there will be no
#'          name mcol in the returned object, to avoid annoying BED3 issues.
#'
#' @return a DuckDBGranges 
#'
#' @import nanoparquet 
#' @import DuckDBGRanges
#'
#' @export 
#' 
ddbGRangesFromParquet <- function(parquet, ...) { 

  cols <- nanoparquet::read_parquet_schema(parquet)$name  
  chrom <- intersect(c("chrom", "seqnames"), cols)[1]
  chromStart <- intersect(c("chromStart", "start"), cols)[1]
  chromEnd <- intersect(c("chromEnd", "end"), cols)[1]
  name <- intersect(c("name", "names"), cols)[1]

  if (is.na(name)) { 
    DuckDBGRanges(conn=parquet, 
                  seqnames=chrom, start=chromStart, end=chromEnd)
  } else { 
    DuckDBGRanges(conn=parquet, 
                  seqnames=chrom, start=chromStart, end=chromEnd, mcols=name)
  }

}
