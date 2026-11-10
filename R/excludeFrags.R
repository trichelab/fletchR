#' blacklist/excludelist "tainting" of fragments with overlaps 
#' 
#' @param parq    parquet filename holding fragments 
#' @param lst     excludelist (blacklist) of regions (BED, GR, or table)
#' @param con     optional connection; will create a new one if needed
#'
#' @return        a result without tainted fragments, for further processing
#'
#' @import DBI 
#' @import Rduckhts
#'
#' @export
#'
excludeFrags <- function(parq, lst, con=NULL) { 

  stop("Not finished yet") 

}
