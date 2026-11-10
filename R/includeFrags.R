#' bound the scope of (perhaps already-blacklisted) fragments to tally
#' 
#' @param parq    parquet filename holding fragments 
#' @param lst     includelist (blacklist) of regions (BED, GR, or table)
#' @param con     optional connection; will create a new one if needed
#'
#' @return        a result limited to the includelist of regions for processing
#'
#' @import DBI 
#' @import Rduckhts
#'
#' @export
#'
includeFrags <- function(parq, lst, con=NULL) { 

  stop("Not finished yet") 

}

