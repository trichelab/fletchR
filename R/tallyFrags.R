#' tally fragments from a parquet table over bins into a sparse matrix
#' 
#' @param parquet     parquet filename holding fragments, or result of a chain
#' @param bins        a basepair size, BED file, GR, or REMO.hg38 (see Details) 
#' @param assayType   if bins is NULL, an assay (coCUTnTag or ATAC) for defaults
#' @param genome      if bins is NULL, a genome (default hg38) for binning
#' @param excludeChrs if bins is NULL, stick to standard chroms only? (TRUE) 
#' @param excludeList if bins is NULL, an excludelist to flag tainted fragments
#' @param standardChr drop nonstandard (ALT) contigs in addition to excludeChrs?
#' @param binSize     if bins and binSize are both NULL, default to 500 or 50000
#' @param binFile     optional filename for the bins parquet (default noted)
#' @param con         optional connection; will create a new one if needed
#' @param ...         optional arguments to rduckhts::connect()
#'
#' @return            a sparse dgCMatrix (or, eventually, a DuckDBArray)
#'
#' @details duckhts supports cgranges indices directly in the database, and 
#'          these indices can be created directly from a BED file, from a table,
#'          or persistent across a session (so if you're tallying a lot of 
#'          fragment files, it may be best to persist the index and name it). 
#'          The long and the short of present circumstances is that I'm trying
#'          to figure out the most elegant way to let users specify this while
#'          still having a non-infuriating default (500b for ATAC, 50kb for 
#'          coCUT&Tag) and also supporting arbitrary genomes and excludelists. 
#'          I also don't know whether a DuckDBArray is going to be a better 
#'          long-term choice for this task, or whether dgCMatrix is enough.
#'          'genome' can be pretty much anything that is known to 
#'          Seqinfo::registered_UCSC_genomes$genome or 
#'          Seqinfo::registered_NCBI_assemblies$assembly, as that is where the 
#'          chromosome names and lengths will be obtained for binning. Try 
#'          tileGenome(Seqinfo(genome="apiMel1"), tilewidth=500, cut.last=TRUE)
#'          or tileGenome(Seqinfo(genome="hg38"), tilewidth=5e4, cut.last=TRUE)
#'          for examples default binning schemes for ATAC and coCUTnTag resp.
#'          At the moment it seems easiest to use a DuckDBGRanges to set up the
#'          cgranges index for fragment filtering and counting. Also, the output
#'          may change to a DuckDBMatrix in the very near future. Last but not
#'          least, we might switch from zstd compression to snappy compression
#'          to trade file size for speed. At the moment, we use zstd everywhere.
#'
#' @import GenomicRanges
#' @import DuckDBArray
#' @import nanoparquet 
#' @import Rduckhts
#' @import Seqinfo
#' @import Matrix
#' @import DBI 
#'
#' @export
#'
tallyFrags <- function(parquet, bins=NULL, assayType=c("coCUTnTag", "ATAC"), genome="hg38", excludeChrs=c('chrX', 'chrY', 'chrM'), excludeList=NULL, binSize=NULL, binFile=NULL, standardChr=TRUE, con=NULL, ...) {

  stopifnot(file.exists(parquet))

  # refactor this out ?
  if (is.null(bins)) { # {{{
    availableGenomes <- c(registered_UCSC_genomes()$genome,
                          registered_NCBI_assemblies()$assembly)
    if (!genome %in% availableGenomes) {
      stop("Unrecognized genome. See Seqinfo::registered_UCSC_genomes().")
    }
    assayType <- match.arg(assayType)
    if (assayType == 'coCUTnTag') { 
      if (is.null(binSize)) binSize <- 5e4 
    } else if (assayType == 'ATAC') { 
      if (is.null(binSize)) binSize <- 500 
    }
    message("Tiling ", genome, " with ", binSize, "kb bins...")
    bins <- tileGenome(Seqinfo(genome=genome), tilewidth=binSize, cut.last=TRUE)
    if (standardChr) {
      message("Pruning to standard chromosomes...")
      bins <- keepStandardChromosomes(bins, pruning.mode="coarse") 
    }
    if (length(excludeChrs) > 0) {
      message("Dropping excluded chromosomes...") 
      bins <- dropSeqlevels(bins, excludeChrs, pruning.mode="coarse")
    }
  } # }}}

  # this may or may not be necessary
  pq_tmp <- tempdir()
  pq_tbl <- paste0("bins.", binSize, "bp.", genome)
  pq_bins <- paste0(pq_tbl, ".parquet")
  pq_path <- file.path(pq_tmp, pq_bins) 
  
  # helper fn
  bed3 <- function(gr) { 
    df <- as.data.frame(gr)[, 1:3]
    names(df)[1:3] <- c("chrom", "chromStart", "chromEnd")
    return(df)
  }

  message("Writing bins to ", pq_path, "...") 
  pq_scm <- parquet_schema(chrom="STRING", chromStart="INT64", chromEnd="INT64")
  write_parquet(bed3(bins), file=pq_path, schema=pq_scm)
  stopifnot(file.exists(pq_path))
  message("Done.")

  stop("Not finished yet :-(") 
  
  # start here (hopefully always...)
  newcon <- is.null(con)
  if (newcon) con <- rduckhts_connect(...)

  # create bin index 
  pq_tblname <- gsub("\\.", "_", pq_tbl)
  dbExecute(con, 
            paste0("CREATE OR REPLACE TABLE ", pq_tblname, 
                   " AS SELECT * FROM '", pq_path, "'"))

  pq_idx_name <- paste0(pq_tblname, "_tmp_idx")
  pq_tbl_colnames <- detectChromStartEnd(pq_path)
  names(pq_tbl_colnames) <- paste0("idx_", names(pq_tbl_colnames))
  dbGetQuery(con,
             with(pq_tbl_colnames, 
                  paste0("SELECT * FROM duckhts_cgranges_from_table(",
                         "'", pq_idx_name, "', ",
                         "'", pq_tblname, "', ",
                         "'", idx_chrom, "', ",
                         "'", idx_chromStart, "', ", 
                         "'", idx_chromEnd, "')")))

  pq_pqt_colnames <- detectChromStartEnd(parquet)
  names(pq_pqt_colnames) <- paste0("pqt_", names(pq_pqt_colnames))
  both_colnames <- c(pq_tbl_colnames, pq_pqt_colnames)

  chroms <- paste0("('", paste(seqlevels(bins), collapse="', '"), "')")

  sql <-with(both_colnames, 
             paste0("SELECT ",
                    pqt_name, " AS name, ",
                    pqt_chrom, " AS chrom, ",
                    pqt_chromStart, " AS chromStart, ",
                    pqt_chromEnd, " AS chromEnd, ",
                    "duckhts_cgranges_count_overlaps('", pq_idx_name, "', 
                                                     ", idx_chrom, ",
                                                     ", idx_chromStart, ", 
                                                     ", idx_chromEnd, ") ",
                   "AS frags ",
                   "FROM read_parquet('", parquet, "') ",
                   "WHERE chrom IN ", chroms,
                   "GROUP BY (name, chrom, chromStart, chromEnd) ",
                   "ORDER BY chrom, chromStart, chromEnd, name"))
  res <- dbGetQuery(con, sql)


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

  if (newcon) DBI::dbDisconnect(con, shutdown = TRUE)

}

