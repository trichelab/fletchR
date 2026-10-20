
# testing
bed <- "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_31726.bed.gz"
checkTabix <- function(bed) { 
  tbi <- paste(bed, "tbi", sep=".") 
  stopifnot(file.exists(tbi))
}
# if a tabix index exists, the file is sorted:

library(rtracklayer)
system.time(gr <- import(bed))
#    user  system elapsed
#  53.852   1.048  55.275

# we can also consult the tabix headers to at least guess the genome used
library(Rsamtools)
tbi <- paste(bed, "tbi", sep=".") 
hdr <- headerTabix(bed) 
contigs <- hdr$seqnames
# there's a trick to obtain the chrom.sizes which I forget at the moment
# anyways

system.time(tbxbed <- TabixFile(bed))
#    user  system elapsed 
#   0.001   0.000   0.001 
#
# The above is obviously not the same as computing overlaps 

# current fragment BEDs are actually a weird BED5-ish format 
# e.g. 
# GL000008.2      71      369     ATGGGAAC_TTGCCTAA_H05-17        1
# GL000008.2      109     200     TGGATCTG_CAGCAACG_G05-17        2

# BED4 columns: chrom, chromStart, chromEnd, name
# BED5 columns: chrom, chromStart, chromEnd, name, score
if (FALSE) { 
  system.time(
    bed_data <- 
      read.table(bed, 
        sep="\t", header=FALSE,
        col.names=c("chrom", "chromStart", "chromEnd", "name", "score"),
        colClasses=c("character","integer","integer","character","numeric")
      )
  )
} # pointless now that we have duckhts

# screw that, do it in parallel 
library(vroom)
bed5_col_types <- list(chrom = "c",
		      chromStart = "i",
		      chromEnd = "i",
    		      name = "c", 
	    	      score = "i")
if (FALSE) { 
  system.time(
    bed_tbl <- 
      vroom::vroom(
       bed,
       delim = "\t",
       col_names = names(bed5_col_types), 
       col_types = bed5_col_types, 
       col_select = 1:4 # for now
       )
  )
  #    user  system elapsed                                                       
  #   8.953   3.702  12.038 
  # vroom can also automatically read from remote compressed files 
  # can we coerce this directly to a GRanges?
  # Of course:
  system.time(bed_gr2  <- as(bed_tbl, "GRanges"))

  #    user  system elapsed 
  #  27.290   1.511  28.957 
  #
  # observe that reading with vroom and then coercing to a GR is faster
  # than read.table or similar (!) 
}

library(Seqinfo)
library(GenomeInfoDb)
bed_genome <- "hg38" 
seqinfo <- Seqinfo(genome=bed_genome)
shared <- intersect(seqlevels(seqinfo), hdr$seqnames)
tiles <- tileGenome(seqinfo[shared], tilewidth=50000, cut.last=TRUE)
#system.time(tiles$frags <- countOverlaps(tiles, bed_gr2))
#    user  system elapsed 
#   2.890   0.237   3.142 


# obviously need to do this by cell barcode, which implies selection, so..
# system.time(cells <- unique(bed_tbl$name))
#    user  system elapsed 
#   2.655   0.107   2.775 

# length(cells)
# [1] 81729

# this would be a COUNT operation in duckdb obvs
# frags_per_cell <- function(cell, frags) length(which(frags$name == cell))
# system.time(fpc <- sapply(cells, frags_per_cell, frags=bed_tbl))
# this is slow AF 
# ^C
# Timing stopped at: 1609 47.08 1667

if (FALSE) { 
  # write to parquet
  library(nanoparquet) 
  stub <- sub("\\.bed\\.gz$", "", bed) 
  bed_pqt <- paste(stub, "parquet", sep=".")
  system.time(write_parquet(bed_tbl, bed_pqt))
  #    user  system elapsed 
  #  11.648   0.166  12.081 
}

# open with duckdb
# which takes FOREVER to install btw 
# see https://duckdb.org/community_extensions/extensions/duckhts
# and https://cran.r-project.org/web/packages/Rduckhts/refman/Rduckhts.html
library(duckdb)
library(Rduckhts) 
con <- rduckhts_connect()

if (FALSE) { 
  stub <- sub("\\.bed\\.gz$", "", bed)
  system.time(bed_hts <- rduckhts_bed(con, stub, bed)) # must be tabixed, duh
  #    user  system elapsed 
  #  35.688   2.741  38.677 

  bed_tbxpqt <- paste(stub, "tabix", "parquet", sep=".")
  system.time(rduckhts_tabix_convert_parquet(con, path=bed, output=bed_tbxpqt))
  #    user  system elapsed 
  #  27.141   0.285  27.658 

  dbGetQuery(con=con, "SHOW TABLES") 
  # data frame with 1 row and 1 column
  #                     name
  #              <character>
  # 	     1 HA_Hs_6L2BM_Rb_H3K27..
}

# Let's compute some overlaps. 
# It turns out we can tile the genome and then blacklist, save those blacklisted
# tiles as a BED, and then use that on the fly to query the tiles across BEDs.
#
tilesBed <- "tiles.50kb.nonBlacklisted.hg38.bed"
if (!file.exists(tilesBed)) {
  tiles <- tileGenome(seqinfo[shared], tilewidth=50000, cut.last=TRUE)
  encodeV2bb <- "https://hgdownload.soe.ucsc.edu/gbdb/hg38/problematic/encBlacklist.bb"
  encV2blacklist.hg38 <- import(encodeV2bb) # bigbed, don't look at me 
  badTiles <- queryHits(findOverlaps(tiles, encV2blacklist.hg38))
  goodTiles <- keepSeqlevels(tiles[-badTiles], paste0("chr", 1:22), 
                    			   pruning.mode="coarse")
  export(goodTiles, tilesBed)
  # tabix it internally or externally 
}
tilesTbx <- paste0(tilesBed, ".gz") 
stopifnot(file.exists(tilesTbx))

system.time(tiles_tbx <- rduckhts_tabix(con, "tiles", tilesTbx))
#    user  system elapsed 
#   0.062   0.004   0.134 

coCnT_beds <- list.files(patt="^HA.*bed.gz$") 
coCnT_tbx <- list.files(patt="^HA.*bed.gz.tbi$") 
overall <- length(coCnT_beds)
indexed <- length(coCnT_tbx)
unindexed <- overall - indexed
if (unindexed > 0) { 
  message(paste0(overall, " fragment BED files"))
  message(paste0(unindexed, " unindexed BED files"))
  message(paste0(indexed, " files have indices. Using those."))
  coCnT_beds <- sub(".tbi", "", coCnT_tbx)
} else { 
  message("All ", overall, " files indexed. Ingesting all fragment files.") 
}

stub <- "coCnT_Alex_MLL"
system.time(rduckhts_tabix_multi(con, 
                                 table=stub,
                                 files=coCnT_beds, 
                                 overwrite=TRUE))
# 
# takes about 10 minutes to ingest all of them  
#

dbGetQuery(con=con, 
	   paste0("CREATE INDEX IF NOT EXISTS cell_idx ON ", stub, "(name)"))
dbGetQuery(con=con, 
	   paste0("CREATE INDEX IF NOT EXISTS file_idx ON ", stub, "(filename)"))
# this takes a while fwiw, maybe 20 seconds 

(fpc <- dbGetQuery(con=con, 
 		   paste0("SELECT name AS cell, ",
	      		   "COUNT(*) AS frags ",
		           "FROM ", stub, " ",
	  		   "GROUP BY name, filename")))
# instantaneous
# data frame with 81729 rows and 2 columns
#                        cell     frags
#                 <character> <numeric>
# 		 1     AGGCAGAA_CAGCAACG_H0..     17760
# 		 2     TTGCTAAG_TCCATCAA_G0..      1924
# 		 3     CCTCGCAG_GCATTAAG_E0..      4546
# 		 4     GCGTTAAA_TGGAAATC_C0..     12365
# 		 5     TGGATCTG_ATCGAATG_B0..      5383
# 		 ...                      ...       ...
# 		 81725 CTCATGGG_TATTTGCG_F0..         1
# 		 81726 CGTACTAG_CGATAGGG_B0..         2
# 		 81727 TAAGGCGA_GGTGAAGG_F0..         1
# 		 81728 GTATTCGG_AGAGTAGA_F0..         1
# 		 81729 GGAGTAAG_GTAAGGAG_D0..         1
# 
# Somewhat insane: Rduckhts has a WASM hook already built-in

# see also https://github.com/Genentech/DuckDBGRanges for more 
# and https://github.com/Genentech/DuckDBDataFrame for more-er
# unfortunately both of these require arrow which takes FOREVER-ER to install

# Test out duckhts with remote files 
bigwig <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
gr <- as('chr22:20000000-20099999', "GRanges")
genome(gr) <- "hg38"

system.time(phylop_chunk_rt <- import(bigwig, which=gr))
#
# Timing stopped at: 9.237 31.13 1175

# hmm
bw_query <- 
  paste0("SELECT count(*) AS stored_intervals,         
                 round(min(VALUE)::DOUBLE, 3) AS minimum,
                 round(max(VALUE)::DOUBLE, 3) AS maximum
            FROM read_bigwig('", bigwig, "',
                             region := 'chr22:20000000-20099999')")
phylop_chunk <- dbGetQuery(con=con, bw_query) # problem 
# Error in `duckdb_result()`:
# ! Invalid Input Error: read_bigwig: failed to open a valid BigWig file
# ℹ Context: rapi_execute
# ℹ Error type: INVALID_INPUT
# Run `rlang::last_trace()` to see where the error occurred.

# New new new: figure out whether we can choose the top features beforehand
# so that we don't suffer from index starvation when doing LSI or NMF 
maxIdx <- 2147483647

# suppose we start with the top 500000 features for step 1 and 25000 overall
# let's suppose that we have 250000 cells but we want some headroom, so 5e5 
# to support more cells, let's chop the number of features in the first pass 
maxFeat <- 1e5 
maxCells <- 5e5 
( maxFeat * maxCells ) / maxIdx
#
# [1] 23.28
#
# so, we need to chop this down a bit before it hits R.  
# we end up ditching the top and bottom 2% of cells, so 
# 
(( maxFeat * maxCells ) * 0.96) / maxIdx
#
# [1] 22.35
#
# that's not going to do it.  How about subsampling?
#
subsampled <- 0.02
(( maxFeat * maxCells * subsampled ) * 0.96) / maxIdx
#
# [1] 0.447034836
#
# Now we're getting somewhere.  As long as we don't need more than about 
# 100K features in the first pass, we can get away with subsampling down 
# to 10K cells (repeatedly, one would hope) for even a fully dense array.
# 
# If our data is sparse, then of course we can do a lot better than this.
# But the conclusion here is that choosing features for the "full" decomp 
# based on a subsample of (up to 10000) cells can be good enough to determine
# which dimensions are being driven entirely by technical artifacts. 
#
# Once we get down to the "final" dimension reduction we can use more:
#
final
finalFeat <- 2.5e4
(( finalFeat * maxCells ) * 0.96) / maxIdx
# 
# [1] 5.58793545
#
# alright maybe not ALL of them but let's say
#
(( 1e4 * 2e5 ) * 0.96) / maxIdx 
# 
# [1] 0.894069672
# 
# And that should usually be enough dense features to work with. 
#
# Note also: https://bioconductor.org/packages/3.23/bioc/vignettes/ScaledMatrix/inst/doc/ScaledMatrix.html
# 

