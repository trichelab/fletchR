
# testing
bed <- "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_31726.bed.gz"
tbi <- paste(bed, "tbi", sep=".") 

# if a tabix index exists, the file is sorted:
stopifnot(file.exists(tbi))

library(rtracklayer)
system.time(gr <- import(bed))
#    user  system elapsed
#  53.852   1.048  55.275

# we can also consult the tabix headers to at least guess the genome used
library(Rsamtools)
hdr <- headerTabix(tbi) 
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
system.time(
  bed_data <- 
    read.table(bed, 
      sep="\t", header=FALSE,
      col.names=c("chrom", "chromStart", "chromEnd", "name", "score"),
      colClasses=c("character","integer","integer","character","numeric")
    )
)

# screw that, do it in parallel 
library(vroom)
bed5_col_types <- list(chrom = "c",
		      chromStart = "i",
		      chromEnd = "i",
    		      name = "c", 
	    	      score = "i")
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

library(GenomeInfoDb)
bed_genome <- "hg38" 
seqinfo <- SeqinfoForUCSCGenome(bed_genome)
shared <- intersect(seqlevels(seqinfo), hdr$seqnames)
tiles <- tileGenome(seqinfo[shared], tilewidth=50000, cut.last=TRUE)
system.time(tiles$frags <- countOverlaps(tiles, bed_gr2))
#    user  system elapsed 
#   2.890   0.237   3.142 


# obviously need to do this by cell barcode, which implies selection, so...

system.time(cells <- unique(bed_tbl$name))
#    user  system elapsed 
#   2.655   0.107   2.775 

length(cells)
# [1] 81729

# this would be a COUNT operation in duckdb obvs
frags_per_cell <- function(cell, frags) length(which(frags$name == cell))
system.time(fpc <- sapply(cells, frags_per_cell, frags=bed_tbl))
# this is slow AF 
# ^C
# Timing stopped at: 1609 47.08 1667

# write to parquet
library(nanoparquet) 
stub <- sub("\\.bed\\.gz$", "", bed) 
bed_pqt <- paste(stub, "parquet", sep=".")
system.time(write_parquet(bed_tbl, bed_pqt))
#    user  system elapsed 
#  11.648   0.166  12.081 

# open with duckdb
# which takes FOREVER to install btw 
# see https://duckdb.org/community_extensions/extensions/duckhts
# and https://cran.r-project.org/web/packages/Rduckhts/refman/Rduckhts.html
library(duckdb)
library(Rduckhts) 
con <- rduckhts_connect()
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

dbGetQuery(con=con, 
	   paste0("CREATE INDEX IF NOT EXISTS cell_idx ON ", stub, "(name)"))
# this takes a while fwiw, maybe 20 seconds 
(fpc <- dbGetQuery(con=con, 
 		   paste0("SELECT name AS cell, ",
	      		   "COUNT(*) AS frags ",
		           "FROM ", stub, " ",
			   "GROUP BY name")))
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

# see https://bwlewis.github.io/duckdb_and_r/ranges/ranges_redux.html
