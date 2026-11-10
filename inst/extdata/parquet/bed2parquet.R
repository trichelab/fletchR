# for testing on various nodes and laptops
hostname <- function() Sys.getenv("HOSTNAME") 

# on HPC, we need `sinfo` for accurate results: 
# 
getSinfo <- function() { 
  res <- system('sinfo -o "%P %n %m %c"', intern=TRUE)
  host <- grep(Sys.getenv("HOSTNAME"), res, value=TRUE)
  mem_cpu <- sapply(strsplit(host, " "), "[", 3:4)
  rownames(mem_cpu) <- c("RAM", "CPU")
  resources <- apply(mem_cpu, 1, unique)
  sapply(resources, as.numeric)
}

# testing
#
# big one
# bed <- "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_31726.bed.gz"
#
# local one
bed<-"HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026.bed.gz"

checkTabix <- function(bed) file.exists(paste(bed, "tbi", sep="."))
# if a tabix index exists, the file is sorted
stopifnot(checkTabix(bed))

# nb. we could also use rduckhts_tabix_index(bed) to create indices as needed.


# test whether to go bed-by-bed or a multiview table (probably the former)
bed_files <- 
  read.csv("bed_files_metadata.csv")[,c("filename", "antibody_a", "antibody_b")]
bed_files$tabixed <- sapply(bed_files$filename, checkTabix)
stopifnot(all(bed_files$tabixed))

bf <- subset(bed_files, tabixed) 

# we can consult the tabix headers to at least guess the genome used
library(Rsamtools)
tbi <- paste(bed, "tbi", sep=".") 
hdr <- headerTabix(bed) 
contigs <- hdr$seqnames

# stock up one of the tables 
getContigs <- function(bed) headerTabix(bed)$seqnames
contigsSeen <- sort(Reduce(base::union, sapply(bf$filename, getContigs)))

# standard human autosomes are contigsSeen[ 1:22 ]
# chrM is usually contigsSeen[ 23 ] 
# chrX is usually contigsSeen[ 24 ] 
# chrY is usually contigsSeen[ 25 ] 
grep("chr", contigsSeen, value=T)
#  [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
# [10] "chr18" "chr19" "chr2"  "chr20" "chr21" "chr22" "chr3"  "chr4"  "chr5"
# [19] "chr6"  "chr7"  "chr8"  "chr9"  "chrM"  "chrX"  "chrY"

# there's a trick to obtain the chrom.sizes which I forget at the moment
# that trick is the one that allows us to infer which genome it came from

library(Rsamtools)
library(rtracklayer)
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
#
# The BED5-ish format records the number of PCR dupes in the score column
#
library(Seqinfo)
library(GenomeInfoDb)
bed_genome <- "hg38" 
seqinfo <- Seqinfo(genome=bed_genome)
hdr <- headerTabix(bed) 
contigs <- hdr$seqnames

# use full catalog instead
shared <- intersect(seqlevels(seqinfo), contigsSeen)

# drop chrM, chrX, chrY 
shared <- setdiff(shared, c("chrM", "chrX", "chrY")) 

tiles_500b <- tileGenome(seqinfo[shared], tilewidth=500, cut.last=TRUE)
tiles_5kb <- tileGenome(seqinfo[shared], tilewidth=5000, cut.last=TRUE)
tiles_50kb <- tileGenome(seqinfo[shared], tilewidth=50000, cut.last=TRUE)

#
#system.time(tiles$frags <- countOverlaps(tiles, bed_gr2))
#    user  system elapsed 
#   2.890   0.237   3.142 
#
tiles <- GRangesList(
  "500b"=tiles_500b, 
  "5kb"=tiles_5kb, 
  "50kb"=tiles_50kb
)
seqlevels(tiles)
#  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
# [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
# [19] "chr19" "chr20" "chr21" "chr22"

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
con <- rduckhts_connect(dbdir = "fletchR.duckdb", read_only = FALSE)
dbExecute(con=con, "SET memory_limit = '16GB'")  # for Alex, mostly 

#
# let's use contigsSeen for keying:
#
system.time(dbExecute(con=con,
                       paste0(" CREATE OR REPLACE TABLE contigs (", 
                              "   contig VARCHAR PRIMARY KEY, ",
                              "   standard BOOLEAN NOT NULL DEFAULT FALSE ", 
                              ")")))
#
#    user  system elapsed 
#   0.001   0.001   0.005 
#

dbGetQuery(con, "DESCRIBE contigs")[, 1:5]
# data frame with 2 rows and 6 columns
#   column_name column_type        null         key              default
#   <character> <character> <character> <character>          <character>
# 1      contig     VARCHAR          NO         PRI                   NA
# 2    standard     BOOLEAN          NO          NA CAST('f' AS BOOLEAN)

contigsSeen <- c(seqlevels(tiles), setdiff(contigsSeen, seqlevels(tiles)))
contigsStandard <- c(rep(TRUE, length(seqlevels(tiles))),
                     rep(FALSE, length(contigsSeen) - length(seqlevels(tiles))))
contigsOrdered <- data.frame(contig=contigsSeen, 
                             standard=contigsStandard)
# 
#    user  system elapsed 
#   0.002   0.001   0.004 
# 
dbWriteTable(con, "contigs", contigsOrdered, append=TRUE)

dbGetQuery(con, "SELECT * FROM contigs WHERE standard = TRUE")
# data frame with 22 rows and 2 columns
#          contig  standard
#     <character> <logical>
# 1          chr1      TRUE
# ...         ...       ...
# 22        chr22      TRUE

dbExecute(con=con,
          paste0("CREATE OR REPLACE TYPE contig ",
                 "AS ENUM (SELECT contig FROM contigs)"))

head(dbGetQuery(con, 
                "SELECT enum_range(NULL::contig) 
                     AS contigs")[1,1][[1]])
# [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"

tail(dbGetQuery(con, 
                "SELECT enum_range(NULL::contig) 
                     AS contigs")[1,1][[1]])
# [1] "KI270752.1" "KI270753.1" "KI270754.1" "KI270755.1" "KI270756.1"
# [6] "KI270757.1"

if (FALSE) { # {{{
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
} # }}}

# Can we smash it down?  Probably:
system.time(dbExecute(con=con,
                       paste0("CREATE OR REPLACE TABLE fragment_files (", 
                              " filename VARCHAR PRIMARY KEY,", 
                              " antibody_a VARCHAR NOT NULL DEFAULT 'unknown',",
                              " antibody_b VARCHAR NOT NULL DEFAULT 'unknown',",
                              " tabixed BOOLEAN NOT NULL DEFAULT FALSE", 
                              ")")))
dbGetQuery(con, "SHOW TABLES")
# data frame with 3 rows and 1 column
#             name
#      <character>
# 1  cell_barcodes
# 2        contigs
# 3 fragment_files

dbGetQuery(con, "DESCRIBE fragment_files")[, 1:5]
#
# data frame with 4 rows and 5 columns
#   column_name column_type        null         key              default
#   <character> <character> <character> <character>          <character>
# 1    filename     VARCHAR          NO         PRI                   NA
# 2  antibody_a     VARCHAR          NO          NA            'unknown'
# 3  antibody_b     VARCHAR          NO          NA            'unknown'
# 4     tabixed     BOOLEAN          NO          NA CAST('f' AS BOOLEAN)
#
dbWriteTable(con, "fragment_files", bed_files, append=TRUE)
dbExecute(con=con,
          paste0("CREATE OR REPLACE TYPE fragment_file ",
                 "AS ENUM (SELECT filename FROM fragment_files)"))

head(dbGetQuery(con, 
                "SELECT enum_range(NULL::fragment_file) 
                     AS fragment_files")[1,1][[1]])
# [1] "AE_Hs_5Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_112525.bed.gz"
# [2] "AE_Hs_5Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_112625.bed.gz"
# [3] "AE_Hs_5Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_112725.bed.gz"
# [4] "AE_Hs_5Leuk_Rb_DEM_PA1241_Ms_H3K27me3_H19_112525.bed.gz"
# [5] "AE_Hs_5Leuk_Rb_DEM_PA1241_Ms_H3K27me3_H19_112625.bed.gz"
# [6] "AE_Hs_5Leuk_Rb_DEM_PA1241_Ms_H3K27me3_H19_112725.bed.gz"

tail(dbGetQuery(con, 
                "SELECT enum_range(NULL::fragment_file) 
                     AS fragment_files")[1,1][[1]])
# [1] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt1_31926.bed.gz"
# [2] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt1_32026.bed.gz"
# [3] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_31726.bed.gz"
# [4] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_31826.bed.gz"
# [5] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_31926.bed.gz"
# [6] "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026.bed.gz"

# this one is trickier, may need to do a second pass to get it 
system.time(dbExecute(con=con,
                       paste0(" CREATE OR REPLACE TABLE cell_barcodes (", 
                              "   barcode VARCHAR PRIMARY KEY ",
                              ")")))

# Let's try loading fragment files one by one, 
# squashing the key-able columns into foreign keys,
# and then computing overlaps against tiles. 
#
bedfile <- bf[1, 1] # just for a test; it should be big enough
stub <- sub("\\.bed\\.gz$", "", bedfile)

# FDFT1 DMR: chr8:11808264-11808510
FDFT1_region <- "chr8:11808264-11808510"
FDFT1_tiles <- GRangesList(lapply(tiles, 
                                  subsetByOverlaps, 
                                  ranges=as(FDFT1_region, "GRanges")))

# regional import (for overlap computation)
# note that we may need to import BED files to individual parquet files
# then combine them all to compress out the redundant information
# (contigs, filenames, cell_barcodes) 
system.time(rduckhts_tabix(con, 
                           table_name = paste0(stub, "_FDFT1"), 
                           path = bedfile,
                           region = FDFT1_region, 
                           header_names = c("chrom",
                                            "start",
                                            "end",
                                            "name",
                                            "score"),
                           column_types = c("contig",
                                            "BIGINT",
                                            "BIGINT",
                                            "VARCHAR", 
                                            "BIGINT")))

#    user  system elapsed 
#   0.057   0.004   0.074 

# unrestricted: 

system.time(rduckhts_tabix(con, 
                           table_name = stub,
                           path = bedfile,
                           header_names = c("chrom",
                                            "start",
                                            "end",
                                            "name",
                                            "score"),
                           column_types = c("contig",
                                            "BIGINT",
                                            "BIGINT",
                                            "VARCHAR", 
                                            "BIGINT"),
                           overwrite = TRUE))
#    user  system elapsed 
#   1.947   0.061   2.043 
#

system.time(dbExecute(con=con,             
                      paste0("CREATE INDEX IF NOT EXISTS cell_idx ON ",
                              stub, "(name)")))
#    user  system elapsed 
#   6.169   0.074   1.908 

tail(dbGetQuery(con=con, paste("SELECT name, COUNT(*) AS fragments",
                               "FROM ", stub,
                               "GROUP BY name",
                               "HAVING fragments >= 100", 
                               "ORDER BY fragments DESC")))
head(dbGetQuery(con=con, paste("SELECT name, COUNT(*) AS fragments",
                               "FROM ", stub,
                               "GROUP BY name",
                               "HAVING fragments >= 100", 
                               "ORDER BY fragments DESC")))

dbGetQuery(con=con, "SHOW TABLES") 

library(rtracklayer)
names(tiles) <- paste("tiles", unique(genome(tiles)), names(tiles), sep=".")
for (tn in names(tiles)) { 
  export.bed(object=tiles[[tn]], con=paste0(tn, ".bed"), index=TRUE)
}

system.time(rduckhts_tabix(con, 
                           table_name = "tiles_50kb_hg38", 
                           path = "tiles.hg38.50kb.bed.bgz",
                           header_names = c("chrom",
                                            "start",
                                            "end",
                                            "barcode",
                                            "dupes"),
                           column_types = c("VARCHAR",
                                            "BIGINT",
                                            "BIGINT",
                                            "NULL",
                                            "NULL"),
                           overwrite = TRUE))
dbExecute(con, "ALTER TABLE tiles_50kb_hg38 DROP name")
dbExecute(con, "ALTER TABLE tiles_50kb_hg38 DROP score")
dbExecute(con, "UPDATE tiles_50kb_hg38 SET start = start + 1")
dbGetQuery(con, "SELECT * FROM tiles_50kb_hg38 LIMIT 6") 
#
#   chrom  start    end
# 1  chr1      1  50000
# 2  chr1  50001 100000
# 3  chr1 100001 150000
# 4  chr1 150001 200000
# 5  chr1 200001 250000
# 6  chr1 250001 300000
#

#### THIS IS THE IMPORTANT STUFF 
####
#### fragmentsToParquet already handles ingestion & indexing
####
#### optimizing blacklisted counting is the real trick here 
####

# convert all non-converted fragment files to parquet:
#
library(fletchR) 
packageVersion("fletchR")
# [1] '0.50.2'

# AE_Hs_5Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_112525.parquet written to disk.
# ...
# this happens at varying speeds depending on the file size 
# ...
checkParquet <- function(bed) {
  file.exists(sub("\\.(bed|tsv)(\\.gz)?", ".parquet", bed))
}
bf$parquet <- sapply(bf$filename, checkParquet)
sapply(subset(bf, !parquet)$filename, fragmentsToParquet)
bf$parquet <- sapply(bf$filename, checkParquet)
stopifnot(all(bf$parquet))

if (FALSE) { 
  # index this 
  dbGetQuery(con, "SELECT duckhts_cgranges_create('tiles_50kb_hg38_idx')")
  #   duckhts_cgranges_create('tiles_50kb_hg38_idx')
  #1                                           TRUE

  dbGetQuery(con, "SELECT duckhts_cgranges_destroy('tiles_50kb_hg38_idx')")
  #   duckhts_cgranges_destroy('tiles_50kb_hg38_idx')
  #1                                           TRUE
}

dbGetQuery(con, 
           paste0("SELECT * FROM duckhts_cgranges_from_table(",
                  "'tiles_50kb_hg38_idx', 'tiles_50kb_hg38', ",
                  "'chrom', 'start', 'end')"))
#   indexed
# 1    TRUE

dbGetQuery(con, 
           paste0("SELECT
                   chrom, start, end,
                   duckhts_cgranges_count_overlaps('tiles_50kb_hg38_idx', 
                                                   chrom, start, end) 
                   AS frags
                   FROM read_bed('HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026.bed.gz')"))
# failure

# let's switch to DuckDBGRanges
library(DuckDBGRanges) 
library(Rduckhts)
stub <- "HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026"
rduckhts_tabix_convert_parquet(
  con = con,
  path = paste0(stub, ".bed.gz"), 
  output = paste0(stub, ".output.parquet"), 
  header_names = c("chrom",
                   "start",
                   "end",
                   "name",
                   "score"),
  column_types = c("VARCHAR",
                   "BIGINT",
                   "BIGINT",
                   "VARCHAR",
                   "BIGINT"),
  compression = "zstd"
)

ddb_gr <- DuckDBGRanges(paste0(stub, ".output.parquet"), 
                        seqnames="chrom", start="start", end="end", 
                        mcols=c("name", "score"))

unique(ddb_gr$name)
# DuckDBColumn of length 77002
#                          1                          2 
# "ACAAAGTG_TGCATACT_G05-20" "CGAGGCTG_ACGATAAC_E05-20" 
#                          3                          4 
# "TAGGCATG_ACTGCATA_H05-20" "AAACGGTC_CTTTAACA_C06-20" 
#                          5                        ... 
# "TTGCTAAG_GATCATTC_G05-20"                        ... 

system.time(head(rev(sort(table(ddb_gr$name)))))
# 
#    user  system elapsed 
#   0.550   0.329   0.507 

head(rev(sort(table(ddb_gr$name))))
# name                 
# AAGAGGCA_GCAAATTT_H05-20 ACCGGCTA_GCATTAAG_A06-20 CAGAGAGG_TGCATACT_G05-20 
#                    35690                    11664                     7186 
# TGGATCTG_GATCATTC_G05-20 CAATCTTA_TAACCAAG_G05-20 CCTCGCAG_AAGGAGTA_D05-20 
#                     6284                     5891                     5737 

# 
# We can drop fragments that overlap blacklists:
#
library(rtracklayer) 
encodeV2bb <- 
  "https://hgdownload.soe.ucsc.edu/gbdb/hg38/problematic/encBlacklist.bb"
encV2blacklist.hg38 <- import(encodeV2bb) # bigbed, don't look at me 
encV2blacklist.hg38$score <- as.numeric(factor(encV2blacklist.hg38$name))
encV2blacklist.hg38$name <- NULL
export.bedGraph(encV2blacklist.hg38, 
                con="encodeV2blacklist.hg38.bed", 
                index=TRUE) 

system.time(
  rduckhts_tabix_convert_parquet(
    con = con,                                        
    path = "encodeV2blacklist.hg38.bed.bgz",
    output = "encodeV2blacklist.hg38.parquet",
    header_names = c("chrom",
                     "start",
                     "end"),
    column_types = c("VARCHAR",
                     "BIGINT",
                     "BIGINT"),
    compression = "zstd"
  )
)
#    user  system elapsed 
#   0.062   0.004   0.134 

blacklist_gr <- DuckDBGRanges("encodeV2blacklist.hg38.parquet", 
                              seqnames="chrom", start="start", end="end")

system.time(excludefrags <- queryHits(findOverlaps(ddb_gr, blacklist_gr)))
# 
#    user  system elapsed 
#   0.536   0.010   0.547 

system.time(includefrags <- setdiff(seq_along(ddb_gr), excludefrags))
# 
#    user  system elapsed 
#   0.041   0.011   0.052 

# Take inventory before doing conversions
#
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
  message("All ", overall, " files indexed.")
}

# stub <- "coCnT_Alex_MLL"
# system.time(rduckhts_tabix_multi(con, 
#                                  table=stub,
#                                  files=coCnT_beds, 
#                                  overwrite=TRUE))
# 
# takes about 832 seconds to ingest all of them on a shortq node
# and about 69GB of resident memory, 
# plus another 220GB in a blackhole somewhere
#
# system.time(dbGetQuery(con=con, 
#                        paste0("CREATE INDEX IF NOT EXISTS cell_idx ON ", 
#                               stub, "(name)")))
#
#     user   system  elapsed 
# 3729.071  122.319  991.213 
#
# system.time(dbGetQuery(con=con,
#                        paste0("CREATE INDEX IF NOT EXISTS file_idx ON ", 
#                               stub, "(filename)")))
#
#    user  system elapsed 
# 715.548  36.318 121.214 
# 

# Can we smash it down?  Probably:
system.time(dbExecute(con=con,
                       paste0(" CREATE OR REPLACE TABLE cell_barcodes (", 
                              "   barcode VARCHAR PRIMARY KEY ",
                              ")")))
system.time(dbExecute(con=con,
                       paste0(" CREATE OR REPLACE TABLE fragment_files (", 
                              "   filename VARCHAR PRIMARY KEY, ", 
                              "   antibody1 VARCHAR DEFAULT 'unknown', ", 
                              "   antibody2 VARCHAR DEFAULT 'unknown' ", 
                              ")")))



# so.  Alex's point about BED-by-BED is likely the way to go, 
# BUT we need to maintain the barcode, fragment file, and contig tables
# if we want these per-fragment-file tables and parquets to be small 


system.time(fpc <- dbGetQuery(con=con, 
                              paste0("SELECT column3 AS cell, ",
                                     "       filename AS fragment_file, ", 
                                     "       COUNT(*) AS frags ",
                                     "FROM ", stub, " ",
                                     "GROUP BY cell, fragment_file")))
#
#    user  system elapsed 
# 290.245  70.087 118.720 
# 
# a couple minutes for all cells x all fragment files 
# 
# not quite as instantaneous as it was per-fragment-file, 
# which supports Alex's scheme of per-fragment-file tables 
# materialized across an experiment/collection by a VIEW
# 

#
# can we shrink these down from the moby fragment table by splitting into 
# smaller per-fragment-file tables that we save as parquet tables or as 
# a full database (minus the original moby fragment table) in a .db file?
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


# Fun with overlap queries and filtration:
rduckhts_tabix(con,
               table_name = "tiles_50kb_hg38",
               path = "tiles.hg38.50kb.bed.bgz",
               header_names = c("chrom",
                                "start",
                                "end"),
               column_types = c("contig",
                                "BIGINT",
                                "BIGINT"),
               overwrite = TRUE)

dbGetQuery(con,
           paste0("SELECT
                   name, chrom, start, \"end\",
                   duckhts_cgranges_count_overlaps('tiles_50kb_hg38_idx', 
                                                   chrom, start, \"end\") 
                   AS frags
                   FROM read_parquet('HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026.parquet')
                   WHERE chrom IN ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
                                   'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                   'chr21', 'chr22')
                   GROUP BY (name, 
                            chrom, 
                            start, 
                            \"end\")
                   ORDER BY frags DESC"))


dbGetQuery(con, 
           paste0("SELECT * FROM duckhts_cgranges_from_table(",
                  "'tiles_50kb_hg38_idx', 'tiles_50kb_hg38', ",
                  "'chrom', 'start', 'end')"))
#   indexed
# 1    TRUE
# time it
system.time(res <- 
            dbGetQuery(con,
           paste0("SELECT
                   name, chrom, start, \"end\",
                   duckhts_cgranges_count_overlaps('tiles_50kb_hg38_idx', 
                                                   chrom, start, \"end\") 
                   AS frags
                   FROM read_parquet('HA_Hs_6L2BM_Rb_H3K27me3_H1_wMs_H3K4me2_H13_12xS_12xM_CoCnt2_32026.parquet')
                   WHERE chrom IN ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
                                   'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                   'chr21', 'chr22')
                   GROUP BY (name, 
                            chrom, 
                            start, 
                            \"end\")
                   ORDER BY chrom, start, \"end\", name"))
            )
#    user  system elapsed 
#   2.995   2.407   1.283 


# same but for a very large file
system.time(res2 <- 
            dbGetQuery(con,
           paste0("SELECT
                   name, chrom, start, \"end\",
                   duckhts_cgranges_count_overlaps('tiles_50kb_hg38_idx', 
                                                   chrom, start, \"end\") 
                   AS frags
                   FROM read_parquet('AE_Hs_6Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_101825.parquet')
                   WHERE chrom IN ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
                                   'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                   'chr21', 'chr22')
                   GROUP BY (name, 
                            chrom, 
                            start, 
                            \"end\")
                   ORDER BY chrom, start, \"end\", name"))
            )
#    user  system elapsed 
# 267.436  39.505 118.049 

# how many cells is this?
system.time(res3 <- 
            dbGetQuery(con,
           paste0("SELECT
                   name, count(name) AS frags 
                   FROM read_parquet('AE_Hs_6Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_101825.parquet')
                   WHERE chrom IN ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
                                   'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                   'chr21', 'chr22')
                   GROUP BY name
                   ORDER BY frags DESC")
           )
)
#    user  system elapsed 
#  32.863   7.286   0.576 



res3
# data frame with 482050 rows and 2 columns
#                          name     frags
#                   <character> <numeric>
# 1      TCTTATTA_TATTTGCG_C0..    187204
# 2      GAGCACGT_GTAAGGAG_C0..    170432
# 3      CCTCGCAG_TGCATACT_D0..    151624
# 4      TCCTGAGC_GATTCCCA_D0..    130644
# 5      CCGTTTGT_GCAGCCTC_C0..    122733
# ...                       ...       ...
# 482046 CAAGAATT_AGCCGTTC_E0..         1
# 482047 ATGATGAT_GAGACGTG_F1..         1
# 482048 TAGGCATG_AGAGATCT_B0..         1
# 482049 TAAGGCGA_CCTGCGGG_D0..         1
# 482050 AAATACAG_TGATTGTT_B0..         1


# how many usable cells?   (frags >= 250) 
system.time(res4 <- 
            dbGetQuery(con,
           paste0("SELECT
                   name, count(name) AS frags 
                   FROM read_parquet('AE_Hs_6Leuk_Ms_H3K27me3_H19_Rb_DEM_PA1241_101825.parquet')
                   WHERE chrom IN ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
                                   'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                   'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                   'chr21', 'chr22')
                   GROUP BY name
                   HAVING frags >= 250
                   ORDER BY frags DESC")
           )
)
#    user  system elapsed 
#  32.607   6.701   0.410 

res4
# data frame with 47029 rows and 2 columns
#                         name     frags
#                  <character> <numeric>
# 1     TCTTATTA_TATTTGCG_C0..    187204
# 2     GAGCACGT_GTAAGGAG_C0..    170432
# 3     CCTCGCAG_TGCATACT_D0..    151624
# 4     TCCTGAGC_GATTCCCA_D0..    130644
# 5     CCGTTTGT_GCAGCCTC_C0..    122733
# ...                      ...       ...
# 47025 TTTCGTCA_TTTCTAGC_D0..       250
# 47026 AAACGGTC_TAACCAAG_C0..       250
# 47027 ATTGCCGG_CTCTCTAT_C0..       250
# 47028 AATTCGTT_CTCTCTAT_D0..       250
# 47029 CCAACACG_ACATTCGT_G1..       250


# OK so let's see whether we can blacklist fragments before doing this 


