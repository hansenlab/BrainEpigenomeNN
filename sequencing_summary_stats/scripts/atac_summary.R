# ATAC-seq summary statistics
# Peter Hickey
# 2017-01-30

library(dplyr)
library(purrr)
library(readr)
library(Rsamtools)

#-------------------------------------------------------------------------------
# Parse files
#

# Parse Bowtie2 output and import into R

bowtie2_logs <-
  list.files("/dcl01/hansen/data/flow-sorted-brain-atac/logs/bowtie2",
             pattern = glob2rx("*.log"),
             full.names = TRUE)
bowtie2_logs <- grep("rep2|25K|prefix|test", bowtie2_logs, value = TRUE,
                     invert = TRUE)
bowtie2_logs_df <- map_df(bowtie2_logs, function(x) {
  total_reads <- system(paste0("grep 'reads; of these:' ", x,
                               "| awk '{ print $1 }'"), intern = TRUE) %>%
    as.numeric()
  alignment_rate <- system(paste0("grep 'overall alignment rate' ", x,
                                  "| awk '{ print $1 }' "), intern = TRUE) %>%
    sub("%", "", .) %>%
    as.numeric()
  data_frame(ID = sub("\\.rep1\\.bowtie2\\.log", "", basename(x)),
             `Number sequenced PE reads` = total_reads,
             `Alignment rate (%)` = alignment_rate)
})

# Parse SAMtools flagstat and import into R

flagstat <- list.files("/dcl01/hansen/data/flow-sorted-brain-atac/data/bam",
                       pattern = glob2rx("*.flagstat"),
                       full.names = TRUE)
flagstat <- grep("rep2|25K|prefix|test", flagstat, value = TRUE,
                 invert = TRUE)
flagstat_df <- map_df(flagstat, function(x) {
  properly_paired <- system(paste0("grep 'properly paired' ", x,
                                   "| awk '{ print $6 }'"), intern = TRUE) %>%
    sub("\\(", "", .) %>%
    sub("%:-nan%\\)", "", .)
  total <- system(paste0("grep 'in total' ", x,
                         "| awk '{ print $1 }'"), intern = TRUE) %>%
    as.numeric()
  duplicates <- system(paste0("grep 'duplicates' ", x,
                              "| awk '{ print $1 }'"), intern = TRUE) %>%
    as.numeric()
  duplicate_rate <- round(100 * duplicates / total, 2)

  data_frame(ID = sub("\\.rep1\\.markdup\\.bam\\.flagstat", "", basename(x)),
             total = total,
             `Properly paired (%)` = properly_paired,
             `Duplicate rate (%)` = duplicate_rate)
})

# Mitochondria rate

bams <- list.files("/dcl01/hansen/data/flow-sorted-brain-atac/data/bam",
                   pattern = glob2rx("*.markdup.bam"),
                   full.names = TRUE)
bams <- grep("rep2|25K|prefix|test", bams, value = TRUE,
             invert = TRUE)
chrM_df <- map_df(bams, function(x) {
  chrM <- GRanges("chrM", IRanges(1, 16571))
  cb <- countBam(file = x,
                 index = sub("\\.bam", "\\.bai", x),
                 param = ScanBamParam(which = chrM))
  data_frame(ID = sub("\\.rep1\\.markdup\\.bam", "", basename(x)),
             chrM = cb$records)
})

#-------------------------------------------------------------------------------
# Make final data frame
#

df <- inner_join(bowtie2_logs_df, flagstat_df) %>%
  inner_join(chrM_df) %>%
  mutate(ID = gsub("-", "_", ID),
         ID = sub("NA", "NAcc", ID),
    `Mitochondrial contamination (%)` = round(100 * chrM / total, 2)) %>%
  select(-chrM, -total)
write_csv(df, "../tables/Summary_of_ATAC-seq.csv")

