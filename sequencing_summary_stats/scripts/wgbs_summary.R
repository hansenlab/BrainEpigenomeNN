# WGBS summary statistics
# Peter Hickey
# 2017-01-26

library(dplyr)
library(purrr)
library(bsseq)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)

#-------------------------------------------------------------------------------
# Parse files
#

# Parse Bismark output and import into R
sequenced <- system("grep 'Sequence pairs analysed in total' /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq*/*/*PE_report.txt", intern = TRUE)
aligned <- system("grep 'Number of paired-end alignments with a unique best hit' /dcl01/feinberg/data/gtex/flow_sorted_nuclei/HiSeq*/*/*PE_report.txt", intern = TRUE)

# More parsing of Bismark output (now within R)
sequenced_df <- map_df(strsplit(sequenced, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 8L),
             `Number sequenced PE reads` = as.numeric(x[[2]]))
})
aligned_df <- map_df(strsplit(aligned, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 8L),
             `Number aligned PE reads` = as.numeric(x[[2]]))
})

# Clean sample IDs
clean_ID <- function(x) {
  x %>%
    mutate(ID = sub("Sample_", "", ID),
           ID = sub("-", "_", ID),
           ID = sub("NA", "NAcc", ID),
           ID = sub("PFC", "BA9", ID),
           ID = ifelse(grepl("_", ID), ID,
                       paste0(substr(ID, 1, 4), "_",
                              substr(ID, 5, nchar(ID)))),
           ID = sub("BA9_5086", "5086_BA9", ID))
}

sequenced_df <- sequenced_df %>%
  clean_ID()
aligned_df <- aligned_df %>%
  clean_ID()

#-------------------------------------------------------------------------------
# Aggregate at sample level
#

sequenced_df <- sequenced_df %>%
  group_by(ID) %>%
  summarise(`Number sequenced PE reads` = sum(`Number sequenced PE reads`))
aligned_df <- aligned_df %>%
  group_by(ID) %>%
  summarise(`Number aligned PE reads` = sum(`Number aligned PE reads`))

#-------------------------------------------------------------------------------
# Alignment rate
#

df <- inner_join(sequenced_df, aligned_df, by = c("ID" = "ID"))
df <- df %>%
  mutate(`Alignment rate (%)` =
           round(100 * `Number aligned PE reads` / `Number sequenced PE reads`,
                 0))

#-------------------------------------------------------------------------------
# Bisulfite-conversion rate
#

# Load objects
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/notsorted/BS.unsorted.fit.small.rda")
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.rda")

# Clean sample IDs
colnames(BS.unsorted.fit.small) <-
  data_frame(ID = colnames(BS.unsorted.fit.small)) %>%
  clean_ID() %>%
  .$ID
colnames(BS.fit.small.sorted) <-
  data_frame(ID = colnames(BS.fit.small.sorted)) %>%
  clean_ID() %>%
  .$ID

# Bisulfite-conversion rates
lambda <- GRanges("lambda",
                  IRanges(1, 48502))
bs_conversion_rate_df <-
  data_frame(ID = c(colnames(BS.unsorted.fit.small),
                    colnames(BS.fit.small.sorted)),
             `Bisulfite conversion rate (%)` =
               round(100 *
                       c(1 - colSums(getCoverage(BSseq = BS.unsorted.fit.small,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS.unsorted.fit.small,
                                               regions = lambda,
                                               type = "Cov")[[1L]]),
                         1 - colSums(getCoverage(BSseq = BS.fit.small.sorted,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS.fit.small.sorted,
                                               regions = lambda,
                                               type = "Cov")[[1L]])),
                     2))
df <- inner_join(df, bs_conversion_rate_df)

#-------------------------------------------------------------------------------
# Coverage-related metrics
#

# Keep autosomal data only
BS.unsorted.fit.small <-
  keepSeqlevels(BS.unsorted.fit.small, paste0("chr", 1:22))
BS.fit.small.sorted <-
  keepSeqlevels(BS.fit.small.sorted, paste0("chr", 1:22))

# Number of CpGs on autosomes in hg19
n_CpGs <- Reduce(sum, bsapply(BSParams = new("BSParams",
                                             X = BSgenome.Hsapiens.UCSC.hg19,
                                             FUN = countPattern,
                                             exclude = c("M", "_", "X", "Y")),
                              pattern = "CG"))

# Compute coverage-related metrics
coverage_df <-
  data_frame(ID = c(colnames(BS.unsorted.fit.small),
                    colnames(BS.fit.small.sorted)),
             `Number covered CpGs` =
               c(colSums(getCoverage(BS.unsorted.fit.small) > 0),
                 colSums(getCoverage(BS.fit.small.sorted) > 0)),
             `Covered CpGs (%)` =
               round(100 * `Number covered CpGs` / n_CpGs, 0),
             `Mean depth` = c(colMeans(getCoverage(BS.unsorted.fit.small)),
                              colMeans(getCoverage(BS.fit.small.sorted))))

df <- inner_join(df, coverage_df, c("ID" = "ID"))

#-------------------------------------------------------------------------------
# Arrange columns and split into sorted and unsorted data
#

df <- select(df, ID, `Number sequenced PE reads`, `Number aligned PE reads`,
             `Alignment rate (%)`, `Number covered CpGs`, `Covered CpGs (%)`,
             `Mean depth`, `Bisulfite conversion rate (%)`)

df_list <- split(df, grepl("pos|neg", df$ID))
write_csv(x = df_list[[1]],
          path = "../tables/Summary_of_Unsorted_WGBS.csv",
          col_names = TRUE)
write_csv(x = df_list[[2]],
          path = "../tables/Summary_of_Sorted_WGBS.csv",
          col_names = TRUE)

