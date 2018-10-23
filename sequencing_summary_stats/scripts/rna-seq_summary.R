# RNA-seq summary statistics
# Peter Hickey
# 2017-01-30

library(dplyr)
library(purrr)
library(readr)
library(jsonlite)

#-------------------------------------------------------------------------------
# Parse files
#

salmon <- list.dirs("/dcl01/hansen/data/flow-sorted-brain-rna-seq/data/salmon",
                    full.names = TRUE,
                    recursive = FALSE)
salmon_df <- map_df(salmon, function(x) {
  json <- read_json(paste0(x, "/aux_info/meta_info.json"))
  num_processed <- json$num_processed
  num_mapped <- json$num_mapped
  percent_mapped <- json$percent_mapped

    data_frame(ID = sub(".rep1.transcripts_quant", "", basename(x)),
             `Number sequenced PE reads` = num_processed,
             `Number quasi-mapped PE reads` = num_mapped,
             `Quasi-mapping rate (%)` = round(percent_mapped, 2))
})

#-------------------------------------------------------------------------------
# Make final data frame
#

df <- salmon_df %>%
  mutate(ID = gsub("-", "_", ID),
         ID = sub("NA", "NAcc", ID))

write_csv(df, "../tables/Summary_of_RNA-seq.csv")

