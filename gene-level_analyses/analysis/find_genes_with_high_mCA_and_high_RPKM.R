# Find genes with high mCA that are also highly expressed
# Peter Hickey
# 2018-03-07

### ----------------------------------------------------------------------------
### Setup
###

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "mC_around_scaled_gene.rds"))

### ----------------------------------------------------------------------------
### Data wrangling
###

sx <- filter(x, bin >= 200, bin <= 300) %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  group_by(condition, context, gene) %>%
  summarise(mC = weighted.mean(mC, nC, names = TRUE)) %>%
  inner_join(genes_df)

### ----------------------------------------------------------------------------
### Find genes with high mCA and high RPKM
###

e <- new.env()
load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda", e)

z <- ungroup(sx) %>%
  filter(condition == "BA9_pos",
         context == "CpA (same strand)") %>%
  mutate(r1 = length(rpkm_NAcc_pos) - rank(rpkm_NAcc_pos, na.last = FALSE),
         r2 = length(mC) - rank(mC, na.last = FALSE),
         mr = (r1 + r2) / 2) %>%
  select(gene, rpkm_NAcc_pos, r1, mC, r2, mr) %>%
  arrange(mr) %>%
  inner_join(select(e$rna_atac_meth, gene, gene_symbol))

# TODO: It's not clear where to cut this table (WIP)
write.table(select(z, gene_symbol),
            "~/kraken.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
