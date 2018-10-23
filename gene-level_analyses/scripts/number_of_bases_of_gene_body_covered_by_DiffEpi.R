# Number of bases of gene body covered by DiffEpi
# Peter Hickey
# 2018-03-08

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))

CH_DMRs <- reduce(unlist(GRangesList(list_of_CH_DMRs)), ignore.strand = TRUE)

### ----------------------------------------------------------------------------
### Functions
###

NumberOfBasesOfGeneBodyCoveredByDiffEpi <- function(genes, dmrs) {
  op <- findOverlapPairs(genes, dmrs)
  w <- width(pintersect(op))
  tw <- tapply(w, names(S4Vectors::first(op)), sum)
  val <- setNames(rep(0L, length(genes)), names(genes))
  val[names(tw)] <- tw
  data_frame(gene = names(val),
             n_bases_covered = val)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, dmrs_NAvsBA9pos) %>%
    mutate(DiffEpi = "CG-DMRs"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, blocks_pos) %>%
    mutate(DiffEpi = "CG-blocks"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, dars_pos) %>%
    mutate(DiffEpi = "DARs"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, big_dars_pos) %>%
    mutate(DiffEpi = "bigDARs"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, list_of_CH_DMRs[["mCA (+)"]]) %>%
    mutate(DiffEpi = "CA_pos-DMRs"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, list_of_CH_DMRs[["mCA (-)"]]) %>%
    mutate(DiffEpi = "CA_neg-DMRs"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, invertStrand(list_of_CH_DMRs[["mCA (+)"]])) %>%
    mutate(DiffEpi = "CA_pos-DMRs (opposite strand)"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, invertStrand(list_of_CH_DMRs[["mCA (-)"]])) %>%
    mutate(DiffEpi = "CA_neg-DMRs (opposite strand)"),
  NumberOfBasesOfGeneBodyCoveredByDiffEpi(
    genes, CH_DMRs) %>%
    mutate(DiffEpi = "CH-DMRs"))

saveRDS(
  val,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "number_of_bases_of_gene_body_covered_by_DiffEpi.rds"))
