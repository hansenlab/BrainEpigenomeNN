# DiffEpi measure around genes using equal number of bins along gene body +/- 2
# gene body equivalents (GBE)
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
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "gene-level_analyses_data.rda"))

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
DiffEpiMeasureAroundScaledGene <- function(genes, feature, yname, gbe = 2,
                                           ngbtile = 100) {
  message("Padding and tiling genes")
  # NOTE: Drop genes that are too short to generate `ngbtile` tiles
  genes <- genes[width(genes) >= ngbtile]
  # Pad the genes with GBE upstream and downstream of the gene
  # NOTE: This can create 'out-of-bounds' ranges; I suppress these warnings
  padded_genes <- suppressWarnings(
    punion(flank(genes, gbe * width(genes), start = TRUE),
           flank(genes, gbe * width(genes), start = FALSE),
           fill.gap = TRUE))
  ntile <- ngbtile * (1 + 2 * gbe)
  tiled_padded_genes <- tile(padded_genes, ntile)
  names(tiled_padded_genes) <- names(padded_genes)
  # NOTE: Reverse order of bins for negative strand genes so that bins go in
  #       ascending order from upstream to downstream
  tiled_padded_genes[strand(padded_genes) == "-"] <-
    endoapply(tiled_padded_genes[strand(padded_genes) == "-"], rev)
  unlisted_tiled_padded_genes <- unlist(tiled_padded_genes)
  message("Finding overlaps")
  ol <- findOverlaps(unlisted_tiled_padded_genes, feature)
  val <- rep(NA_real_, length(unlisted_tiled_padded_genes))
  message("Computing median DiffEpi_measure")
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  data_frame(bin = rep(seq_len(ntile),
                       times = length(tiled_padded_genes)),
             gene = rep(names(tiled_padded_genes), each = ntile),
             DiffEpi_measure = val)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  DiffEpiMeasureAroundScaledGene(
    genes, dmrs_NAvsBA9pos, "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CG-DMRs"),
  DiffEpiMeasureAroundScaledGene(
    genes, blocks_pos, "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CG-blocks"),
  DiffEpiMeasureAroundScaledGene(
    genes, dars_pos, "logFC", gbe = 20) %>%
    mutate(DiffEpi = "DARs"),
  DiffEpiMeasureAroundScaledGene(
    genes, big_dars_pos, "logFC", gbe = 20) %>%
    mutate(DiffEpi = "bigDARs"),
  DiffEpiMeasureAroundScaledGene(
    genes, list_of_CH_DMRs[["mCA (+)"]], "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CA_pos-DMRs"),
  DiffEpiMeasureAroundScaledGene(
    genes, list_of_CH_DMRs[["mCA (-)"]], "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CA_neg-DMRs"),
  DiffEpiMeasureAroundScaledGene(
    genes, invertStrand(list_of_CH_DMRs[["mCA (+)"]]), "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CA_pos-DMRs (opposite strand)"),
  DiffEpiMeasureAroundScaledGene(
    genes, invertStrand(list_of_CH_DMRs[["mCA (-)"]]), "meanDiff", gbe = 20) %>%
    mutate(DiffEpi = "CA_neg-DMRs (opposite strand)"))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "DiffEpi_measure_around_scaled_gene.rds"))
