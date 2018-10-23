# Dinucleotide proportions around genes using equal number of bins along gene
# body +/- 2 gene body equivalents (GBE)
# Peter Hickey
# 2018-02-18

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)
library(Biostrings)

options("mc.cores" = 16)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
dinucs <- levels(interaction(DNA_BASES, DNA_BASES, sep = ""))
names(dinucs) <- dinucs

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
DinucleotideCountAroundScaledGene <- function(genes, dinuc,
                                              opposite_strand = FALSE, gbe = 2,
                                              ngbtile = 100) {
  message("Loading GRanges: ", dinuc)
  feature <- readRDS(file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                               paste0(dinuc, ".GRanges.rds")))
  if (opposite_strand) {
    feature <- invertStrand(feature)
  }
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
  message("Counting overlaps")
  val <- countOverlaps(unlisted_tiled_padded_genes, feature)
  message("Constructing result")
  data_frame(dinuc = dinuc,
             bin = rep(seq_len(ntile),
                       times = length(tiled_padded_genes)),
             gene = rep(names(tiled_padded_genes), each = ntile),
             n_dinucleotide = val,
             width = width(unlisted_tiled_padded_genes))
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(mclapply(dinucs, function(dinuc) {
  DinucleotideCountAroundScaledGene(genes, dinuc)
}))

saveRDS(
  val,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "dinucleotide_count_around_around_scaled_gene.rds"))

val <- bind_rows(mclapply(dinucs, function(dinuc) {
  DinucleotideCountAroundScaledGene(genes, dinuc, opposite_strand = TRUE)
}))
saveRDS(
  val,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "dinucleotide_count_around_around_scaled_gene.opposite_strand.rds"))
