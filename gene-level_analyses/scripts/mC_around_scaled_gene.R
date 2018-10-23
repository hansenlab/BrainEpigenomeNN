# mC levels around genes using equal number of bins along gene body +/- 2 gene
# body equivalents (GBE)
# Peter Hickey
# 2018-02-16

### ----------------------------------------------------------------------------
### Setup
###

library(bsseq)
library(HDF5Array)
library(dplyr)
library(tidyr)
library(DelayedMatrixStats)

options("mc.cores" = 2)
options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))

CG_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "pos_CA_small-flow-sorted-brain-wgbs"))
CA_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "neg_CA_small-flow-sorted-brain-wgbs"))

conditions_df <- data_frame(condition = c("NAcc_pos", "BA9_pos"),
                            tissue = c("NA", "BA9"),
                            NeuN = c("pos", "pos"))

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
mCAroundScaledGene <- function(genes, BSseq, conditions_df, gbe = 2,
                               ngbtile = 100, min_cov = 1) {
  stopifnot(is(BSseq, "BSseq"))

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

  list_of_df <- mclapply(seq_len(nrow(conditions_df)), function(i) {
    message("Subsetting BSseq: ", conditions_df$condition[i])
    BSseq <- BSseq[, BSseq$Tissue %in% conditions_df$tissue[i] &
                     BSseq$NeuN %in% conditions_df$NeuN[i]]
    message("Computing meth: ", conditions_df$condition[i])
    meth <- as.matrix(getMeth(BSseq, type = "smooth"))
    keep <- rowMins(getCoverage(BSseq)) >= min_cov
    # NOTE: Constructing a SE is not strictly necessary, but does help
    #       keep GRanges and meth parallel
    se <- SummarizedExperiment(
      assay = as.matrix(rowMeans2(meth[keep, , drop = FALSE])),
      rowRanges = rowRanges(BSseq)[keep],
      colData = DataFrame(row.names = conditions_df$condition[i]))
    # Compute average methylation in each bin
    ol <- findOverlaps(unlisted_tiled_padded_genes, se)
    ol_list <- as(ol, "List")
    m <- assay(se, withDimnames = FALSE)
    val <- vapply(ol_list, function(i) {
      mean2(m, idxs = i)
    }, numeric(1))
    data_frame(condition = conditions_df$condition[i],
               bin = rep(seq_len(ntile),
                         times = length(tiled_padded_genes)),
               gene = rep(names(tiled_padded_genes), each = ntile),
               nC = countOverlaps(unlisted_tiled_padded_genes, se),
               mC = val)
  })
  bind_rows(list_of_df)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  mCAroundScaledGene(
    genes, CG_BSseq, conditions_df) %>%
    mutate(context = "CpG"),
  mCAroundScaledGene(
    genes, CA_pos_BSseq, conditions_df) %>%
    mutate(context = "CpA (+)"),
  mCAroundScaledGene(
    genes, CA_neg_BSseq, conditions_df) %>%
    mutate(context = "CpA (-)"),
  mCAroundScaledGene(
    genes, invertStrand(CA_pos_BSseq), conditions_df) %>%
    mutate(context = "CpA (+) (opposite strand)"),
  mCAroundScaledGene(
    genes, invertStrand(CA_neg_BSseq), conditions_df) %>%
    mutate(context = "CpA (-) (opposite strand)"))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "mC_around_scaled_gene.rds"))
