# Compute condition-specific average mCA and mCT in 1kb bins along genome
# Peter Hickey
# 2018-02-19

### ----------------------------------------------------------------------------
### Setup
###

library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
library(dplyr)

options("mc.cores" = 8)
options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "pos_CA_small-flow-sorted-brain-wgbs"))
CA_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "neg_CA_small-flow-sorted-brain-wgbs"))
CT_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "pos_CT_small-flow-sorted-brain-wgbs"))
CT_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "neg_CT_small-flow-sorted-brain-wgbs"))

bins <- tileGenome(seqlengths(CA_pos_BSseq), tilewidth = 1000,
                   cut.last.tile.in.chrom = TRUE)
strand(bins) <- "+"
bins <- c(bins, invertStrand(bins))

conditions_df <- data_frame(condition = c("NAcc_pos", "BA9_pos"),
                            tissue = c("NA", "BA9"),
                            NeuN = c("pos", "pos"))

### ----------------------------------------------------------------------------
### Functions
###

# Returns a SummarizedExperiment
AveragemCInBins <- function(BSseq, bins, conditions_df, min_cov = 1) {
  list_of_se <- mclapply(seq_len(nrow(conditions_df)), function(i) {
    message("Subsetting BSseq: ", conditions_df$condition[i])
    BSseq <- BSseq[, BSseq$Tissue %in% conditions_df$tissue[i] &
                     BSseq$NeuN %in% conditions_df$NeuN[i]]
    message("Computing meth: ", conditions_df$condition[i])
    meth <- as.matrix(getMeth(BSseq, type = "smooth"))
    keep <- rowMins(getCoverage(BSseq)) >= min_cov
    # TODO: Constructing a SE is not strictly necessary, but does help
    #       keep GRanges and meth parallel
    se <- SummarizedExperiment(
      assay = as.matrix(rowMeans2(meth[keep, , drop = FALSE])),
      rowRanges = rowRanges(BSseq)[keep],
      colData = DataFrame(row.names = conditions_df$condition[i]))
    # Compute average methylation in each bin
    ol <- findOverlaps(bins, se)
    ol_list <- as(ol, "List")
    m <- assay(se, withDimnames = FALSE)
    val <- vapply(ol_list, function(i) {
      mean2(m, idxs = i, na.rm = TRUE)
    }, numeric(1))
    SummarizedExperiment(
      assay = as.matrix(val),
      rowRanges = bins,
      colData = DataFrame(row.names = conditions_df$condition[i]))
  }, mc.cores = 2)
  do.call(cbind, list_of_se)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

list_of_BSseq <- list("CpA (+)" = CA_pos_BSseq,
                      "CpA (-)" = CA_neg_BSseq,
                      "CpT (+)" = CT_pos_BSseq,
                      "CpT (-)" = CT_neg_BSseq)
val <- mclapply(list_of_BSseq,
                AveragemCInBins,
                bins = bins,
                conditions_df = conditions_df,
                mc.cores = 4)
val[c("CpA (+)", "CpT (+)")] <- lapply(
  val[c("CpA (+)", "CpT (+)")], function(x) {
    x[strand(x) == "+", ]
  })
val[c("CpA (-)", "CpT (-)")] <- lapply(
  val[c("CpA (-)", "CpT (-)")], function(x) {
    x[strand(x) == "-", ]
  })
val <- do.call(cbind, lapply(names(val), function(n) {
  colnames(val[[n]]) <- paste0(n, " ", colnames(val[[n]]))
  unstrand(val[[n]])
}))
saveRDS(val, "../objects/average_mCA_and_mCT_in_1kb_bins.rds")
