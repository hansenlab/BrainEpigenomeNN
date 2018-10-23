# mC levels at bins along genome
# Peter Hickey
# 2018-03-22

### ----------------------------------------------------------------------------
### Setup
###

library(bsseq)
library(HDF5Array)
library(dplyr)
library(DelayedMatrixStats)
library(BSgenome.Hsapiens.UCSC.hg19)

options("mc.cores" = 2)
options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

extdir <- "../extdata"

si <- keepSeqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg19),
                    paste0("chr", 1:22))

### ----------------------------------------------------------------------------
### Load data
###

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
# Each row is the (median) DiffEpi measure for bin
mCInBins <- function(BSseq, conditions_df, seqlengths, binsize,
                         min_cov = 1) {
  stopifnot(is(BSseq, "BSseq"))

  message("Tiling genome")
  bins <- tileGenome(seqlengths,
                     tilewidth = binsize,
                     cut.last.tile.in.chrom = TRUE)
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
    ol <- findOverlaps(bins, se)
    ol_list <- as(ol, "List")
    m <- assay(se, withDimnames = FALSE)
    val <- vapply(ol_list, function(i) {
      mean2(m, idxs = i)
    }, numeric(1))
    data_frame(condition = conditions_df$condition[i],
               bin = seq_along(bins),
               nC = countOverlaps(bins, se),
               mC = val) %>%
      bind_cols(as.data.frame(bins))
  })
  bind_rows(list_of_df)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  mCInBins(CG_BSseq, conditions_df, seqlengths(si), 10 ^ 4) %>%
    mutate(context = "mCG"),
  mCInBins(CA_pos_BSseq, conditions_df, seqlengths(si), 10 ^ 4) %>%
    mutate(context = "mCA (+)"),
  mCInBins(CA_neg_BSseq, conditions_df, seqlengths(si), 10 ^ 4) %>%
    mutate(context = "mCA (-)"))

saveRDS(val, "../objects/mC_in_bins.rds")

val <- bind_rows(
  mCInBins(CG_BSseq, conditions_df, seqlengths(si), 10 ^ 3) %>%
    mutate(context = "mCG"),
  mCInBins(CA_pos_BSseq, conditions_df, seqlengths(si), 10 ^ 3) %>%
    mutate(context = "mCA (+)"),
  mCInBins(CA_neg_BSseq, conditions_df, seqlengths(si), 10 ^ 3) %>%
    mutate(context = "mCA (-)"))

saveRDS(val, "../objects/mC_in_1kb_bins.rds")
