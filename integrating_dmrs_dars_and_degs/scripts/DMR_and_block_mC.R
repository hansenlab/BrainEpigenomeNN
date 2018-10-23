# Average mC over CG-DMRs, CG-blocks, and CH-DMRs
# Peter Hickey
# 2017-12-07

library(bsseq)
library(matrixStats)
library(limma)

options("mc.cores" = 6)

### ============================================================================
### Load data
###

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

load("../objects/assays-and-features.rda")

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
                                 ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq[, BSseq$NeuN == "pos"]
})

# ------------------------------------------------------------------------------
# mCG
#

CG_small_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.small.sorted.somatic.all"))
CG_small_BSseq <- CG_small_BSseq[, CG_small_BSseq$NeuN == "pos"]

CG_large_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.large.sorted.somatic.all"))
CG_large_BSseq <- CG_large_BSseq[, CG_large_BSseq$NeuN == "pos"]

list_of_CG_BSseq <- list("mCG (S)" = CG_small_BSseq,
                         "mCG (L)" = CG_large_BSseq)

# ------------------------------------------------------------------------------
# mC
#

list_of_BSseq <- c(list_of_CG_BSseq, list_of_CH_BSseq)

# ------------------------------------------------------------------------------
# DMRs and blocks
#

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.rds")
# NOTE: CG-DMRs and CG-blocks have used FWER <= 0.05 rather than FWER < 0.05
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x[x$fwer <= 50]
})
DMRs_and_blocks <- c(list("mCG (S)" = dmrs_pos,
                          "mCG (L)" = makeGRangesFromDataFrame(sig_block_dmrs)),
                     list_of_CH_DMRs)

### ============================================================================
### Compute mean methylation over DMRs and blocks
###

# NOTE: Rather than directly using getMeth(), use a faster hack that laods all
#       data into memory
mC <- mclapply(list_of_BSseq, function(BSseq, list_of_regions, fac) {
  fac <- BSseq[[fac]]
  flattened_regions <- reduce(
    unlist(GRangesList(lapply(list_of_regions, granges)),
           use.names = TRUE))
  BSseq <- subsetByOverlaps(BSseq, flattened_regions)
  meth <- getMeth(BSseq)
  meth_all <- as.matrix(meth)
  # Some of next steps pinched from internals of bsseq::getMeth()
  grBSseq <- granges(BSseq)
  x <- lapply(names(list_of_regions), function(n) {
    regions <- list_of_regions[[n]]
    ov <- findOverlaps(grBSseq, regions)
    meth <- meth_all[queryHits(ov), , drop = FALSE]
    # # Base-level condition-specific average mC
    # base_level_ave_mC <- sapply(unique(fac), function(condition) {
    #   rowMeans2(meth, col = grep(condition, colnames(meth)))
    # })
    # Region-level average mC
    out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
    out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)), ] <- out
    rownames(outMatrix) <- paste0(n, ".", seq_len(nrow(outMatrix)))
    # Region-level condition-specific average mC
    ave_mC <- t(apply(outMatrix, 1,
                      function(xx) tapply(xx, fac, mean, na.rm = TRUE)))
    list(sample_level = outMatrix,
         condition_level = ave_mC)
  })
  names(x) <- names(list_of_regions)
  x
}, list_of_regions = DMRs_and_blocks, fac = "Tissue")

list_of_SEs <- lapply(names(DMRs_and_blocks), function(n) {
  regions <- DMRs_and_blocks[[n]]
  ave_mC_sample_level <- lapply(lapply(mC, "[[", n), "[[", "sample_level")
  ave_mC_condition_level <- lapply(lapply(mC, "[[", n), "[[", "condition_level")
  list(sample_level = SummarizedExperiment(assays = ave_mC_sample_level,
                                           rowRanges = unname(regions)),
       condition_level = SummarizedExperiment(assays = ave_mC_condition_level,
                                              rowRanges = unname(regions)))
})
names(list_of_SEs) <- names(DMRs_and_blocks)

saveRDS(list_of_SEs, "../objects/list_of_SEs.DMRs_and_blocks.rds",
        compress = "xz")
