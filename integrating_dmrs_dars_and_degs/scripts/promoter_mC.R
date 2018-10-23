# Average mC over gene promoters
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

### ============================================================================
### Compute mean methylation over gene promoters
###

# NOTE: Rather than directly using getMeth(), use a faster hack that laods all
#       data into memory
promoter_mC <- mclapply(list_of_BSseq, function(BSseq, regions, fac) {
  fac <- BSseq[[fac]]
  unlisted_regions <- unlist(regions)
  BSseq <- subsetByOverlaps(BSseq, regions, ignore.strand = TRUE)
  meth <- getMeth(BSseq)
  meth_all <- as.matrix(meth)
  # Some of next steps pinched from internals of bsseq::getMeth()
  grBSseq <- granges(BSseq)
  # Promoter-level condition-specific average mC
  ov <- findOverlaps(grBSseq, unlisted_regions, ignore.strand = TRUE)
  meth <- meth_all[queryHits(ov), , drop = FALSE]
  out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
  out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
  outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(unlisted_regions))
  colnames(outMatrix) <- sampleNames(BSseq)
  outMatrix[as.integer(rownames(out)), ] <- out
  promoter_level_ave_mC <- t(apply(outMatrix, 1,
                                   function(xx) tapply(xx, fac, mean,
                                                       na.rm = TRUE)))
  promoter_level_ave_mC <- relist(promoter_level_ave_mC, regions)
  # Average over all promoters in a gene
  ov <- findOverlaps(grBSseq, regions, ignore.strand = TRUE)
  meth <- meth_all[queryHits(ov), , drop = FALSE]
  out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
  out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
  outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
  colnames(outMatrix) <- sampleNames(BSseq)
  outMatrix[as.integer(rownames(out)), ] <- out
  ave_promoter_level_ave_mC <- t(apply(outMatrix, 1,
                                       function(xx) tapply(xx, fac, mean,
                                                           na.rm = TRUE)))
  rownames(ave_promoter_level_ave_mC) <- names(regions)

  list(promoter_level = promoter_level_ave_mC,
       ave_promoter_level = ave_promoter_level_ave_mC)
}, regions = promoters_by_gene, fac = "Tissue")

promoter_mC_base_level <- lapply(promoter_mC, "[[", "promoter_level")
saveRDS(promoter_mC_base_level,
        file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "promoter_mC_base_level.rds"))
promoter_mC_ave_promoter_level <-
  lapply(promoter_mC, "[[", "ave_promoter_level")
saveRDS(promoter_mC_ave_promoter_level,
        "../objects/promoter_mC_ave_promoter_level.rds")
