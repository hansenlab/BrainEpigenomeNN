# Average mC over gene bodies
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
### Compute mean methylation over gene bodies
###

# NOTE: Rather than directly using getMeth(), use a faster hack that laods all
#       data into memory
gene_body_mC <- mclapply(list_of_BSseq, function(BSseq, regions, fac) {
  fac <- BSseq[[fac]]
  BSseq <- subsetByOverlaps(BSseq, regions, ignore.strand = TRUE)
  meth <- getMeth(BSseq)
  meth_all <- as.matrix(meth)
  # Some of next steps pinched from internals of bsseq::getMeth()
  grBSseq <- granges(BSseq)
  ov <- findOverlaps(grBSseq, regions, ignore.strand = TRUE)
  meth <- meth_all[queryHits(ov), , drop = FALSE]
  # Base-level condition-specific average mC
  base_level_ave_mC <- sapply(unique(fac), function(condition) {
    rowMeans2(meth, col = grep(condition, colnames(meth)), na.rm = TRUE)
  })
  # Gene-level condition-specific average mC
  out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
  out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
  outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
  colnames(outMatrix) <- sampleNames(BSseq)
  outMatrix[as.integer(rownames(out)), ] <- out
  gene_level_ave_mC <- t(apply(outMatrix, 1,
                               function(xx) tapply(xx, fac, mean, na.rm = TRUE)))
  rownames(gene_level_ave_mC) <- names(regions)
  list(base_level = base_level_ave_mC,
       gene_level = gene_level_ave_mC)
}, regions = gencode_features$genes, fac = "Tissue")

gene_body_mC_base_level <- lapply(gene_body_mC, "[[", "base_level")
saveRDS(gene_body_mC_base_level,
        file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "gene_body_mC_base_level.rds"))
gene_body_mC_gene_level <- lapply(gene_body_mC, "[[", "gene_level")
saveRDS(gene_body_mC_gene_level,
        "../objects/gene_body_mC_gene_level.rds")

