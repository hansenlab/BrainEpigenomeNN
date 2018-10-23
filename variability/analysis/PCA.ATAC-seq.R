# Look at PCA of ATAC-seq between samples
# Peter Hickey
# 2017-08-10

### ============================================================================
### Setup
###

library(SummarizedExperiment)
library(matrixStats)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(edgeR)

extdir <- "../extdata"

ATAC_SE <- readRDS(
  file.path("..", "..", "ATAC-seq", "objects",
            "flow-sorted-brain-atac.union_narrowPeak_reduced.se.rds"))
ATAC_SE <- ATAC_SE[, ATAC_SE$REPLICATE == "rep1"]
ATAC_SE <- keepSeqlevels(ATAC_SE, paste0("chr", 1:22), pruning.mode = "coarse")
ATAC_SE <- ATAC_SE[rowSums(cpm(assay(ATAC_SE)) > 1) >= 5, ]
counts <- assay(ATAC_SE)
cpm <- cpm(counts,
           log = TRUE)
colnames(ATAC_SE) <- gsub("NA", "NAcc", colnames(ATAC_SE))
ATAC_SE$TISSUE <- ifelse(ATAC_SE$TISSUE == "NA", "NAcc", ATAC_SE$TISSUE)

### ============================================================================
### Compute PCA
###

# ------------------------------------------------------------------------------
# Using all bins tested for differential accessibility
#

# All samples
ATAC_row_centered <- cpm - rowMeans2(cpm)
ATAC_row_centered_cp <- crossprod(ATAC_row_centered)
ATAC_row_centered_cp_svd <- svd(ATAC_row_centered_cp)
rm(ATAC_row_centered)
ATAC_row_centered_cp_pve <- ATAC_row_centered_cp_svd$d /
  sum(ATAC_row_centered_cp_svd$d)

# NeuN+ samples
ATAC_pos_row_centered <- cpm[, ATAC_SE$NEUN == "pos"] -
  rowMeans2(cpm, cols = which(ATAC_SE$NEUN == "pos"))
ATAC_pos_row_centered_cp <- crossprod(ATAC_pos_row_centered)
ATAC_pos_row_centered_cp_svd <- svd(ATAC_pos_row_centered_cp)
rm(ATAC_pos_row_centered)
ATAC_pos_row_centered_cp_pve <- ATAC_pos_row_centered_cp_svd$d /
  sum(ATAC_pos_row_centered_cp_svd$d)

# NeuN- samples
ATAC_neg_row_centered <- cpm[, ATAC_SE$NEUN == "neg"] -
  rowMeans2(cpm, cols = which(ATAC_SE$NEUN == "neg"))
ATAC_neg_row_centered_cp <- crossprod(ATAC_neg_row_centered)
ATAC_neg_row_centered_cp_svd <- svd(ATAC_neg_row_centered_cp)
rm(ATAC_neg_row_centered)
ATAC_neg_row_centered_cp_pve <- ATAC_neg_row_centered_cp_svd$d /
  sum(ATAC_neg_row_centered_cp_svd$d)

# ------------------------------------------------------------------------------
# Save objects
#

flow_sorting_project.ATAC.logcpm.pca <- ATAC_row_centered_cp_svd
flow_sorting_project.pos.ATAC.logcpm.pca <- ATAC_pos_row_centered_cp_svd
flow_sorting_project.neg.ATAC.logcpm.pca <- ATAC_neg_row_centered_cp_svd
flow_sorting_project_colData <- colData(ATAC_SE)
flow_sorting_project.ATAC.logcpm.pve <- ATAC_row_centered_cp_pve
flow_sorting_project.pos.ATAC.logcpm.pve <- ATAC_pos_row_centered_cp_pve
flow_sorting_project.neg.ATAC.logcpm.pve <- ATAC_neg_row_centered_cp_pve

save(flow_sorting_project.ATAC.logcpm.pca,
     flow_sorting_project.pos.ATAC.logcpm.pca,
     flow_sorting_project.neg.ATAC.logcpm.pca,
     flow_sorting_project_colData,
     flow_sorting_project.ATAC.logcpm.pve,
     flow_sorting_project.pos.ATAC.logcpm.pve,
     flow_sorting_project.neg.ATAC.logcpm.pve,
     file = "../objects/flow_sorting_project.ATAC.logcpm.PCA.RData")

# ------------------------------------------------------------------------------
# Plots
#

pdf("../figures/ATAC-seq.PCA.pdf")
plot(flow_sorting_project.ATAC.logcpm.pca$u[, 1:2],
     col = flow_sorting_project_colData$TISSUE_COLOR,
     pch = flow_sorting_project_colData$NEUN,
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.ATAC.logcpm.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.ATAC.logcpm.pve[2], 1),
                   "%)"),
     main = "ATAC_row_centered_cp")
plot(flow_sorting_project.ATAC.logcpm.pca$u[, 1:2],
     col = flow_sorting_project_colData$DONOR_COLOR,
     pch = flow_sorting_project_colData$NEUN,
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.ATAC.logcpm.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.ATAC.logcpm.pve[2], 1),
                   "%)"),
     main = "ATAC_row_centered_cp")

plot(flow_sorting_project.pos.ATAC.logcpm.pca$u[, 1:2],
     col = flow_sorting_project_colData$TISSUE_COLOR[
       flow_sorting_project_colData$NEUN == "pos"],
     pch = flow_sorting_project_colData$NEUN[
       flow_sorting_project_colData$NEUN == "pos"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.pos.ATAC.logcpm.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.pos.ATAC.logcpm.pve[2], 1),
                   "%)"),
     main = "ATAC_pos_row_centered")

plot(flow_sorting_project.neg.ATAC.logcpm.pca$u[, 1:2],
     col = flow_sorting_project_colData$TISSUE_COLOR[
       flow_sorting_project_colData$NEUN == "neg"],
     pch = flow_sorting_project_colData$NEUN[
       flow_sorting_project_colData$NEUN == "neg"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.neg.ATAC.logcpm.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.neg.ATAC.logcpm.pve[2], 1),
                   "%)"),
     main = "ATAC_neg_row_centered")

plot(flow_sorting_project.neg.ATAC.logcpm.pca$u[, 1:2],
     col = flow_sorting_project_colData$DONOR_COLOR[
       flow_sorting_project_colData$NEUN == "neg"],
     pch = flow_sorting_project_colData$NEUN[
       flow_sorting_project_colData$NEUN == "neg"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.neg.ATAC.logcpm.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.neg.ATAC.logcpm.pve[2], 1),
                   "%)"),
     main = "ATAC_neg_row_centered")
dev.off()
