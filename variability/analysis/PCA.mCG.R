# Look at PCA of mCG
# Peter Hickey
# 2017-07-03

### ============================================================================
### Setup
###

library(bsseq)
library(matrixStats)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

mCG_small_sorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
colnames(mCG_small_sorted) <- gsub("NA", "NAcc", colnames(mCG_small_sorted))
mCG_small_sorted$Tissue <- gsub("NA", "NAcc", mCG_small_sorted$Tissue)
mCG_small_bulk <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))
# NOTE: Not interested in caudate so removing
mCG_small_bulk <- mCG_small_bulk[, !grepl("caudate", colnames(mCG_small_bulk))]
colnames(mCG_small_bulk) <- gsub("NA", "NAcc", colnames(mCG_small_bulk))
mCG_small_bulk$Tissue <- gsub("NA", "NAcc", mCG_small_bulk$Tissue)
# NOTE: Add bulk as NeuN type
colnames(mCG_small_bulk) <- paste0(colnames(mCG_small_bulk), "_bulk")
mCG_small_bulk$NeuN <- "bulk"

# ------------------------------------------------------------------------------
# Construct SummarizedExperiment with methylation data for common CpGs
#

ol <- findOverlaps(mCG_small_sorted, mCG_small_bulk)
meth_sorted <- as.matrix(getMeth(mCG_small_sorted))[queryHits(ol), ]
meth_bulk <- as.matrix(getMeth(mCG_small_bulk))[subjectHits(ol), ]
raw_meth_sorted <-
  as.matrix(getMeth(mCG_small_sorted, type = "raw"))[queryHits(ol), ]
raw_meth_bulk <-
  as.matrix(getMeth(mCG_small_bulk, type = "raw"))[subjectHits(ol), ]

combined <- SummarizedExperiment(
  assays = SimpleList(meth = cbind(meth_sorted, meth_bulk),
                      raw_meth = cbind(raw_meth_sorted, raw_meth_bulk)),
  rowRanges = rowRanges(mCG_small_sorted)[queryHits(ol)],
  colData = combine(as.data.frame(colData(mCG_small_sorted)),
                    as.data.frame(colData(mCG_small_bulk))))
rm(meth_sorted, meth_bulk, raw_meth_sorted, raw_meth_bulk)
seqinfo(combined) <- intersect(seqinfo(combined),
                               seqinfo(BSgenome.Hsapiens.UCSC.hg19))

### ============================================================================
### Compute PCA
###

# NOTE: Useful link on relationship between PCA and SVD
#       (https://stats.stackexchange.com/q/134282)
# NOTE: Relationship between svd(x) and sva(crossprod(x))
#       svd(x)    svd(crossprod(x))
#   U   u         v (up to sign flip)
#   D   d         d^2
#   V   v         v (up to sign flip)
#
# NOTE: PCA of covariance matrix requires centering of data matrix
# NOTE: PCA of correlation matrix requires centering and scaling of data matrix

meth <- assay(combined, "meth", withDimnames = FALSE)
raw_meth <- assay(combined, "raw_meth", withDimnames = FALSE)

# ------------------------------------------------------------------------------
# Using all common CpGs
#

# All samples
meth_row_centered <- meth - rowMeans2(meth)
meth_row_centered_cp <- crossprod(meth_row_centered)
meth_row_centered_cp_svd <- svd(meth_row_centered_cp)
meth_row_centered_cp_pve <- meth_row_centered_cp_svd$d /
  sum(meth_row_centered_cp_svd$d)

# Sorted samples
meth_sorted_row_centered <- meth[, combined$NeuN != "bulk"] -
  rowMeans2(meth, cols = which(combined$NeuN != "bulk"))
meth_sorted_row_centered_cp <- crossprod(meth_sorted_row_centered)
meth_sorted_row_centered_cp_svd <- svd(meth_sorted_row_centered_cp)
rm(meth_sorted_row_centered)
meth_sorted_row_centered_cp_pve <- meth_sorted_row_centered_cp_svd$d /
  sum(meth_sorted_row_centered_cp_svd$d)

# NeuN+ samples
meth_pos_row_centered <- meth[, combined$NeuN == "pos"] -
  rowMeans2(meth, cols = which(combined$NeuN == "pos"))
meth_pos_row_centered_cp <- crossprod(meth_pos_row_centered)
meth_pos_row_centered_cp_svd <- svd(meth_pos_row_centered_cp)
rm(meth_pos_row_centered)
meth_pos_row_centered_cp_pve <- meth_pos_row_centered_cp_svd$d /
  sum(meth_pos_row_centered_cp_svd$d)

# NeuN- samples
meth_neg_row_centered <- meth[, combined$NeuN == "neg"] -
  rowMeans2(meth, cols = which(combined$NeuN == "neg"))
meth_neg_row_centered_cp <- crossprod(meth_neg_row_centered)
meth_neg_row_centered_cp_svd <- svd(meth_neg_row_centered_cp)
rm(meth_neg_row_centered)
meth_neg_row_centered_cp_pve <- meth_neg_row_centered_cp_svd$d /
  sum(meth_neg_row_centered_cp_svd$d)

# Bulk samples
meth_bulk_row_centered <- meth[, combined$NeuN == "bulk"] -
  rowMeans2(meth, cols = which(combined$NeuN == "bulk"))
meth_bulk_row_centered_cp <- crossprod(meth_bulk_row_centered)
meth_bulk_row_centered_cp_svd <- svd(meth_bulk_row_centered_cp)
rm(meth_bulk_row_centered)
meth_bulk_row_centered_cp_pve <- meth_bulk_row_centered_cp_svd$d /
  sum(meth_bulk_row_centered_cp_svd$d)

# ------------------------------------------------------------------------------
# Using average methylation in 1 kb bins
#

bins <- tileGenome(seqlengths(combined),
                   tilewidth = 1000,
                   cut.last.tile.in.chrom = TRUE)
ol <- findOverlaps(bins, combined)

binned_meth <- do.call(rbind, lapply(as.list(ol), function(i) {
  colMeans2(meth, rows = i)
}))
binned_meth <- binned_meth[!rowAlls(binned_meth, value = NaN), ]
binned_meth_row_centered <- binned_meth - rowMeans2(binned_meth)
binned_meth_row_centered_cp <- crossprod(binned_meth_row_centered)
binned_meth_row_centered_cp_svd <- svd(binned_meth_row_centered_cp)
binned_meth_row_centered_cp_pve <- binned_meth_row_centered_cp_svd$d /
  sum(binned_meth_row_centered_cp_svd$d)

binned_raw_meth <- do.call(rbind, lapply(as.list(ol), function(i) {
  colMeans2(raw_meth, rows = i)
}))
binned_raw_meth <- binned_raw_meth[!rowAlls(binned_raw_meth, value = NaN), ]
binned_raw_meth_row_centered <- binned_raw_meth - rowMeans2(binned_raw_meth)
binned_raw_meth_row_centered_cp <- crossprod(binned_raw_meth_row_centered)
binned_raw_meth_row_centered_cp_svd <- svd(binned_raw_meth_row_centered_cp)
binned_raw_meth_row_centered_cp_pve <-  binned_raw_meth_row_centered_cp_svd$d /
  sum(binned_raw_meth_row_centered_cp_svd$d)

# ------------------------------------------------------------------------------
# Save objects
#

flow_sorting_project.mCG.pca <- meth_row_centered_cp_svd
flow_sorting_project.sorted.mCG.pca <- meth_sorted_row_centered_cp_svd
flow_sorting_project.pos.mCG.pca <- meth_pos_row_centered_cp_svd
flow_sorting_project.neg.mCG.pca <- meth_neg_row_centered_cp_svd
flow_sorting_project.bulk.mCG.pca <- meth_bulk_row_centered_cp_svd
flow_sorting_project.binned_smooth_mCG.pca <- binned_meth_row_centered_cp_svd
flow_sorting_project.binned_raw_mCG.pca <- binned_raw_meth_row_centered_cp_svd
flow_sorting_project.mCG.pve <- meth_row_centered_cp_pve
flow_sorting_project.sorted.mCG.pve <- meth_sorted_row_centered_cp_pve
flow_sorting_project.pos.mCG.pve <- meth_pos_row_centered_cp_pve
flow_sorting_project.neg.mCG.pve <- meth_neg_row_centered_cp_pve
flow_sorting_project.bulk.mCG.pve <- meth_bulk_row_centered_cp_pve
flow_sorting_project.binned_smooth_mCG.pve <- binned_meth_row_centered_cp_pve
flow_sorting_project.binned_raw_mCG.pve <- binned_raw_meth_row_centered_cp_pve
flow_sorting_project_colData <- colData(combined)
save(flow_sorting_project.mCG.pca,
     flow_sorting_project.sorted.mCG.pca,
     flow_sorting_project.pos.mCG.pca,
     flow_sorting_project.neg.mCG.pca,
     flow_sorting_project.bulk.mCG.pca,
     flow_sorting_project.binned_smooth_mCG.pca,
     flow_sorting_project.binned_raw_mCG.pca,
     flow_sorting_project.mCG.pve,
     flow_sorting_project.sorted.mCG.pve,
     flow_sorting_project.pos.mCG.pve,
     flow_sorting_project.neg.mCG.pve,
     flow_sorting_project.bulk.mCG.pve,
     flow_sorting_project.binned_smooth_mCG.pve,
     flow_sorting_project.binned_raw_mCG.pve,
     flow_sorting_project_colData,
     file = "../objects/flow_sorting_project.mCG.PCA.RData")

# ------------------------------------------------------------------------------
# Plots
#

pdf("../figures/mCG.PCA.pdf")
plot(flow_sorting_project.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color,
     pch = flow_sorting_project_colData$NeuN,
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.mCG.pve[2], 1),
                   "%)"),
     main = "meth_row_centered_cp")

plot(flow_sorting_project.binned_smooth_mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color,
     pch = flow_sorting_project_colData$NeuN,
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.binned_smooth_mCG.pve[1],
                         1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.binned_smooth_mCG.pve[2],
                         1),
                   "%)"),
     main = "binned_smooth_meth_row_centered_cp")

plot(flow_sorting_project.binned_raw_mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color,
     pch = flow_sorting_project_colData$NeuN,
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.binned_raw_mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.binned_raw_mCG.pve[2], 1),
                   "%)"),
     main = "binned_raw_meth_row_centered_cp")

plot(flow_sorting_project.sorted.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color[
       flow_sorting_project_colData$NeuN != "bulk"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN != "bulk"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.sorted.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.sorted.mCG.pve[2], 1),
                   "%)"),
     main = "meth_sorted_row_centered_cp")

plot(flow_sorting_project.pos.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color[
       flow_sorting_project_colData$NeuN == "pos"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN == "pos"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.pos.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.pos.mCG.pve[2], 1),
                   "%)"),
     main = "meth_pos_row_centered_cp")

plot(flow_sorting_project.neg.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color[
       flow_sorting_project_colData$NeuN == "neg"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN == "neg"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.neg.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.neg.mCG.pve[2], 1),
                   "%)"),
     main = "meth_neg_row_centered_cp")
plot(flow_sorting_project.neg.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Individual_color[
       flow_sorting_project_colData$NeuN == "neg"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN == "neg"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.neg.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.neg.mCG.pve[2], 1),
                   "%)"),
     main = "meth_neg_row_centered_cp")

plot(flow_sorting_project.bulk.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Tissue_color[
       flow_sorting_project_colData$NeuN == "bulk"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN == "bulk"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.bulk.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.bulk.mCG.pve[2], 1),
                   "%)"),
     main = "meth_bulk_row_centered_cp")
plot(flow_sorting_project.bulk.mCG.pca$u[, 1:2],
     col = flow_sorting_project_colData$Individual_color[
       flow_sorting_project_colData$NeuN == "bulk"],
     pch = flow_sorting_project_colData$NeuN[
       flow_sorting_project_colData$NeuN == "bulk"],
     xlab = paste0("PC1 (",
                   round(100 * flow_sorting_project.bulk.mCG.pve[1], 1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * flow_sorting_project.bulk.mCG.pve[2], 1),
                   "%)"),
     main = "meth_bulk_row_centered_cp")

dev.off()
