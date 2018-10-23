# Explore variability of ATAC data
# Peter Hickey
# 2017-11-28

library(SummarizedExperiment)
library(matrixStats)
library(edgeR)

extdir <- "../extdata"

# ==============================================================================
# Load data
#

ATAC_SE <- readRDS(
  file.path(extdir, "flow-sorted-brain-atac", "objects",
            "flow-sorted-brain-atac.union_narrowPeak_reduced.se.rds"))
# NOTE: Tidy up sample colData
colnames(ATAC_SE) <- gsub("NA", "NAcc", colnames(ATAC_SE))
ATAC_SE$TISSUE <- gsub("NA", "NAcc", ATAC_SE$TISSUE)
ATAC_SE$NEUN <- factor(ATAC_SE$NEUN,
                       levels = c("bulk", "neg", "pos"))
ATAC_SE_pos <- ATAC_SE[, ATAC_SE$NEUN == "pos"]
ATAC_SE_neg <- ATAC_SE[, ATAC_SE$NEUN == "neg"]

# Load ATAC data into memory for easier/faster manipulation and convert to cpm
atac_pos <- cpm(assay(ATAC_SE_pos), log = TRUE)
atac_neg <- cpm(assay(ATAC_SE_neg), log = TRUE)

# ==============================================================================
# Compute mean and sd of ATAC in each Tissue:NeuN combination
#

tissues <- c("BA9", "NAcc")
names(tissues) <- tissues

mean_atac_pos <- do.call(cbind, lapply(tissues, function(tissue) {
  rowMeans2(atac_pos, cols = grep(tissue, colnames(atac_pos)))
}))
mean_atac_neg <- do.call(cbind, lapply(tissues, function(tissue) {
  rowMeans2(atac_neg, cols = grep(tissue, colnames(atac_neg)))
}))

sd_atac_pos <- do.call(cbind, lapply(tissues, function(tissue) {
  sqrt(rowSds(atac_pos, cols = grep(tissue, colnames(atac_pos))))
}))
sd_atac_neg <- do.call(cbind, lapply(tissues, function(tissue) {
  sqrt(rowSds(atac_neg, cols = grep(tissue, colnames(atac_neg))))
}))

# ==============================================================================
# Fit loess to mean-sd relationships
#

fits_pos <- lapply(tissues, function(tissue) {
  lowess(mean_atac_pos[, tissue], sd_atac_pos[, tissue], f = 1 / 10)
})
fits_neg <- lapply(tissues, function(tissue) {
  lowess(mean_atac_neg[, tissue], sd_atac_neg[, tissue], f = 1 / 10)
})
fits <- c(fits_pos, fits_neg)
approx_fits <- lapply(fits, approxfun)

tissue_col <- c("BA9" = "deepskyblue",
                "NAcc" = "chocolate1")
neun_col <- c("neg" = "purple",
              "pos" = "darkgreen")

Cairo::CairoPNG("../figures/atac.smooth_scatter.png",
                width = 960,
                height = 960)
par(mfcol = c(2, 3))
spwt(x = mean_atac_pos[, 1],
     y = sd_atac_pos[, 1],
     xlim = c(-6, 8),
     ylim = c(0, 3))
plot(function(x) approx_fits[[1]](x),
     from = -6,
     to = 8,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2)
spwt(x = mean_atac_pos[, 2],
     y = sd_atac_pos[, 2],
     xlim = c(-6, 8),
     ylim = c(0, 3))
plot(function(x) approx_fits[[2]](x),
     from = -6,
     to = 8,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2)
spwt(x = mean_atac_neg[, 1],
     y = sd_atac_neg[, 1],
     xlim = c(-6, 8),
     ylim = c(0, 3))
plot(function(x) approx_fits[[3]](x),
     from = -6,
     to = 8,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 2)
spwt(x = mean_atac_neg[, 2],
     y = sd_atac_neg[, 2],
     xlim = c(-6, 8),
     ylim = c(0, 3))
plot(function(x) approx_fits[[4]](x),
     from = -6,
     to = 8,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 2)
# Or all on the same plot
plot(function(x) approx_fits[[1]](x), from = -6, to = 8,
     col = tissue_col[names(approx_fits)[1]],
     lwd = 2,
     ylim = c(0, 3),
     xlim = c(-6, 8))
plot(function(x) approx_fits[[2]](x), from = -6, to = 8,
     col = tissue_col[names(approx_fits)[2]],
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[3]](x), from = -6, to = 8,
     col = tissue_col[names(approx_fits)[3]],
     lwd = 2,
     add = TRUE,
     lty = 2)
plot(function(x) approx_fits[[4]](x), from = -6, to = 8,
     col = tissue_col[names(approx_fits)[4]],
     lwd = 2,
     add = TRUE,
     lty = 2)
dev.off()

# Loess fits on same plot faceted by brain region
pdf("../figures/atac.mean_var_loess.NeuN_facet.pdf",
    width = 7,
    height = 4,
    useDingbats = FALSE)
par(mfrow = c(1, 2))
plot(function(x) approx_fits[[1]](x), from = -6, to = 8,
     col = neun_col["pos"],
     lwd = 2,
     ylim = c(0, 3),
     xlim = c(-6, 8),
     main = "BA9")
plot(function(x) approx_fits[[3]](x), from = -6, to = 8,
     col = neun_col["neg"],
     lty = 2,
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[2]](x), from = -6, to = 8,
     col = neun_col["pos"],
     lwd = 2,
     ylim = c(0, 3),
     xlim = c(-6, 8),
     main = "NAcc")
plot(function(x) approx_fits[[4]](x), from = -6, to = 8,
     col = neun_col["neg"],
     lty = 2,
     lwd = 2,
     add = TRUE)
dev.off()

