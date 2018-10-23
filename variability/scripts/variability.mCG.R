# Explore variability of mCG data
# Peter Hickey
# 2017-11-27

library(bsseq)
library(matrixStats)
library(ggplot2)

extdir <- "../extdata"

# ==============================================================================
# Load data
#

mCG_small_sorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
# NOTE: Tidy up sample colData
colnames(mCG_small_sorted) <- gsub("NA", "NAcc", colnames(mCG_small_sorted))
mCG_small_sorted$Tissue <- gsub("NA", "NAcc", mCG_small_sorted$Tissue)
mCG_small_sorted$NeuN <- factor(mCG_small_sorted$NeuN,
                                levels = c("bulk", "neg", "pos"))
mCG_small_pos <- mCG_small_sorted[, mCG_small_sorted$NeuN == "pos"]
mCG_small_neg <- mCG_small_sorted[, mCG_small_sorted$NeuN == "neg"]

mCG_small_bulk <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))
# NOTE: Tidy up sample colData
colnames(mCG_small_bulk) <- gsub("NA", "NAcc", colnames(mCG_small_bulk))
mCG_small_bulk$Tissue <- gsub("NA", "NAcc", mCG_small_bulk$Tissue)
mCG_small_bulk$NeuN <- factor(rep("bulk", ncol(mCG_small_bulk)),
                              levels = c("bulk", "neg", "pos"))
# NOTE: Not interested in caudate so removing
mCG_small_bulk <- mCG_small_bulk[, !grepl("caudate", colnames(mCG_small_bulk))]

# Load methylation data into memory for easier/faster manipulation
mCG_bulk <- as.matrix(getMeth(mCG_small_bulk))
mCG_pos <- as.matrix(getMeth(mCG_small_pos))
mCG_neg <- as.matrix(getMeth(mCG_small_neg))

# ==============================================================================
# Compute mean and sd of methylation in each Tissue:NeuN combination
#

tissues <- c("BA9", "BA24", "HC", "NAcc")
names(tissues) <- tissues

mean_mCG_bulk <- do.call(cbind, lapply(tissues, function(tissue) {
  rowMeans2(mCG_bulk, cols = grep(tissue, colnames(mCG_bulk)))
}))
mean_mCG_pos <- do.call(cbind, lapply(tissues, function(tissue) {
  rowMeans2(mCG_pos, cols = grep(tissue, colnames(mCG_pos)))
}))
mean_mCG_neg <- do.call(cbind, lapply(tissues, function(tissue) {
  rowMeans2(mCG_neg, cols = grep(tissue, colnames(mCG_neg)))
}))

sd_mCG_bulk <- do.call(cbind, lapply(tissues, function(tissue) {
  rowSds(mCG_bulk, cols = grep(tissue, colnames(mCG_bulk)))
}))
sd_mCG_pos <- do.call(cbind, lapply(tissues, function(tissue) {
  rowSds(mCG_pos, cols = grep(tissue, colnames(mCG_pos)))
}))
sd_mCG_neg <- do.call(cbind, lapply(tissues, function(tissue) {
  rowSds(mCG_neg, cols = grep(tissue, colnames(mCG_neg)))
}))

# ==============================================================================
# Fit loess to mean-sd relationships
#

fits_bulk <- lapply(tissues, function(tissue) {
  i <- sample(nrow(mean_mCG_bulk), 10 ^ 6)
  lowess(mean_mCG_bulk[i, tissue], sd_mCG_bulk[i, tissue], f = 1 / 10)
})
fits_pos <- lapply(tissues, function(tissue) {
  i <- sample(nrow(mean_mCG_pos), 10 ^ 6)
  lowess(mean_mCG_pos[i, tissue], sd_mCG_pos[i, tissue], f = 1 / 10)
})
fits_neg <- lapply(tissues, function(tissue) {
  i <- sample(nrow(mean_mCG_neg), 10 ^ 6)
  lowess(mean_mCG_neg[i, tissue], sd_mCG_neg[i, tissue], f = 1 / 10)
})
fits <- c(fits_bulk, fits_pos, fits_neg)

approx_fits <- lapply(fits, approxfun)

# ==============================================================================
# Plot
#

tissue_col <- c("BA24" = "deeppink",
                "BA9" = "deepskyblue",
                "HC" = "darkgrey",
                "NAcc" = "chocolate1")
neun_col <- c("bulk" = "black",
              "neg" = "purple",
              "pos" = "darkgreen")

Cairo::CairoPNG("../figures/mCG.smooth_scatter.png",
                width = 960,
                height = 960)
par(mfrow = c(4, 4))
# NOTE: Downsample points to reduce file size
i <- sample(nrow(mean_mCG_bulk), 10 ^ 5)
spwt(x = mean_mCG_bulk[i, 1],
     y = sd_mCG_bulk[i, 1],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2,
     lty = 3)
spwt(x = mean_mCG_bulk[i, 2],
     y = sd_mCG_bulk[i, 2],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 3)
spwt(x = mean_mCG_bulk[i, 3],
     y = sd_mCG_bulk[i, 3],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 3)
spwt(x = mean_mCG_bulk[i, 4],
     y = sd_mCG_bulk[i, 4],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 3)

i <- sample(nrow(mean_mCG_pos), 10 ^ 5)
spwt(x = mean_mCG_pos[i, 1],
     y = sd_mCG_pos[i, 1],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2,
     lty = 1)
spwt(x = mean_mCG_pos[i, 2],
     y = sd_mCG_pos[i, 2],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 1)
spwt(x = mean_mCG_pos[i, 3],
     y = sd_mCG_pos[i, 3],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 1)
spwt(x = mean_mCG_pos[i, 4],
     y = sd_mCG_pos[i, 4],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 1)
spwt(x = mean_mCG_neg[i, 1],
     y = sd_mCG_neg[i, 1],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2,
     lty = 2)
spwt(x = mean_mCG_neg[i, 2],
     y = sd_mCG_neg[i, 2],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 2)
spwt(x = mean_mCG_neg[i, 3],
     y = sd_mCG_neg[i, 3],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 2)
spwt(x = mean_mCG_neg[i, 4],
     y = sd_mCG_neg[i, 4],
     xlim = c(0, 1),
     ylim = c(0, 0.11))
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 2)
# Or all on the same plot
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     xlim = c(0, 1),
     ylim = c(0, 0.11),
     lwd = 2,
     lty = 3)
plot(function(x) approx_fits[[2]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 3)
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 3)
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 3)
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2,
     lty = 1)
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 1)
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 1)
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 1)
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[1]],
     add = TRUE,
     lwd = 2,
     lty = 2)
plot(function(x) approx_fits[[1]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[2]],
     add = TRUE,
     lwd = 2,
     lty = 2)
plot(function(x) approx_fits[[3]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[3]],
     add = TRUE,
     lwd = 2,
     lty = 2)
plot(function(x) approx_fits[[4]](x),
     from = 0,
     to = 1,
     col = tissue_col[names(approx_fits)[4]],
     add = TRUE,
     lwd = 2,
     lty = 2)
dev.off()

pdf("../figures/mCG.mean_var_loess.Tissue_facet.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
par(mfrow = c(2, 2))
# BA9
plot(function(x) approx_fits[[1]](x), from = 0, to = 1,
     col = neun_col["bulk"],
     lwd = 2,
     xlim = c(0, 1),
     ylim = c(0, 0.1),
     xlab = "mean(mCG)",
     ylab = "sd(mCG)",
     main = "BA9")
plot(function(x) approx_fits[[5]](x), from = 0, to = 1,
     col = neun_col["pos"],
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[9]](x), from = 0, to = 1,
     col = neun_col["neg"],
     lwd = 2,
     add = TRUE)
# BA24
plot(function(x) approx_fits[[2]](x), from = 0, to = 1,
     col = neun_col["bulk"],
     lwd = 2,
     xlim = c(0, 1),
     ylim = c(0, 0.1),
     xlab = "mean(mCG)",
     ylab = "sd(mCG)",
     main = "BA24")
plot(function(x) approx_fits[[6]](x), from = 0, to = 1,
     col = neun_col["pos"],
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[10]](x), from = 0, to = 1,
     col = neun_col["neg"],
     lwd = 2,
     add = TRUE)
# HC
plot(function(x) approx_fits[[3]](x), from = 0, to = 1,
     col = neun_col["bulk"],
     lwd = 2,
     xlim = c(0, 1),
     ylim = c(0, 0.1),
     xlab = "mean(mCG)",
     ylab = "sd(mCG)",
     main = "HC")
plot(function(x) approx_fits[[7]](x), from = 0, to = 1,
     col = neun_col["pos"],
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[11]](x), from = 0, to = 1,
     col = neun_col["neg"],
     lwd = 2,
     add = TRUE)
# NAcc
plot(function(x) approx_fits[[4]](x), from = 0, to = 1,
     col = neun_col["bulk"],
     lwd = 2,
     xlim = c(0, 1),
     ylim = c(0, 0.1),
     xlab = "mean(mCG)",
     ylab = "sd(mCG)",
     main = "NAcc")
plot(function(x) approx_fits[[8]](x), from = 0, to = 1,
     col = neun_col["pos"],
     lwd = 2,
     add = TRUE)
plot(function(x) approx_fits[[12]](x), from = 0, to = 1,
     col = neun_col["neg"],
     lwd = 2,
     add = TRUE)
dev.off()

# ------------------------------------------------------------------------------
# Plot proportion of meanDiffs > 0.1 for each pairwise comparison in each NeuN
#

comparisons <- list(c("NAcc", "BA9"),
                    c("NAcc", "BA24"),
                    c("NAcc", "HC"),
                    c("HC", "BA9"),
                    c("HC", "BA24"),
                    c("BA24", "BA9"))
meanDiff_df <- do.call(
  rbind, lapply(comparisons, function(comparison, cutoff = 0.1) {
    bulk <- mean_mCG_bulk[, comparison[1]] - mean_mCG_bulk[, comparison[2]]
    pos <- mean_mCG_pos[, comparison[1]] - mean_mCG_pos[, comparison[2]]
    neg <- mean_mCG_neg[, comparison[1]] - mean_mCG_neg[, comparison[2]]
    data.frame(Comparison = paste0(comparison[1], " vs.", comparison[2]),
               NeuN = c("bulk", "pos", "neg"),
               Prop_CpGs = c(sum(abs(bulk) > cutoff) / length(bulk),
                             sum(abs(pos) > cutoff) / length(pos),
                             sum(abs(neg) > cutoff) / length(neg)))
  })
)

g <- ggplot(meanDiff_df,
            aes(x = Comparison, y = Prop_CpGs, fill = NeuN)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = neun_col) +
  theme_bw() +
  ylab("Proportion CpGs with |meanDiff| > 0.1")
ggsave("../figures/mCG.proportion_CpGs_with_large_meanDiff.pdf", g)



