# Plot condition-specific average mCA and mCT in 1kb bins along genome
# Peter Hickey
# 2018-02-19

### ----------------------------------------------------------------------------
### Setup
###

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)

### ----------------------------------------------------------------------------
### Load data
###

se <- readRDS("../objects/average_mCA_and_mCT_in_1kb_bins.rds")
# NOTE: Only using autosomal bins
se <- keepSeqlevels(se, paste0("chr", 1:22), pruning.mode = "coarse")

### ----------------------------------------------------------------------------
### Scatterplots
###

# ------------------------------------------------------------------------------
# meanDiff(mCA+) vs meanDiff(mCT+) for each strand
#

x <- data_frame(
  meanDiff_mCA_pos = assay(se)[, "CpA (+) NAcc_pos"] -
    assay(se)[, "CpA (+) BA9_pos"],
  meanDiff_mCA_neg = assay(se)[, "CpA (-) NAcc_pos"] -
    assay(se)[, "CpA (-) BA9_pos"],
  meanDiff_mCT_pos = assay(se)[, "CpT (+) NAcc_pos"] -
    assay(se)[, "CpT (+) BA9_pos"],
  meanDiff_mCT_neg = assay(se)[, "CpT (-) NAcc_pos"] -
    assay(se)[, "CpT (-) BA9_pos"])

lim <- max(sapply(x, function(xx) max(abs(xx), na.rm = TRUE)))

pdf("~/tmp.pdf", height = 7, width = 14)
par(mfrow = c(2, 2))
smoothScatter(
  x = x$meanDiff_mCA_pos,
  y = x$meanDiff_mCT_pos,
  main = paste0("cor = ",
                round(cor(x$meanDiff_mCA_pos, x$meanDiff_mCT_pos,
                          use = "complete.obs"), 2)),
  xlim = c(-lim, lim),
  ylim = c(-lim, lim))
abline(a = 0, b = 1, col = "red")
lm <- lm(x$meanDiff_mCT_pos ~ 0 + x$meanDiff_mCA_pos)
abline(a = 0, b = coef(lm), col = "blue")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

smoothScatter(
  x = x$meanDiff_mCA_neg,
  y = x$meanDiff_mCT_neg,
  main = paste0("cor = ",
                round(cor(x$meanDiff_mCA_neg, x$meanDiff_mCT_neg,
                          use = "complete.obs"), 2)),
  xlim = c(-lim, lim),
  ylim = c(-lim, lim))
abline(a = 0, b = 1, col = "red")
lm <- lm(x$meanDiff_mCT_neg ~ 0 + x$meanDiff_mCA_neg)
abline(a = 0, b = coef(lm), col = "blue")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

smoothScatter(
  x = x$meanDiff_mCA_pos,
  y = x$meanDiff_mCA_neg,
  main = paste0("cor = ",
                round(cor(x$meanDiff_mCA_pos, x$meanDiff_mCA_neg,
                          use = "complete.obs"), 2)),
  xlim = c(-lim, lim),
  ylim = c(-lim, lim))
abline(a = 0, b = 1, col = "red")
lm <- lm(x$meanDiff_mCA_neg ~ 0 + x$meanDiff_mCA_pos)
abline(a = 0, b = coef(lm), col = "blue")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

smoothScatter(
  x = x$meanDiff_mCT_pos,
  y = x$meanDiff_mCT_neg,
  main = paste0("cor = ",
                round(cor(x$meanDiff_mCT_pos, x$meanDiff_mCT_neg,
                          use = "complete.obs"), 2)),
  xlim = c(-lim, lim),
  ylim = c(-lim, lim))
abline(a = 0, b = 1, col = "red")
lm <- lm(x$meanDiff_mCT_neg ~ 0 + x$meanDiff_mCT_pos)
abline(a = 0, b = coef(lm), col = "blue")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
dev.off()
