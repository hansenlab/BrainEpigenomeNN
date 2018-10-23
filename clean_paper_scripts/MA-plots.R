# MA-plots for paper
# Peter Hickey
# 2017-09-15

library(limma)
library(scales)
rna <- readRDS("../RNA-seq/objects/fit_with_sv.rds")
atac <- readRDS("../ATAC-seq/objects/fit_with_sv.rds")

set.seed(666)

#-------------------------------------------------------------------------------
# RNA-seq
#

pdf("NA_posvsBA9_pos.MD-plot.RNA-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
plotMD(rna,
       column = 1,
       status = ifelse(
         topTable(rna, coef = 1, n = Inf,
                  sort.by = "none")$adj.P.Val < 0.05, "DEG", "non-DEG"),
       hl.col = alpha("orange", 0.3),
       hl.cex = 1,
       bg.col = alpha("black", 0.3),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       xlab = "log2(expression)",
       ylab = "log2(FC)",
       ylim = c(-8, 10))
legend(x = "topright",
       legend = c("non-DEG", "DEG"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()

pdf("NA_negvsBA9_neg.MD-plot.RNA-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
plotMD(rna, column = 2,
       status = ifelse(
         topTable(rna, coef = 2, n = Inf,
                  sort.by = "none")$adj.P.Val < 0.05, "DEG", "non-DEG"),
       hl.col = alpha("orange", 0.3),
       hl.cex = 1,
       bg.col = alpha("black", 0.3),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       xlab = "log2(expression)",
       ylab = "log2(FC)",
       ylim = c(-8, 10))
legend(x = "topright",
       legend = c("non-DEG", "DEG"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()

pdf("ave_pos_vs_ave_neg.MD-plot.RNA-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
plotMD(rna,
       column = 3,
       status = ifelse(
         topTable(rna, coef = 3, n = Inf,
                  sort.by = "none")$adj.P.Val < 0.05, "DEG", "non-DEG"),
       hl.col = alpha("orange", 0.3),
       hl.cex = 1,
       bg.col = alpha("black", 0.3),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       xlab = "log2(expression)",
       ylab = "log2(FC)",
       ylim = c(-8, 10))
legend(x = "topright",
       legend = c("non-DEG", "DEG"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()

#-------------------------------------------------------------------------------
# ATAC-seq
#

pdf("NA_posvsBA9_pos.MD-plot.ATAC-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
# A subset of the ATAC-seq MArrayLM object where a random 20% of 'significant'
# rows and a random 10% of the 'non-significant' rows are used for plotting
P <- topTable(atac, coef = 1, n = Inf, sort.by = "none")$adj.P.Val
logFC <- topTable(atac, coef = 1, n = Inf, sort.by = "none")$logFC
i <- sort(c(sample(which(P < 0.05),
                   floor(sum(P < 0.05) / 5)),
            sample(which(P > 0.05),
                   floor(sum(P > 0.05) / 10))))
plotMD(atac[i, ],
       column = 1,
       status = ifelse(P[i] < 0.05, "DAR", "non-DAR"),
       hl.col = alpha("orange", 0.1),
       hl.cex = 1,
       bg.col = alpha("black", 0.1),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       ylim = c(-3, 3),
       xlab = "log2(abundance)",
       ylab = "log2(FC)")
abline(h = 1, lty = 2)
abline(h = -1, lty = 2)
legend(x = "topright",
       legend = c("non-DAR", "DAR"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()

pdf("NA_negvsBA9_neg.MD-plot.ATAC-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
# A subset of the ATAC-seq MArrayLM object where all 'significant' rows and a
# random 10% of the 'non-significant' rows are used for plotting
# NOTE: No subsampling of DARs since there are so few points
P <- topTable(atac, coef = 2, n = Inf, sort.by = "none")$adj.P.Val
logFC <- topTable(atac, coef = 2, n = Inf, sort.by = "none")$logFC
i <- sort(c(which(P < 0.05),
            sample(which(P > 0.05),
                   floor(sum(P > 0.05) / 10))))
plotMD(atac[i, ],
       column = 2,
       status = ifelse(P[i] < 0.05, "DAR", "non-DAR"),
       # NOTE: No transparency of DARs since there are so few points
       hl.col = alpha("orange", 1),
       hl.cex = 1,
       bg.col = alpha("black", 0.1),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       ylim = c(-3, 3),
       xlab = "log2(abundance)",
       ylab = "log2(FC)")
abline(h = 1, lty = 2)
abline(h = -1, lty = 2)
legend(x = "topright",
       legend = c("non-DAR", "DAR"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()

pdf("ave_pos_vs_ave_neg.MD-plot.ATAC-seq.pdf", width = 5, height = 5,
    useDingbats = FALSE)
# A subset of the ATAC-seq MArrayLM object where all 'significant' rows and a
# random 10% of the 'non-significant' rows are used for plotting
P <- topTable(atac, coef = 3, n = Inf, sort.by = "none")$adj.P.Val
logFC <- topTable(atac, coef = 3, n = Inf, sort.by = "none")$logFC
i <- sort(c(sample(which(P < 0.05),
                   floor(sum(P < 0.05) / 5)),
            sample(which(P > 0.05),
                   floor(sum(P > 0.05) / 10))))
plotMD(atac[i, ],
       column = 3,
       status = ifelse(P[i] < 0.05, "DAR", "non-DAR"),
       values = c("DAR"),
       hl.col = alpha("orange", 0.1),
       hl.cex = 1,
       bg.col = alpha("black", 0.1),
       bg.cex = 0.3,
       legend = FALSE,
       main = "",
       bty = "l",
       ylim = c(-3, 3),
       xlab = "log2(abundance)",
       ylab = "log2(FC)")
abline(h = 1, lty = 2)
abline(h = -1, lty = 2)
legend(x = "topright",
       legend = c("non-DAR", "DAR"),
       pch = c(16, 16),
       col = c("black", "orange"),
       cex = 0.9,
       bty = "n",
       pt.cex = c(1, 1))
dev.off()
