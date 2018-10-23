# Correlation of CH-DMRs with RNA-seq profiles (adapted from
# https://f1000research.com/articles/6-2055/v1)
# Peter Hickey
# 2018-01-29

# NOTE: Using mCA (+) as representative mCH
# NOTE: A CH-DMR may not necessarily be differential between NAcc_pos and
#       BA9_pos but it's a very safe bet

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(limma)

load("../objects/assays-and-features.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_candidate_CH_DMRs <- lapply(
  X = list_of_candidate_CH_DMRs,
  FUN = function(x) x[x$fwer <= 50])
CA_pos_DMRs <- list_of_candidate_CH_DMRs[["mCA (+)"]]
elist <- readRDS("../../RNA-seq/objects/elist_with_sv.rds")

### ----------------------------------------------------------------------------
### Set up xlim and ylim
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC)
lfc$meanDiff <- NA_real_

ol <- findOverlaps(CA_pos_DMRs,
                   reduce(c(unlist(promoters_by_gene[lfc$gene]),
                            gencode_features$genes[lfc$gene],
                            unlist(fantom5_enhancers_by_gene_all_genes[lfc$gene]),
                            ignore.mcols = TRUE),
                          ignore.strand = TRUE))
meanDiffs <- tapply(
  X = (CA_pos_DMRs$NA_pos - CA_pos_DMRs$BA9_pos)[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$meanDiff[as.numeric(names(meanDiffs))] <- meanDiffs

xlim <- range(lfc$meanDiff, na.rm = TRUE)
ylim <- range(lfc$expLogFC, na.rm = TRUE)

### ----------------------------------------------------------------------------
### Promoters
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC)
lfc$meanDiff <- NA_real_
ol <- findOverlaps(CA_pos_DMRs, promoters_by_gene[lfc$gene],
                   ignore.strand = TRUE)
# NOTE: If multiple hits, take the one with the biggest meanDiff
meanDiffs <- tapply(
  X = (CA_pos_DMRs$NA_pos - CA_pos_DMRs$BA9_pos)[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$meanDiff[as.numeric(names(meanDiffs))] <- meanDiffs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. meanDiff(mCH)
#

pdf("../figures/logFC_vs_CH-DMRs.promoters.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ meanDiff,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$meanDiff, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ meanDiff,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$meanDiff, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$meanDiff)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.CH-DMRs.promoters.pdf",
    height = 7,
    width = 7)
gw <- lfc$meanDiff
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = expression(paste(Delta, "mCH", sep = "")),
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("topright",
       col = c("red", "blue"),
       lty = 1,
       lwd = 2,
       legend = c("NAcc_pos hypermethylated", "NAcc_pos hypomethylated"))
dev.off()

### ----------------------------------------------------------------------------
### Gene body
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC)
lfc$meanDiff <- NA_real_
ol <- findOverlaps(CA_pos_DMRs, gencode_features$genes[lfc$gene],
                   ignore.strand = TRUE)
# NOTE: If multiple hits, take the one with the biggest meanDiff
meanDiffs <- tapply(
  X = (CA_pos_DMRs$NA_pos - CA_pos_DMRs$BA9_pos)[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$meanDiff[as.numeric(names(meanDiffs))] <- meanDiffs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. meanDiff(mCH)
#

pdf("../figures/logFC_vs_CH-DMRs.gene_body.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ meanDiff,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$meanDiff, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ meanDiff,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$meanDiff, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$meanDiff)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.CH-DMRs.gene_body.pdf",
    height = 7,
    width = 7)
gw <- lfc$meanDiff
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = expression(paste(Delta, "mCH", sep = "")),
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("topright",
       col = c("red", "blue"),
       lty = 1,
       lwd = 2,
       legend = c("NAcc_pos hypermethylated", "NAcc_pos hypomethylated"))
dev.off()

### ----------------------------------------------------------------------------
### Gene body excluding promoters
###

gene_excluding_promoters <- psetdiff(
  gencode_features$genes,
  promoters_by_gene[names(gencode_features$genes)])
lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC)
lfc$meanDiff <- NA_real_
ol <- findOverlaps(CA_pos_DMRs, gene_excluding_promoters[lfc$gene],
                   ignore.strand = TRUE)
# NOTE: If multiple hits, take the one with the biggest meanDiff
meanDiffs <- tapply(
  X = (CA_pos_DMRs$NA_pos - CA_pos_DMRs$BA9_pos)[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$meanDiff[as.numeric(names(meanDiffs))] <- meanDiffs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. meanDiff(mCH)
#

pdf("../figures/logFC_vs_CH-DMRs.gene_body_excluding_promoter.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ meanDiff,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$meanDiff, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ meanDiff,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$meanDiff, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$meanDiff)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

# NOTE: Need an object where for each gene there is *one* meanDiff(mCH)

# NOTE: CH-DMR-based barcodeplot
pdf("../figures/barcodeplot.CH-DMRs.gene_body_excluding_promoter.pdf",
    height = 7,
    width = 7)
gw <- lfc$meanDiff
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = expression(paste(Delta, "mCH", sep = "")),
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("topright",
       col = c("red", "blue"),
       lty = 1,
       lwd = 2,
       legend = c("NAcc_pos hypermethylated", "NAcc_pos hypomethylated"))
dev.off()

### ----------------------------------------------------------------------------
### FANTOM5 enhancers
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC)
lfc$meanDiff <- NA_real_
ol <- findOverlaps(CA_pos_DMRs,
                   fantom5_enhancers_by_gene_all_genes[lfc$gene],
                   ignore.strand = TRUE)
# NOTE: If multiple hits, take the one with the biggest meanDiff
meanDiffs <- tapply(
  X = (CA_pos_DMRs$NA_pos - CA_pos_DMRs$BA9_pos)[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$meanDiff[as.numeric(names(meanDiffs))] <- meanDiffs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. meanDiff(mcG)
#

pdf("../figures/logFC_vs_CH-DMRs.FANTOM5_enhancers.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ meanDiff,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$meanDiff, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ meanDiff,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$meanDiff, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "meanDiff(mCH)",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$meanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$meanDiff)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.CH-DMRs.FANTOM5_enhancers.pdf",
    height = 7,
    width = 7)
gw <- lfc$meanDiff
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = expression(paste(Delta, "mCH", sep = "")),
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("topright",
       col = c("red", "blue"),
       lty = 1,
       lwd = 2,
       legend = c("NAcc_pos hypermethylated", "NAcc_pos hypomethylated"))
dev.off()
