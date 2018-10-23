# Correlation of bigDARs with RNA-seq profiles (adapted from
# https://f1000research.com/articles/6-2055/v1)
# Peter Hickey
# 2018-01-30

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(limma)

load("../objects/assays-and-features.rda")
elist <- readRDS("../../RNA-seq/objects/elist_with_sv.rds")
big_dars_pos <- dars_pos[abs(dars_pos$logFC) > 1]
rna_atac_meth <- mutate(rna_atac_meth,
       CGI_promoter = gene %in% unlist(cgi_promoters$GENEID))

### ----------------------------------------------------------------------------
### Set up xlim and ylim
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC, CGI_promoter)
lfc$logFC <- NA_real_

ol <- findOverlaps(big_dars_pos,
                   reduce(c(unlist(promoters_by_gene[lfc$gene]),
                           gencode_features$genes[lfc$gene],
                           unlist(fantom5_enhancers_by_gene_all_genes[lfc$gene]),
                           ignore.mcols = TRUE),
                   ignore.strand = TRUE))
logFCs <- tapply(
  X = big_dars_pos$logFC[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$logFC[as.numeric(names(logFCs))] <- logFCs

xlim <- range(lfc$logFC, na.rm = TRUE)
ylim <- range(lfc$expLogFC, na.rm = TRUE)

### ----------------------------------------------------------------------------
### Promoters
###

lfc <- rna_atac_meth %>%
  dplyr::select(gene, gene_symbol, DE, expLogFC, CGI_promoter)
lfc$logFC <- NA_real_

# NOTE: If multiple hits, take the one with the biggest logFC
ol <- findOverlaps(big_dars_pos, promoters_by_gene[lfc$gene],
                   ignore.strand = TRUE)
logFCs <- tapply(
  X = big_dars_pos$logFC[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$logFC[as.numeric(names(logFCs))] <- logFCs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. logFC(ATAC)
#

pdf("../figures/logFC_vs_bigDARs.promoters.pdf",
    height = 7,
    width = 7)
par(mfrow = c(1, 1))
plot(expLogFC ~ logFC,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$logFC, lfc$expLogFC,
                             use = "complete.obs"), 2),
                   "\n"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)

par(mfrow = c(1, 2))
lfc_cgi <- lfc[lfc$CGI_promoter, ]
plot(expLogFC ~ logFC,
     data = lfc_cgi,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_cgi$logFC, lfc_cgi$expLogFC,
                             use = "complete.obs"), 2),
                   "\nCGI promoter"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_cgi$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_cgi))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_cgi$expLogFC ~ 0 + lfc_cgi$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
lfc_no_cgi <- lfc[!lfc$CGI_promoter, ]
plot(expLogFC ~ logFC,
     data = lfc_no_cgi,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_no_cgi$logFC, lfc_no_cgi$expLogFC,
                             use = "complete.obs"), 2),
                   "\nNo CGI promoter"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_no_cgi$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_no_cgi))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_no_cgi$expLogFC ~ 0 + lfc_no_cgi$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
par(mfrow = c(1, 1))
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ logFC,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$logFC, lfc_degs$expLogFC,
                             use = "complete.obs"), 2),
                   "\nCGI promoter"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
par(mfrow = c(1, 2))
lfc_degs_cgi <- lfc_degs[lfc_degs$CGI_promoter, ]
plot(expLogFC ~ logFC,
     data = lfc_degs_cgi,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs_cgi$logFC, lfc_degs_cgi$expLogFC,
                             use = "complete.obs"), 2),
                   "\nCGI promoter"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs_cgi$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs_cgi))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs_cgi$expLogFC ~ 0 + lfc_degs_cgi$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
lfc_degs_no_cgi <- lfc_degs[!lfc_degs$CGI_promoter, ]
plot(expLogFC ~ logFC,
     data = lfc_degs_no_cgi,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs_no_cgi$logFC, lfc_degs_no_cgi$expLogFC,
                             use = "complete.obs"), 2),
                   "\nNo CGI promoter"),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs_no_cgi$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs_no_cgi))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs_no_cgi$expLogFC ~ 0 + lfc_degs_no_cgi$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$logFC)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.bigDARs.promoters.pdf",
    height = 7,
    width = 7)
gw <- lfc$logFC
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = "bigDAR",
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("bottomright",
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
lfc$logFC <- NA_real_
# NOTE: If multiple hits, take the one with the biggest logFC
ol <- findOverlaps(big_dars_pos, gencode_features$genes[lfc$gene],
                   ignore.strand = TRUE)
logFCs <- tapply(
  X = big_dars_pos$logFC[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$logFC[as.numeric(names(logFCs))] <- logFCs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. logFC(ATAC)
#

pdf("../figures/logFC_vs_bigDARs.gene_body.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ logFC,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$logFC, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ logFC,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$logFC, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$logFC)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.bigDARs.gene_body.pdf",
    height = 7,
    width = 7)
gw <- lfc$logFC
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = "bigDAR",
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("bottomright",
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
lfc$logFC <- NA_real_
# NOTE: If multiple hits, take the one with the biggest logFC
ol <- findOverlaps(big_dars_pos, gene_excluding_promoters[lfc$gene],
                   ignore.strand = TRUE)
logFCs <- tapply(
  X = big_dars_pos$logFC[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$logFC[as.numeric(names(logFCs))] <- logFCs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. logFC(ATAC)
#

pdf("../figures/logFC_vs_bigDARs.gene_body_excluding_promoters.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ logFC,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$logFC, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ logFC,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$logFC, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$logFC)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.bigDARs.gene_body_excluding_promoters.pdf",
    height = 7,
    width = 7)
gw <- lfc$logFC
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = "DAR",
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("bottomright",
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
lfc$logFC <- NA_real_
# NOTE: If multiple hits, take the one with the biggest logFC
ol <- findOverlaps(big_dars_pos, fantom5_enhancers_by_gene_all_genes[lfc$gene],
                   ignore.strand = TRUE)
logFCs <- tapply(
  X = big_dars_pos$logFC[queryHits(ol)],
  INDEX = subjectHits(ol),
  FUN = function(x) x[which.max(abs(x))])
lfc$logFC[as.numeric(names(logFCs))] <- logFCs

# -----------------------------------------------------------------------------
# Scatterplots of logFC(RNA) vs. logFC(ATAC)
#

pdf("../figures/logFC_vs_bigDARs.FANTOM5_enhancers.pdf",
    height = 7,
    width = 7)
plot(expLogFC ~ logFC,
     data = lfc,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc$logFC, lfc$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$expLogFC ~ 0 + lfc$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# NOTE: Only looking at DEGs
lfc_degs <- lfc[lfc$DE, ]
plot(expLogFC ~ logFC,
     data = lfc_degs,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(lfc_degs$logFC, lfc_degs$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "log-FC bigDAR",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(lfc_degs$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(lfc_degs))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc_degs$expLogFC ~ 0 + lfc_degs$logFC)
coef(summary(u))
abline(u, col = "red", lwd = 2)
dev.off()

# ------------------------------------------------------------------------------
# limma::fry()
#

fry_df <- data.frame(gene = lfc$gene, weights = lfc$logFC)
fry(y = elist,
    index = fry_df[!is.na(fry_df$weights), ],
    design = elist$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = elist$design))

# ------------------------------------------------------------------------------
# limma::barcodeplot()
#

pdf("../figures/barcodeplot.bigDARs.FANTOM5_enhancers.pdf",
    height = 7,
    width = 7)
gw <- lfc$logFC
gw[is.na(gw)] <- 0
barcodeplot(lfc$expLogFC,
            gene.weights = gw,
            weights.label = "bigDAR",
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos",
            xlab = "log-FC Gene Expression")
legend("bottomright",
       col = c("red", "blue"),
       lty = 1,
       lwd = 2,
       legend = c("NAcc_pos hypermethylated", "NAcc_pos hypomethylated"))
dev.off()
