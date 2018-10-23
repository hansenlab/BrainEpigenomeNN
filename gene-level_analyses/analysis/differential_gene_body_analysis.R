# Differential analysis of gene body mCA and mCG
# Peter Hickey
# 2018-02-28

### ----------------------------------------------------------------------------
### Setup
###

library(limma)
library(SummarizedExperiment)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

se <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses",  "objects",
            "mC_at_gene_level.rds"))

rna_DE <- read.csv(
  "../../RNA-seq/extdata/topTable.NA_posvsBA9_pos.RNA-seq.csv.gz")
rownames(rna_DE) <- rna_DE$X
rna_DE$X <- NULL
y_rna <- readRDS("../../RNA-seq/objects/elist_with_sv.rds")

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")

### ----------------------------------------------------------------------------
### EDA
###

assay_name <- "mCA"

# Create variables for ID and group
id <- factor(se$Individual)
group <- factor(paste0(se$Tissue, "_", se$NeuN),
                levels = c("BA9_pos", "BA24_pos", "HC_pos", "NA_pos"))
col <- se$Tissue_color

# Extract beta values and M-values
# NOTE: Filter out genes with NA values
beta <- assay(se, assay_name, withDimnames = TRUE)
keep <- !rowAnys(beta, value = NA)
beta <- beta[keep, ]
mval <- log2(beta / (1 - beta))

# Compute quantile normalize beta and M-values
# NOTE: Don't create an EList out of beta-values because this will take the
#       log2-transform!
beta_qn <- normalizeQuantiles(beta)
mval_qn <- normalizeQuantiles(mval)

# Density plots
par(mfrow = c(2, 2))
plotDensities(beta, group = group, col = col,
              legend = ifelse(assay_name == "mCG", "topleft", "topright"),
              main = "Raw beta")
plotDensities(beta_qn, group = group, col = col,
              legend = ifelse(assay_name == "mCG", "topleft", "topright"),
              main = "QN beta")
plotDensities(mval, group = group, col = col, legend = "topleft",
              main = "Raw M-value")
plotDensities(mval_qn, group = group, col = col, legend = "topleft",
              main = "QN M-value")

# MDS plots
# NOTE: QN seems to help for mCA but not so much for mCG
par(mfrow = c(2, 2))
plotMDS(beta, group = group, col = col, legend = "topleft",
        main = "Raw beta")
plotMDS(beta_qn, group = group, col = col, legend = "topleft",
        main = "QN beta")
plotMDS(mval, group = group, col = col, legend = "topleft",
        main = "Raw M-value")
plotMDS(mval_qn, group = group, col = col, legend = "topleft",
        main = "QN M-value")

# Construct design matrices
# NOTE: Blocking on Individual and treating Tissue_pos as 'Treatment'
# NOTE: Any gain in power by using duplicateCorrelation over blocking is small
#       (https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html)
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
blocked_design <- model.matrix(~ id + group )
colnames(blocked_design) <- gsub("group", "", colnames(blocked_design))

# Construct contrasts
contrasts <- makeContrasts(NA_posvsBA9_pos = "NA_pos - BA9_pos",
                           BA24_posvsBA9_pos = "BA24_pos - BA9_pos",
                           HC_posvsBA9_pos = "HC_pos - BA9_pos",
                           levels = design)
blocked_contrasts <- makeContrasts(NA_posvsBA9_pos = "NA_pos",
                                   BA24_posvsBA9_pos = "BA24_pos",
                                   HC_posvsBA9_pos = "HC_pos",
                                   levels = blocked_design)

# Fit linear models and apply contrasts
beta_lm <- contrasts.fit(lmFit(beta, design), contrasts)
blocked_beta_lm <- contrasts.fit(lmFit(beta, blocked_design),
                                 blocked_contrasts)
beta_qn_lm <- contrasts.fit(lmFit(beta_qn, design),
                            contrasts)
blocked_beta_qn_lm <- contrasts.fit(lmFit(beta_qn, blocked_design),
                                    blocked_contrasts)
mval_lm <- contrasts.fit(lmFit(mval, design),
                         contrasts)
blocked_mval_lm <- contrasts.fit(lmFit(mval, blocked_design),
                                 blocked_contrasts)
mval_qn_lm <- contrasts.fit(lmFit(mval_qn, design),
                            contrasts)
blocked_mval_qn_lm <- contrasts.fit(lmFit(mval_qn, blocked_design),
                                    blocked_contrasts)

# Summarise results without use of eBayes()
# Summary:
#   - Choice of beta or M-value doesn't make a huge difference to number of DMGs
#   - Quantile normalization increases the number of DMGs (NA_pos vs BA9_pos)
#     and makes a more 50:50 hypo:hyper ratio
#   - Blocking slightly increases the number of DMGs
summary(decideTests(beta_lm))
summary(decideTests(beta_qn_lm))
summary(decideTests(blocked_beta_lm))
summary(decideTests(blocked_beta_qn_lm))
summary(decideTests(mval_lm))
summary(decideTests(mval_qn_lm))
summary(decideTests(blocked_mval_lm))
summary(decideTests(blocked_mval_qn_lm))

# Run eByaes()
beta_ebayes <- eBayes(beta_lm)
blocked_beta_ebayes <- eBayes(blocked_beta_lm)
beta_qn_ebayes <- eBayes(beta_qn_lm)
blocked_beta_qn_ebayes <- eBayes(blocked_beta_qn_lm)
mval_ebayes <- eBayes(mval_lm)
blocked_mval_ebayes <- eBayes(blocked_mval_lm)
mval_qn_ebayes <- eBayes(mval_qn_lm)
blocked_mval_qn_ebayes <- eBayes(blocked_mval_qn_lm)

# Plots of beta values for top DM genes
top <- topTable(beta_ebayes, "NA_posvsBA9_pos")
genes <- rownames(top)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  stripchart(beta[rownames(beta) == genes[i], ] ~ group,
             method = "jitter",
             group.names = levels(group),
             pch = 16,
             cex = 1.5,
             col = se$Tissue_color,
             ylab = "Beta values",
             vertical = TRUE,
             cex.axis = 1.5,
             cex.lab = 1.5,
             ylim = c(0, 1))
  title(genes[i], cex.main = 1.5)
}

top <- topTable(mval_qn_ebayes, "NA_posvsBA9_pos")
genes <- rownames(top)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  stripchart(beta[rownames(beta) == genes[i], ] ~ group,
             method = "jitter",
             group.names = levels(group),
             pch = 16,
             cex = 1.5,
             col = se$Tissue_color,
             ylab = "Beta values",
             vertical = TRUE,
             cex.axis = 1.5,
             cex.lab = 1.5,
             ylim = c(0, 1))
  title(genes[i], cex.main = 1.5)
}

# Treat analysis
# NOTE: lfcs correspond to meanDiff when analysing beta values
meanDiffs <- c(seq(0, 0.04, 0.01), seq(0.05, 0.3, 0.05))
names(meanDiffs) <- meanDiffs
sapply(meanDiffs, function(meanDiff) {
  nrow(topTreat(treat(beta_lm, meanDiff), number = Inf, p.value = 0.05))
})
# mCG
#     0  0.01  0.02  0.03  0.04  0.05   0.1  0.15   0.2  0.25   0.3
# 14118  9986  6891  5033  3763  2852   919   354   145    59    24
# mCA
#     0  0.01  0.02  0.03  0.04  0.05   0.1  0.15   0.2  0.25   0.3
# 17642  7889  2456   948   404   193     0     0     0     0     0

fcs <- seq(1, 4, 0.5)
lfcs <- log2(fcs)
names(lfcs) <- fcs
sapply(lfcs, function(lfc) {
  nrow(topTreat(treat(mval_lm, lfc), number = Inf, p.value = 0.05))
})
# mCG
#     1   1.5     2   2.5     3   3.5     4
# 14478  2140   645   258   113    66    35
# mCA
#     1   1.5     2   2.5     3   3.5     4
# 18426  1187   199    53    14     5     1

# NOTE: lfc chosen to have approximately 1000 genes for beta values and
#       the same log(M-value) cutoff
beta_treat <- treat(beta_lm,
                    lfc = ifelse(assay_name == "mCG", 0.1, 0.03))
beta_qn_treat <- treat(beta_qn_lm,
                       lfc = ifelse(assay_name == "mCG", 0.1, 0.03))
mval_treat <- treat(mval_lm,
                    lfc = log2(1.5))

par(mfrow = c(1, 2))
plotMD(beta_treat,
       column = "NA_posvsBA9_pos",
       status = topTreat(beta_treat, "NA_posvsBA9_pos",
                         number = Inf, sort.by = "none")$P.Value < 0.05,
       ylab = "meanDiff",
       xlab = "Ave mC")
abline(h = -ifelse(assay_name == "mCG", 0.1, 0.03), lty = 2)
abline(h = ifelse(assay_name == "mCG", 0.1, 0.03), lty = 2)

plotMD(mval_treat,
       column = "NA_posvsBA9_pos",
       status = topTreat(mval_treat, "NA_posvsBA9_pos",
                         number = Inf, sort.by = "none")$P.Value < 0.05)
abline(h = -log2(1.5), lty = 2)
abline(h = log2(1.5), lty = 2)

# Correlation of gene body mC and RNA
par(mfrow = c(2, 3))

tt <- topTable(beta_ebayes, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "meanDiff",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("Beta: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

tt <- topTable(beta_ebayes, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
tt <- tt[abs(tt$logFC) > ifelse(assay_name == "mCG", 0.1, 0.03), ]
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "meanDiff",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("Beta cutoff: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

tt <- topTreat(beta_treat, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "meanDiff",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("Beta treat: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

tt <- topTable(mval_ebayes, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "log-FC Methylation",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("M-val: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

tt <- topTable(mval_ebayes, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
tt <- tt[abs(tt$logFC) > log2(1.5), ]
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "log-FC Methylation",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("M-val cutoff: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

tt <- topTreat(mval_treat, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
plot(lfc$ME, lfc$RNA,
     xlab = "log-FC Methylation",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = "gray30",
     main = paste0("M-val treat: cor = ",
                   round(cor(lfc$ME, lfc$RNA), 2)),
     sub = paste0("n = ", sum(complete.cases(lfc)), " / ", nrow(tt)))
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(lfc$RNA ~ 0 + lfc$ME)
abline(u, col = "red", lwd = 2)

# Gene set testing
tt <- topTreat(beta_treat, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(tt), rownames(rna_DE))
lfc <- data.frame(
  gene = rownames(tt),
  ME = tt[, "logFC"],
  stringsAsFactors = FALSE)
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
ME <- data.frame(GeneID = lfc$gene, weights = lfc$ME)
fry(y_rna,
    index = ME,
    design = y_rna$design,
    contrast = makeContrasts(NA_pos - BA9_pos, levels = y_rna$design))

# Barcode plot
tt <- topTreat(beta_treat, "NA_posvsBA9_pos", n = Inf, p.value = 0.05)
m <- match(row.names(rna_DE), row.names(tt))
gw <- tt$logFC[m]
gw[is.na(gw)] <- 0
par(mfrow = c(1, 2))
barcodeplot(rna_DE$logFC,
            gene.weights = gw,
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos")
legend("topright",
       col = c("red","blue"),
       lty = 1,
       lwd = 2,
       legend = c("Up-m in NAcc_pos", "Up-m in BA9_pos"))

tt <- rna_DE[rna_DE$adj.P.Val < 0.05, ]
m <- match(row.names(beta_treat), row.names(tt))
gw <- tt$logFC[m]
gw[is.na(gw)] <- 0
barcodeplot(topTreat(beta_treat,
                     "NA_posvsBA9_pos", n = Inf,
                     p.value = 1)$logFC,
            gene.weights = gw,
            labels = c("BA9_pos","NAcc_pos"),
            main = "NAcc_pos vs BA9_pos")
legend("topright",
       col = c("red","blue"),
       lty = 1,
       lwd = 2,
       legend = c("Up-r in NAcc_pos", "Up-r in BA9_pos"))

# How are DMGs that are also DEGs different from those that aren't?
dmg_names <- rownames(
  topTreat(beta_treat, "NA_posvsBA9_pos", n = Inf, p.value = 0.05))
fisher.test(table(dmg_names %in% deg_names,
                  dmg_names %in% names(cgi_promoters_by_gene)))
fisher.test(table(rownames(se) %in% names(cgi_promoters_by_gene),
                  rownames(se) %in% deg_names))

k <- dmg_names %in% deg_names
boxplot(width(rowRanges(se)[dmg_names]) ~ k)

fisher.test(table(dmg_names %in% deg_names,
                  dmg_names %in% pc_genes))
fisher.test(table(rownames(se) %in% pc_genes,
                  rownames(se) %in% deg_names))

fisher.test(table(dmg_names %in% deg_names,
                  dmg_names %in% lnc_genes))
fisher.test(table(rownames(se) %in% lnc_genes,
                  rownames(se) %in% deg_names))

### ----------------------------------------------------------------------------
### TODOs
###

# TODO: Add gene body ATAC-seq?

### ----------------------------------------------------------------------------
### NOTEs
###

# - [x] Summarise beta vs. m-value
#   - mCG
#     - Beta values better than M-value based on MDS plot
#     - Similar correlation from treat (-0.41 for beta, -0.38 for M)
#   - mCA
#     - M-values give visually better separation of tissues
#     - Similar correlation without treat (-0.49 for beta, -0.45 for M)
# - [x] Raw vs. QN
#   - mCG
#     - Minor tidy up of beta value MDS plot, bigger tidy up of M-value MDS plot
#   - mCA
#     - QN tightens tissue-clusters in both beta value and M-value MDS plot
# - [x] Argue for use of treat
#   - mCG
#     - 43% of genes are DM otherwise
#     - Treat improves correlation with logFC
#   - mCA
#     - 53% of genes are DM otherwise
#     - Treat improves correlation with logFC
# - [x] Barcode plots
#   - mCG
#     - DMGs have decent enrichment for DEGs but DEGs no enrichment for DMGs
#   - mCA
#     - DMGs have decent enrichment for DEGs but DEGs no enrichment for DMGs
#
# # mCG
# > summary(decideTests(beta_lm))
# NA_posvsBA9_pos BA24_posvsBA9_pos HC_posvsBA9_pos
# Down              2571                37            2794
# NotSig           18834             30297           28073
# Up               11547              2618            2085
#
#
# # meanDiff = 0.1
# > summary(decideTests(beta_treat))
# NA_posvsBA9_pos BA24_posvsBA9_pos HC_posvsBA9_pos
# Down               102                 0               3
# NotSig           32033             32952           32948
# Up                 817                 0               1
#
# # -----------------------------------------------------
#
# # mCA
#
# > summary(decideTests(beta_lm))
# NA_posvsBA9_pos BA24_posvsBA9_pos HC_posvsBA9_pos
# Down             15180                 5            1347
# NotSig           15257             32891           31279
# Up                2462                 3             273
#
#
# # meanDiff = 0.05
# > summary(decideTests(beta_treat))
#        NA_posvsBA9_pos BA24_posvsBA9_pos HC_posvsBA9_pos
# Down               656                 0               1
# NotSig           31951             32899           32896
# Up                 292                 0               2


# fisher.test(table(names(gencode_features$genes) %in% names(cgi_promoters_by_gene), names(gencode_features$genes) %in% deg_names))
