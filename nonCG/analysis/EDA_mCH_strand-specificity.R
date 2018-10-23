# Strand-specificity of mCH
# Peter Hickey
# 2017-08-28

# TODO: Check scatterplot with Lister's data (only after trying the above)

### ============================================================================
### Setup
###

library(bsseq)
library(matrixStats)
library(viridis)
library(gplots)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

options("scipen" = 10)

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
names(CH_BSseq_names) <- paste0("m", contexts, " (",
                                ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq[, grep("pos", colnames(BSseq))]
})
list_of_CH_BSseq <- lapply(list_of_CH_BSseq, function(bsseq) {
  bsseq$col <- bsseq$Tissue_color
  bsseq
})

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")

### ============================================================================
### Per-sample average methylation in bins
###

w <- 100000 / 10 ^ 3
bins <- tileGenome(
  seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0("chr", 1:22)],
  tilewidth = w,
  cut.last.tile.in.chrom = TRUE)
# Add dinucleotide counts of each bin
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, bins)
dns <- c("AA", "AC", "AG", "AT",
         "CA", "CC", "CG", "CT",
         "GA", "GC", "GG", "GT",
         "TA", "TC", "TG", "TT")
names(dns) <- dns
dn_counts <- lapply(X = dns,
                    FUN = vcountPattern,
                    subject = seq)
mcols(bins) <- as(dn_counts, "DataFrame")
# Save bins
saveRDS(bins, paste0("../objects/bins.", w, "kb.rds"))

list_of_meanMeth <- lapply(list_of_CH_BSseq, function(BSseq) {
  getMeth(BSseq = BSseq, regions = bins, type = "smooth",
          what = "perRegion")
})
list_of_meanMeth <- lapply(list_of_meanMeth, as.matrix)
saveRDS(list_of_meanMeth,
        paste0("../objects/list_of_meanMeth.", w, "kb_bins.rds"))

list_of_meanRawMeth <- lapply(list_of_CH_BSseq, function(BSseq) {
  getMeth(BSseq = BSseq, regions = bins, type = "raw",
          what = "perRegion")
})
list_of_meanRawMeth <- lapply(list_of_meanRawMeth, as.matrix)
saveRDS(list_of_meanRawMeth,
        paste0("../objects/list_of_meanRawMeth.", w, "kb_bins.rds"))

### ============================================================================
### Identify bins with large CH strand bias
###

# ------------------------------------------------------------------------------
# Strand-bias by dinucleotide ratios or differences
#

# TODO: Filter out bins with low CA/CT counts

png(paste0("../figures/mCA_strand-specificity.", w, "kb.%d.png"),
    width = 960,
    height = 960)
par(mfrow = c(4, 4))
lapply(dns, function(dn) {
  fwd <- dn
  rev <- as.character(reverseComplement(DNAString(fwd)))
  df <- data.frame(x = log2((mcols(bins)[[fwd]] + 0.5) /
                              (mcols(bins)[[rev]] + 0.5)),
                   y = list_of_meanMeth[["mCA (+)"]][, "5248_BA24_pos"] -
                     list_of_meanMeth[["mCA (-)"]][, "5248_BA24_pos"])
  repete::spwt(y ~ x, df[complete.cases(df), ],
               xlab = paste0("log2(", fwd, "_pos / ", rev, "_neg)"),
               ylab = "mCA_pos - mCA_neg",
               col = scales::alpha("black", 0.2),
               main = paste0(w, " kb bins"))
})
par(mfrow = c(4, 4))
lapply(dns, function(dn) {
  fwd <- dn
  rev <- as.character(reverseComplement(DNAString(fwd)))
  df <- data.frame(x = mcols(bins)[[fwd]] - mcols(bins)[[rev]],
                   y = list_of_meanMeth[["mCA (+)"]][, "5248_BA24_pos"] -
                     list_of_meanMeth[["mCA (-)"]][, "5248_BA24_pos"])
  repete::spwt(y ~ x, df[complete.cases(df), ],
               xlab = paste0(fwd, "_pos - ", rev, "_neg"),
               ylab = "mCA_pos - mCA_neg",
               col = scales::alpha("black", 0.2),
               main = paste0(w, " kb bins"))
})
dev.off()

png(paste0("../figures/mCT_strand-specificity.", w, "kb.%d.png"),
    width = 960,
    height = 960)
par(mfrow = c(4, 4))
lapply(dns, function(dn) {
  fwd <- dn
  rev <- as.character(reverseComplement(DNAString(fwd)))
  df <- data.frame(x = log2((mcols(bins)[[fwd]] + 0.5) /
                              (mcols(bins)[[rev]] + 0.5)),
                   y = list_of_meanMeth[["mCT (+)"]][, "5248_BA24_pos"] -
                     list_of_meanMeth[["mCT (-)"]][, "5248_BA24_pos"])
  repete::spwt(y ~ x, df[complete.cases(df), ],
               xlab = paste0("log2(", fwd, "_pos / ", rev, "_neg)"),
               ylab = "mCT_pos - mCT_neg",
               col = scales::alpha("black", 0.2),
               main = paste0(w, " kb bins"))
})
par(mfrow = c(4, 4))
lapply(dns, function(dn) {
  fwd <- dn
  rev <- as.character(reverseComplement(DNAString(fwd)))
  df <- data.frame(x = mcols(bins)[[fwd]] - mcols(bins)[[rev]],
                   y = list_of_meanMeth[["mCT (+)"]][, "5248_BA24_pos"] -
                     list_of_meanMeth[["mCT (-)"]][, "5248_BA24_pos"])
  repete::spwt(y ~ x, df[complete.cases(df), ],
               xlab = paste0(fwd, "_pos - ", rev, "_neg"),
               ylab = "mCT_pos - mCT_neg",
               col = scales::alpha("black", 0.2),
               main = paste0(w, " kb bins"))
})
dev.off()

# ------------------------------------------------------------------------------
# Strand-bias by genomic context
#

list_of_meanMeth_genes <- lapply(list_of_CH_BSseq, function(BSseq) {
  getMeth(BSseq = BSseq,
          regions = unstrand(unflattened_features$genes),
          type = "smooth",
          what = "perRegion")
})
list_of_meanMeth_genic <- lapply(list_of_meanMeth_genes, as.matrix)
saveRDS(list_of_meanMeth_genic,
        paste0("../objects/list_of_meanMeth_genic.rds"))

intergenic <- gaps(unstrand(unflattened_features$genes))
intergenic <- intergenic[strand(intergenic) == "*"]
list_of_meanMeth_intergenic <-  lapply(list_of_CH_BSseq, function(BSseq) {
  getMeth(BSseq = BSseq,
          regions = intergenic,
          type = "smooth",
          what = "perRegion")
})
list_of_meanMeth_intergenic <- lapply(list_of_meanMeth_intergenic, as.matrix)
saveRDS(list_of_meanMeth_intergenic,
        paste0("../objects/list_of_meanMeth_intergenic.rds"))

# TODO: Compare strand-bias in genic and intergenic regions


# ------------------------------------------------------------------------------
# TODO
#

# TODO: Figure out what the bins with high mCH strand-bias have in common.
#       Good candidates: gene bodies, transposons, sequencing coverage or
#       mappability bias

### ============================================================================
### Correlations of binned mCH
###

# TODO: Should mCG correlation heatmaps use Spearman rather than Pearson?

# ------------------------------------------------------------------------------
# Within context and strand
#

list_of_cor_matrices <- lapply(list_of_meanMeth, function(meth) {
  cor(meth, use = "complete.obs", method = "spearman")
})

hmcol <- viridis(100)

pdf(paste0("../figures/correlation_heatmaps.", w,
           "kb_bins_mean_mCH.within_context_and_strand.pdf"))
lapply(names(list_of_cor_matrices), function(n) {
  cor_matrix <- list_of_cor_matrices[[n]]
  heatmap.2(cor_matrix,
            Rowv = TRUE,
            Colv = TRUE,
            trace = "none",
            col = hmcol,
            margins = c(4, 7),
            srtCol = 30,
            # density.info = "none",
            main = n,
            RowSideColors =
              dplyr::case_when(grepl("5248", colnames(cor_matrix)) ~
                                 "aquamarine",
                               grepl("5284", colnames(cor_matrix)) ~
                                 "chocolate1",
                               grepl("5552", colnames(cor_matrix)) ~
                                 "darkgrey",
                               grepl("5569", colnames(cor_matrix)) ~
                                 "deepskyblue",
                               grepl("5570", colnames(cor_matrix)) ~
                                 "firebrick3",
                               grepl("5628", colnames(cor_matrix)) ~
                                 "purple"),
            ColSideColors =
              dplyr::case_when(grepl("BA24", colnames(cor_matrix)) ~
                                 "deeppink",
                               grepl("BA9", colnames(cor_matrix)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(cor_matrix)) ~
                                 "darkgrey",
                               grepl("NA", colnames(cor_matrix)) ~
                                 "chocolate1"))
})
dev.off()

list_of_cor_matrices_raw <- lapply(list_of_meanRawMeth, function(meth) {
  cor(meth, use = "complete.obs", method = "spearman")
})

w <- median(width(bins)) / 10 ^ 3
pdf(paste0("../figures/correlation_heatmaps.", w,
           "kb_bins_mean_raw_mCH.within_context_and_strand.pdf"))
lapply(names(list_of_cor_matrices_raw), function(n) {
  cor_matrix <- list_of_cor_matrices_raw[[n]]
  heatmap.2(cor_matrix,
            Rowv = TRUE,
            Colv = TRUE,
            trace = "none",
            col = hmcol,
            margins = c(4, 7),
            srtCol = 30,
            # density.info = "none",
            main = paste0(n, " (raw)"),
            RowSideColors =
              dplyr::case_when(grepl("5248", colnames(cor_matrix)) ~
                                 "aquamarine",
                               grepl("5284", colnames(cor_matrix)) ~
                                 "chocolate1",
                               grepl("5552", colnames(cor_matrix)) ~
                                 "darkgrey",
                               grepl("5569", colnames(cor_matrix)) ~
                                 "deepskyblue",
                               grepl("5570", colnames(cor_matrix)) ~
                                 "firebrick3",
                               grepl("5628", colnames(cor_matrix)) ~
                                 "purple"),
            ColSideColors =
              dplyr::case_when(grepl("BA24", colnames(cor_matrix)) ~
                                 "deeppink",
                               grepl("BA9", colnames(cor_matrix)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(cor_matrix)) ~
                                 "darkgrey",
                               grepl("NA", colnames(cor_matrix)) ~
                                 "chocolate1"))
})
dev.off()

# ------------------------------------------------------------------------------
# Across context and strand
#

meanMeth <- do.call(cbind, list_of_meanMeth)
cn <- paste0(colnames(meanMeth), "_",
             rep(names(list_of_meanMeth),
                 times = sapply(list_of_meanMeth, ncol)))
colnames(meanMeth) <- cn

cor_meanMeth <- cor(meanMeth, use = "complete.obs", method = "spearman")

pdf(paste0("../figures/correlation_heatmaps.", w,
           "kb_bins_mean_mCH.across_context_and_strand.pdf"))
heatmap.2(cor_meanMeth,
          Rowv = TRUE,
          Colv = TRUE,
          trace = "none",
          col = hmcol,
          margins = c(5, 8),
          srtCol = 30,
          # density.info = "none",
          main = paste0("meanMeth (", w, " kb bins)"),
          RowSideColors =
            dplyr::case_when(grepl("mCA \\(\\+", cn) ~
                               "#66a61e",
                             grepl("mCA \\(\\-", cn) ~
                               scales::alpha("#66a61e", 0.5),
                             grepl("mCT \\(\\+", cn) ~
                               "#e6ab02",
                             grepl("mCT \\(\\-", cn) ~
                               scales::alpha("#e6ab02", 0.5)),
          ColSideColors =
            dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                             grepl("BA9", cn) ~ "deepskyblue",
                             grepl("HC", cn) ~ "darkgrey",
                             grepl("NA", cn) ~ "chocolate1"))
dev.off()

meanRawMeth <- do.call(cbind, list_of_meanRawMeth)
cn <- paste0(colnames(meanRawMeth), "_",
             rep(names(list_of_meanRawMeth),
                 times = sapply(list_of_meanRawMeth, ncol)))
colnames(meanRawMeth) <- cn

cor_meanRawMeth <- cor(meanRawMeth, use = "complete.obs", method = "spearman")

pdf(paste0("../figures/correlation_heatmaps.", w,
           "kb_bins_mean_raw_mCH.across_context_and_strand.pdf"))
heatmap.2(cor_meanRawMeth,
          Rowv = TRUE,
          Colv = TRUE,
          trace = "none",
          col = hmcol,
          margins = c(5, 8),
          srtCol = 30,
          # density.info = "none",
          main = paste0("meanRawMeth (", w, " kb bins)"),
          RowSideColors =
            dplyr::case_when(grepl("mCA \\(\\+", cn) ~
                               "#66a61e",
                             grepl("mCA \\(\\-", cn) ~
                               scales::alpha("#66a61e", 0.5),
                             grepl("mCT \\(\\+", cn) ~
                               "#e6ab02",
                             grepl("mCT \\(\\-", cn) ~
                               scales::alpha("#e6ab02", 0.5)),
          ColSideColors =
            dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                             grepl("BA9", cn) ~ "deepskyblue",
                             grepl("HC", cn) ~ "darkgrey",
                             grepl("NA", cn) ~ "chocolate1"))
dev.off()

### ============================================================================
### PCA of binned methylation
###

meanMeth <- do.call(cbind, list_of_meanMeth)
cn <- paste0(colnames(meanMeth), "_",
             rep(names(list_of_meanMeth),
                 times = sapply(list_of_meanMeth, ncol)))
colnames(meanMeth) <- cn
meanMeth <- meanMeth[!rowAlls(meanMeth, value = NA), ]
# NOTE: Need to subtract rowMeans() separately for context and strand
meanMeth_row_centered <- meanMeth -
  rowMeans2(meanMeth, cols = grep("mCA \\(\\+\\)", colnames(meanMeth)),
            na.rm = TRUE)
meanMeth_row_centered[, grep("mCA \\(-\\)", colnames(meanMeth))] <-
  meanMeth[, grep("mCA \\(-\\)", colnames(meanMeth))] -
  rowMeans2(meanMeth, cols = grep("mCA \\(-\\)", colnames(meanMeth)),
            na.rm = TRUE)
meanMeth_row_centered[, grep("mCT \\(\\+\\)", colnames(meanMeth))] <-
  meanMeth[, grep("mCT \\(\\+\\)", colnames(meanMeth))] -
  rowMeans2(meanMeth, cols = grep("mCT \\(\\+\\)", colnames(meanMeth)),
            na.rm = TRUE)
meanMeth_row_centered[, grep("mCT \\(-\\)", colnames(meanMeth))] <-
  meanMeth[, grep("mCT \\(-\\)", colnames(meanMeth))] -
  rowMeans2(meanMeth, cols = grep("mCT \\(-\\)", colnames(meanMeth)),
            na.rm = TRUE)
meanMeth_row_centered <-
  meanMeth_row_centered[!matrixStats::rowAnyNAs(meanMeth_row_centered), ]
meanMeth_row_centered_cp <- crossprod(meanMeth_row_centered)
meanMeth_row_centered_cp_svd <- svd(meanMeth_row_centered_cp)
meanMeth_row_centered_cp_pve <- meanMeth_row_centered_cp_svd$d /
  sum(meanMeth_row_centered_cp_svd$d)

pdf(paste0("../figures/PCA.", w, "kb_bins_mean_mCH.pdf"))
plot(meanMeth_row_centered_cp_svd$u[, 1:2],
     col = dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                            grepl("BA9", cn) ~ "deepskyblue",
                            grepl("HC", cn) ~ "darkgrey",
                            grepl("NA", cn) ~ "chocolate1"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC1 (",
                   round(100 * meanMeth_row_centered_cp_pve[1],
                         1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * meanMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     main = paste0("mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"),
       legend = c("BA24", "BA9", "HC", "NAcc"),
       bty = "n")

plot(meanMeth_row_centered_cp_svd$u[, 1:2],
     col = dplyr::case_when(grepl("5248", cn) ~
                              "aquamarine",
                            grepl("5284", cn) ~
                              "chocolate1",
                            grepl("5552", cn) ~
                              "darkgrey",
                            grepl("5569", cn) ~
                              "deepskyblue",
                            grepl("5570", cn) ~
                              "firebrick3",
                            grepl("5628", cn) ~
                              "purple"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC1 (",
                   round(100 * meanMeth_row_centered_cp_pve[1],
                         1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * meanMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     main = paste0("mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("aquamarine", "chocolate1", "darkgrey", "deepskyblue",
                "firebrick3", "purple"),
       legend = c("5248", "5284", "5552", "5569", "5570", "5628"),
       bty = "n")

plot(meanMeth_row_centered_cp_svd$u[, 2:3],
     col = dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                            grepl("BA9", cn) ~ "deepskyblue",
                            grepl("HC", cn) ~ "darkgrey",
                            grepl("NA", cn) ~ "chocolate1"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC2 (",
                   round(100 * meanMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     ylab = paste0("PC3 (",
                   round(100 * meanMeth_row_centered_cp_pve[3],
                         1),
                   "%)"),
     main = paste0("mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"),
       legend = c("BA24", "BA9", "HC", "NAcc"),
       bty = "n")

plot(meanMeth_row_centered_cp_svd$u[, 2:3],
     col = dplyr::case_when(grepl("5248", cn) ~
                              "aquamarine",
                            grepl("5284", cn) ~
                              "chocolate1",
                            grepl("5552", cn) ~
                              "darkgrey",
                            grepl("5569", cn) ~
                              "deepskyblue",
                            grepl("5570", cn) ~
                              "firebrick3",
                            grepl("5628", cn) ~
                              "purple"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC2 (",
                   round(100 * meanMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     ylab = paste0("PC3 (",
                   round(100 * meanMeth_row_centered_cp_pve[3],
                         1),
                   "%)"),
     main = paste0("mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("aquamarine", "chocolate1", "darkgrey", "deepskyblue",
                "firebrick3", "purple"),
       legend = c("5248", "5284", "5552", "5569", "5570", "5628"),
       bty = "n")
dev.off()

meanRawMeth <- do.call(cbind, list_of_meanRawMeth)
cn <- paste0(colnames(meanRawMeth), "_",
             rep(names(list_of_meanRawMeth),
                 times = sapply(list_of_meanRawMeth, ncol)))
colnames(meanRawMeth) <- cn
meanRawMeth <- meanRawMeth[!rowAlls(meanRawMeth, value = NA), ]
# NOTE: Need to subtract rowMeans() separately for context and strand
meanRawMeth_row_centered <- meanRawMeth -
  rowMeans2(meanRawMeth, cols = grep("mCA \\(\\+\\)", colnames(meanRawMeth)),
            na.rm = TRUE)
meanRawMeth_row_centered[, grep("mCA \\(-\\)", colnames(meanRawMeth))] <-
  meanRawMeth[, grep("mCA \\(-\\)", colnames(meanRawMeth))] -
  rowMeans2(meanRawMeth, cols = grep("mCA \\(-\\)", colnames(meanRawMeth)),
            na.rm = TRUE)
meanRawMeth_row_centered[, grep("mCT \\(\\+\\)", colnames(meanRawMeth))] <-
  meanRawMeth[, grep("mCT \\(\\+\\)", colnames(meanRawMeth))] -
  rowMeans2(meanRawMeth, cols = grep("mCT \\(\\+\\)", colnames(meanRawMeth)),
            na.rm = TRUE)
meanRawMeth_row_centered[, grep("mCT \\(-\\)", colnames(meanRawMeth))] <-
  meanRawMeth[, grep("mCT \\(-\\)", colnames(meanRawMeth))] -
  rowMeans2(meanRawMeth, cols = grep("mCT \\(-\\)", colnames(meanRawMeth)),
            na.rm = TRUE)
meanRawMeth_row_centered <-
  meanRawMeth_row_centered[!matrixStats::rowAnyNAs(meanRawMeth_row_centered), ]
meanRawMeth_row_centered_cp <- crossprod(meanRawMeth_row_centered)
meanRawMeth_row_centered_cp_svd <- svd(meanRawMeth_row_centered_cp)
meanRawMeth_row_centered_cp_pve <- meanRawMeth_row_centered_cp_svd$d /
  sum(meanRawMeth_row_centered_cp_svd$d)

pdf(paste0("../figures/PCA.", w, "kb_bins_mean_raw_mCH.pdf"))
plot(meanRawMeth_row_centered_cp_svd$u[, 1:2],
     col = dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                            grepl("BA9", cn) ~ "deepskyblue",
                            grepl("HC", cn) ~ "darkgrey",
                            grepl("NA", cn) ~ "chocolate1"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC1 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[1],
                         1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     main = paste0("raw mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"),
       legend = c("BA24", "BA9", "HC", "NAcc"),
       bty = "n")

plot(meanRawMeth_row_centered_cp_svd$u[, 1:2],
     col = dplyr::case_when(grepl("5248", cn) ~
                              "aquamarine",
                            grepl("5284", cn) ~
                              "chocolate1",
                            grepl("5552", cn) ~
                              "darkgrey",
                            grepl("5569", cn) ~
                              "deepskyblue",
                            grepl("5570", cn) ~
                              "firebrick3",
                            grepl("5628", cn) ~
                              "purple"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC1 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[1],
                         1),
                   "%)"),
     ylab = paste0("PC2 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     main = paste0("raw mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("aquamarine", "chocolate1", "darkgrey", "deepskyblue",
                "firebrick3", "purple"),
       legend = c("5248", "5284", "5552", "5569", "5570", "5628"),
       bty = "n")

plot(meanRawMeth_row_centered_cp_svd$u[, 2:3],
     col = dplyr::case_when(grepl("BA24", cn) ~ "deeppink",
                            grepl("BA9", cn) ~ "deepskyblue",
                            grepl("HC", cn) ~ "darkgrey",
                            grepl("NA", cn) ~ "chocolate1"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC2 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     ylab = paste0("PC3 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[3],
                         1),
                   "%)"),
     main = paste0("raw mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"),
       legend = c("BA24", "BA9", "HC", "NAcc"),
       bty = "n")

plot(meanRawMeth_row_centered_cp_svd$u[, 2:3],
     col = dplyr::case_when(grepl("5248", cn) ~
                              "aquamarine",
                            grepl("5284", cn) ~
                              "chocolate1",
                            grepl("5552", cn) ~
                              "darkgrey",
                            grepl("5569", cn) ~
                              "deepskyblue",
                            grepl("5570", cn) ~
                              "firebrick3",
                            grepl("5628", cn) ~
                              "purple"),
     pch = dplyr::case_when(grepl("mCA \\(\\+", cn) ~ "A",
                            grepl("mCA \\(\\-", cn) ~ "a",
                            grepl("mCT \\(\\+", cn) ~ "T",
                            grepl("mCT \\(\\-", cn) ~ "t"),
     xlab = paste0("PC2 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[2],
                         1),
                   "%)"),
     ylab = paste0("PC3 (",
                   round(100 * meanRawMeth_row_centered_cp_pve[3],
                         1),
                   "%)"),
     main = paste0("raw mCH (", w, " kb bins)"))
legend("bottomleft",
       pch = c("A", "a", "T", "t"),
       legend = c("mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       bty = "n")
legend("bottomright",
       fill = c("aquamarine", "chocolate1", "darkgrey", "deepskyblue",
                "firebrick3", "purple"),
       legend = c("5248", "5284", "5552", "5569", "5570", "5628"),
       bty = "n")
dev.off()

### ============================================================================
### Illustrative scatterplots
###

meanMeth <- do.call(cbind, list_of_meanMeth)
cn <- paste0(colnames(meanMeth), "_",
             rep(names(list_of_meanMeth),
                 times = sapply(list_of_meanMeth, ncol)))
colnames(meanMeth) <- cn
meanRawMeth <- do.call(cbind, list_of_meanRawMeth)
colnames(meanRawMeth) <- cn

png(paste0("../figures/mCH.scatterplots.", w, "kb.png"),
    width = 480,
    height = 480)
par(mfrow = c(2, 2))
repete::spwt(meanMeth[, "5248_BA24_pos_mCA (+)"],
             meanMeth[, "5248_BA24_pos_mCA (-)"],
             main = paste0("5248 BA24 (NeuN+)\nrho = ",
                           round(cor(meanMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanMeth[, "5248_BA24_pos_mCA (-)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("mCA (+) ", w, " kb bins"),
             ylab = paste0("mCA (-) ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanMeth[, "5248_BA24_pos_mCA (-)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanMeth[, "5248_BA24_pos_mCA (+)"],
             meanMeth[, "5248_BA24_pos_mCT (+)"],
             main = paste0("5248 BA24 (NeuN+)\nrho = ",
                           round(cor(meanMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanMeth[, "5248_BA24_pos_mCT (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("mCA (+) ", w, " kb bins"),
             ylab = paste0("mCT (+) ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanMeth[, "5248_BA24_pos_mCT (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanMeth[, "5248_BA24_pos_mCA (+)"],
             meanMeth[, "5248_NA_pos_mCA (+)"],
             main = paste0("5248 mCA (+) (NeuN+)\nrho = ",
                           round(cor(meanMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanMeth[, "5248_NA_pos_mCA (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("BA24 ", w, " kb bins"),
             ylab = paste0("NAcc ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanMeth[, "5248_NA_pos_mCA (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanMeth[, "5248_BA24_pos_mCA (+)"],
             meanMeth[, "5284_BA24_pos_mCA (+)"],
             main = paste0("BA24 mCA (+) (NeuN+)\nrho = ",
                           round(cor(meanMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanMeth[, "5284_BA24_pos_mCA (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("5248 ", w, " kb bins"),
             ylab = paste0("5284 ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanMeth[, "5284_BA24_pos_mCA (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

dev.off()

png(paste0("../figures/raw_mCH.scatterplots.", w, "kb.png"),
    width = 480,
    height = 480)
par(mfrow = c(2, 2))
repete::spwt(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
             meanRawMeth[, "5248_BA24_pos_mCA (-)"],
             main = paste0("raw 5248 BA24 (NeuN+)\nrho = ",
                           round(cor(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanRawMeth[, "5248_BA24_pos_mCA (-)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("mCA (+) ", w, " kb bins"),
             ylab = paste0("mCA (-) ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanRawMeth[, "5248_BA24_pos_mCA (-)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
             meanRawMeth[, "5248_BA24_pos_mCT (+)"],
             main = paste0("raw 5248 BA24 (NeuN+)\nrho = ",
                           round(cor(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanRawMeth[, "5248_BA24_pos_mCT (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("mCA (+) ", w, " kb bins"),
             ylab = paste0("mCT (+) ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanRawMeth[, "5248_BA24_pos_mCT (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
             meanRawMeth[, "5248_NA_pos_mCA (+)"],
             main = paste0("raw 5248 mCA (+) (NeuN+)\nrho = ",
                           round(cor(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanRawMeth[, "5248_NA_pos_mCA (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("BA24 ", w, " kb bins"),
             ylab = paste0("NAcc ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanRawMeth[, "5248_NA_pos_mCA (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

repete::spwt(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
             meanRawMeth[, "5284_BA24_pos_mCA (+)"],
             main = paste0("raw BA24 mCA (+) (NeuN+)\nrho = ",
                           round(cor(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanRawMeth[, "5284_BA24_pos_mCA (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("5248 ", w, " kb bins"),
             ylab = paste0("5284 ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanRawMeth[, "5284_BA24_pos_mCA (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))

dev.off()

png(paste0("../figures/smoothed_vs_raw_mCH.scatterplots.", w, ".kb.png"),
    width = 480,
    height = 480)
repete::spwt(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
             meanMeth[, "5248_BA24_pos_mCA (+)"],
             main = paste0("5248 BA24 mCA (+) (NeuN+)\nrho = ",
                           round(cor(meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                                     meanMeth[, "5248_BA24_pos_mCA (+)"],
                                     method = "spearman",
                                     use = "complete.obs"), 2)),
             xlab = paste0("raw mCA (+) ", w, " kb bins"),
             ylab = paste0("smoothed mCA (+) ", w, " kb bins"),
             xlim = c(0, 1),
             ylim = c(0, 1),
             col = scales::alpha("black", 0.2))
abline(a = 0, b = 1, lty = 2)
df <- data.frame(x = meanRawMeth[, "5248_BA24_pos_mCA (+)"],
                 y = meanMeth[, "5248_BA24_pos_mCA (+)"])
fit <- lm(y ~ x, data = df[complete.cases(df), ])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(y ~ 0 + x, data = df[complete.cases(df), ])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
legend("top",
       legend = c("y = x", "y ~ 0 + x", "y ~ x"),
       col = c("black", "red", "blue"),
       lty = c(2, 2, 2))
dev.off()

# TODO: Make mCG versions of these scatterplots (actually, could just/also show
#       correlation heatmaps)

### ============================================================================
### mCH boxplots
###

meanMeth <- do.call(cbind, list_of_meanMeth)
cn <- paste0(colnames(meanMeth), "_",
             rep(names(list_of_meanMeth),
                 times = sapply(list_of_meanMeth, ncol)))
colnames(meanMeth) <- cn

pdf(paste0("../figures/mCH.", w, "kb_bins.boxplots.pdf"))
par(mfrow = c(2, 2))
x <- meanMeth[, grep("mCA \\(\\+\\)", colnames(meanMeth))]
boxplot(x,
        outline = FALSE,
        col = dplyr::case_when(grepl("BA24", colnames(x)) ~
                                 "deeppink",
                               grepl("BA9", colnames(x)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(x)) ~
                                 "darkgrey",
                               grepl("NA", colnames(x)) ~
                                 "chocolate1"),
        xaxt = "n",
        main = paste0("mCA (+): ", w, " kb bins"),
        ylim = c(0, 0.25))
x <- meanMeth[, grep("mCA \\(-\\)", colnames(meanMeth))]
boxplot(x,
        outline = FALSE,
        col = dplyr::case_when(grepl("BA24", colnames(x)) ~
                                 "deeppink",
                               grepl("BA9", colnames(x)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(x)) ~
                                 "darkgrey",
                               grepl("NA", colnames(x)) ~
                                 "chocolate1"),
        xaxt = "n",
        main = paste0("mCA (-): ", w, " kb bins"),
        ylim = c(0, 0.25))
x <- meanMeth[, grep("mCT \\(\\+\\)", colnames(meanMeth))]
boxplot(x,
        outline = FALSE,
        col = dplyr::case_when(grepl("BA24", colnames(x)) ~
                                 "deeppink",
                               grepl("BA9", colnames(x)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(x)) ~
                                 "darkgrey",
                               grepl("NA", colnames(x)) ~
                                 "chocolate1"),
        xaxt = "n",
        main = paste0("mCT (+): ", w, " kb bins"),
        ylim = c(0, 0.25))
x <- meanMeth[, grep("mCT \\(-\\)", colnames(meanMeth))]
boxplot(x,
        outline = FALSE,
        col = dplyr::case_when(grepl("BA24", colnames(x)) ~
                                 "deeppink",
                               grepl("BA9", colnames(x)) ~
                                 "deepskyblue",
                               grepl("HC_pos", colnames(x)) ~
                                 "darkgrey",
                               grepl("NA", colnames(x)) ~
                                 "chocolate1"),
        xaxt = "n",
        main = paste0("mCT (-): ", w, " kb bins"),
        ylim = c(0, 0.25))

dev.off()
