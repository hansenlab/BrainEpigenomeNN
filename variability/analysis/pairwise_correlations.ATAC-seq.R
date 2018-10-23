# Look at pairwise correlation of ATAC-seq between samples
# Peter Hickey
# 2017-08-10

### ============================================================================
### Setup
###

library(SummarizedExperiment)
library(dplyr)
library(gtools)
library(ggplot2)
library(tidyr)
library(purrr)
library(coop)
library(gplots)
library(RColorBrewer)
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

### ============================================================================
### Compute correlations
###

# ------------------------------------------------------------------------------
# Using all bins tested for differential accessibility
#

pcors <- pcor(cpm)
x <- pcors
x[upper.tri(x)] <- NA_real_
pcors_df <- data_frame(cor = as.vector(x),
                       S1 = rep(colnames(ATAC_SE), times = ncol(ATAC_SE)),
                       S2 = rep(colnames(ATAC_SE), each = ncol(ATAC_SE))) %>%
  # NOTE: Remove diagonal and upper triangle
  filter(S1 != S2,
         !is.na(cor)) %>%
  mutate(S1_Donor = colData(ATAC_SE)[S1, "DONOR"],
         S2_Donor = colData(ATAC_SE)[S2, "DONOR"],
         S1_Tissue = colData(ATAC_SE)[S1, "TISSUE"],
         S2_Tissue = colData(ATAC_SE)[S2, "TISSUE"],
         S1_NeuN = colData(ATAC_SE)[S1, "NEUN"],
         S2_NeuN = colData(ATAC_SE)[S2, "NEUN"])

# ------------------------------------------------------------------------------
# Save objects
#

ATAC_cor_matrix <- pcors
dimnames(ATAC_cor_matrix) <- list(colnames(ATAC_SE), colnames(ATAC_SE))
ATAC_cor_df <- pcors_df
ATAC_colData <- colData(ATAC_SE)

save(ATAC_cor_matrix, ATAC_cor_df, ATAC_colData,
     file = "../objects/ATAC-seq_pairwise_correlations.RData")

### ============================================================================
### Plot correlations
###

# ------------------------------------------------------------------------------
# Heatmap of correlation matrix
#

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf("../figures/ATAC-seq_union_narrowPeak_reduced.correlation_heatmap.pdf")
heatmap.2(pcors,
          trace = "none",
          ColSideColors = colData(ATAC_SE)[["TISSUE_COLOR"]],
          RowSideColors = ifelse(colData(ATAC_SE)[["NEUN"]] == "pos",
                                 "darkgreen", "purple"),
          col = hmcol,
          cexRow = -0.2 + 1 / log10(ncol(pcors)),
          cexCol = -0.2 + 1 / log10(ncol(pcors)),
          main = "ATAC-seq (cpm) overall peaks")
dev.off()
