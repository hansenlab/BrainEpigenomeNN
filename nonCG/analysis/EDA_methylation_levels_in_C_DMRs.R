# EDA of methylation levels in CG/CH DMRs and blocks
# Peter Hickey
# 2017-08-18

library(SummarizedExperiment)
library(repete)
library(viridis)
library(gplots)
library(Matrix)

# ------------------------------------------------------------------------------
# Pairs plots with correlations
#

list_of_SE <- readRDS("../objects/list_of_SE.DMRs_and_block.rds")

pdf("../figures/pairsWithCor.DMRs_and_blocks.pdf")
lapply(names(list_of_SE), function(n) {
  SE <- list_of_SE[[n]]
  cns <- colnames(SE)
  ans <- assayNames(SE)
  lapply(cns, function(cn) {
    x <- do.call(cbind, lapply(assays(SE), function(xx) xx[, cn]))
    repete::pairsWithCor(x,
                         main = paste0(n, ": ", cn),
                         xlim = c(0, 1),
                         ylim = c(0, 1))
  })
})
dev.off()

# ------------------------------------------------------------------------------
# Correlation heatmaps
#

list_of_cors <- lapply(list_of_SE, function(SE) {
  cns <- colnames(SE)
  ans <- assayNames(SE)
  val <- lapply(cns, function(cn) {
    x <- do.call(cbind, lapply(assays(SE), function(xx) xx[, cn]))
    cor(x, use = "complete.obs", method = "spearman")
  })
  names(val) <- cns
  val
})
list_of_cor_matrices <- lapply(list_of_cors, function(cors) {
  val <- do.call(bdiag, cors)
  rownames(val) <- sapply(cors, rownames)
  colnames(val) <- sapply(cors, colnames)
  val
})

hmcol <- viridis(100)
side_col <- c("BA24_pos" = "deeppink",
              "BA9_pos" = "deepskyblue",
              "HC_pos" = "darkgrey",
              "NA_pos" = "chocolate1")

pdf("../figures/correlation_heatmaps.DMRs_and_blocks.pdf")
lapply(names(list_of_cor_matrices), function(n) {
  cor_matrix <- list_of_cor_matrices[[n]]
  heatmap.2(as.matrix(cor_matrix),
            Rowv = NULL,
            Colv = NULL,
            trace = "none",
            col = hmcol,
            margins = c(4, 5),
            srtCol = 30,
            density.info = "none",
            main = n,
            ColSideColors = side_col[rep(names(list_of_cors[[n]]),
                                         each = ncol(list_of_cors[[n]][[1]]))],
            RowSideColors = side_col[rep(names(list_of_cors[[n]]),
                                         each = ncol(list_of_cors[[n]][[1]]))])
})
dev.off()
