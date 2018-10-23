# Scatterplots of mCH comparing strand and context
# Peter Hickey
# 2017-12-19

library(bsseq)

### ============================================================================
### Load data
###

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "neg_CA", "pos_CT", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
                                 ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq[, BSseq$NeuN == "pos"]
})
# NOTE: Fix up colours
list_of_CH_BSseq <- lapply(list_of_CH_BSseq, function(BSseq) {
  BSseq$col <- BSseq$Tissue_color
  BSseq
})

# ------------------------------------------------------------------------------
# Binned data
#

list_of_binned_mCH <- readRDS("../objects/list_of_meanMeth.1kb_bins.rds")
bins <- readRDS("../objects/bins.1kb.rds")

### ============================================================================
### Plots
###

# Keep only those bins with at least min_p proportion of Cs with at least
# min_cov coverage. E.g., 70% of CTs in bin with coverage >= 1
keeper <- function(cov1, cov2, C1, C2, min_p = 0.7, min_cov = 1) {
  if (is.null(cov1) || is.null(cov2)) {
    return(FALSE)
  }
  (sum(cov1 >= min_cov) >= (min_p * C1)) &&
    (sum(cov2 >= min_cov) >= (min_p * C2))
}

# ------------------------------------------------------------------------------
# Scatterplot of biological replicates
#

mCA_pos_5248_BA9_pos <- list_of_binned_mCH[["mCA (+)"]][, "5248_BA9_pos"]
mCA_pos_5569_BA9_pos <- list_of_binned_mCH[["mCA (+)"]][, "5569_BA9_pos"]
cov_CA_pos_5248_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCA (+)"]][, "5248_BA9_pos"], bins)
cov_CA_pos_5569_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCA (+)"]][, "5569_BA9_pos"], bins)

keep <- mapply(keeper,
               cov1 = cov_CA_pos_5248_BA9_pos,
               cov2 = cov_CA_pos_5569_BA9_pos,
               C1 = bins$CA,
               C2 = bins$CA)
sum(keep) / length(keep) # 89%

pdf("../figures/mCH_scatterplot.biological_replicates.pdf",
    width = 7,
    height = 7)
par(mfrow = c(1, 1))
smoothScatter(x = mCA_pos_5248_BA9_pos[keep],
              y = mCA_pos_5569_BA9_pos[keep],
              main = round(cor(mCA_pos_5248_BA9_pos[keep],
                               mCA_pos_5569_BA9_pos[keep],
                               use = "complete.obs"),
                           2),
              xlab = "5248_BA9_pos: mCA (+)",
              ylab = "5569_BA9_pos: mCA (+)",
              xlim = c(0, 0.4),
              ylim = c(0, 0.4),
              nrpoints = 1000)
abline(a = 0, b = 1, col = "red")
fit <- lm(mCA_pos_5569_BA9_pos[keep] ~ mCA_pos_5248_BA9_pos[keep])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(mCA_pos_5569_BA9_pos[keep] ~ 0 + mCA_pos_5248_BA9_pos[keep])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
dev.off()

# ------------------------------------------------------------------------------
# Scatterplot of strand effect
#

# TODO: A version with geom_density_2d() using all samples
mCA_pos <- list_of_binned_mCH[["mCA (+)"]]
mCA_neg <- list_of_binned_mCH[["mCA (-)"]]
cov_CA_pos <- getCoverage(list_of_CH_BSseq[["mCA (+)"]], bins)
cov_CA_neg <- getCoverage(list_of_CH_BSseq[["mCA (-)"]], bins)



mCA_pos_5248_BA9_pos <- list_of_binned_mCH[["mCA (+)"]][, "5248_BA9_pos"]
mCA_neg_5248_BA9_pos <- list_of_binned_mCH[["mCA (-)"]][, "5248_BA9_pos"]
cov_CA_pos_5248_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCA (+)"]][, "5248_BA9_pos"], bins)
cov_CA_neg_5248_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCA (-)"]][, "5248_BA9_pos"], bins)

keep <- mapply(keeper,
               cov1 = cov_CA_pos_5248_BA9_pos,
               cov2 = cov_CA_neg_5248_BA9_pos,
               C1 = bins$CA,
               C2 = bins$TG)
sum(keep) / length(keep) # 80%

pdf("../figures/mCH_scatterplot.strand_effect.pdf",
    width = 7,
    height = 7)
par(mfrow = c(1, 1))
smoothScatter(x = mCA_pos_5248_BA9_pos[keep],
              y = mCA_neg_5248_BA9_pos[keep],
              main = round(cor(mCA_pos_5248_BA9_pos[keep],
                               mCA_neg_5248_BA9_pos[keep],
                               use = "complete.obs"),
                           2),
              xlab = "5248_BA9_pos: mCA (+)",
              ylab = "5248_BA9_pos: mCA (-)",
              xlim = c(0, 0.4),
              ylim = c(0, 0.4),
              nrpoints = 1000)
abline(a = 0, b = 1, col = "red")
fit <- lm(mCA_neg_5248_BA9_pos[keep] ~ mCA_pos_5248_BA9_pos[keep])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(mCA_neg_5248_BA9_pos[keep] ~ 0 + mCA_pos_5248_BA9_pos[keep])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
dev.off()

# ------------------------------------------------------------------------------
# Scatterplot of context effect
#

mCA_pos_5248_BA9_pos <- list_of_binned_mCH[["mCA (+)"]][, "5248_BA9_pos"]
mCT_pos_5248_BA9_pos <- list_of_binned_mCH[["mCT (+)"]][, "5248_BA9_pos"]
cov_CA_pos_5248_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCA (+)"]][, "5248_BA9_pos"], bins)
cov_CT_pos_5248_BA9_pos <-
  getCoverage(list_of_CH_BSseq[["mCT (+)"]][, "5248_BA9_pos"], bins)

keep <- mapply(keeper,
               cov1 = cov_CA_pos_5248_BA9_pos,
               cov2 = cov_CT_pos_5248_BA9_pos,
               C1 = bins$CA,
               C2 = bins$CT)
sum(keep) / length(keep) # 84%

pdf("../figures/mCH_scatterplot.context_effect.pdf",
    width = 7,
    height = 7)
par(mfrow = c(1, 1))
smoothScatter(x = mCA_pos_5248_BA9_pos[keep],
              y = mCT_pos_5248_BA9_pos[keep],
              main = round(cor(mCA_pos_5248_BA9_pos[keep],
                               mCT_pos_5248_BA9_pos[keep],
                               use = "complete.obs"),
                           2),
              xlab = "5248_BA9_pos: mCA (+)",
              ylab = "5248_BA9_pos: mCT (+)",
              xlim = c(0, 0.4),
              ylim = c(0, 0.4),
              nrpoints = 1000)
abline(a = 0, b = 1, col = "red")
fit <- lm(mCT_pos_5248_BA9_pos[keep] ~ mCA_pos_5248_BA9_pos[keep])
abline(a = coef(fit)[1], b = coef(fit)[2], lty = 2, col = "blue")
fit0 <- lm(mCT_pos_5248_BA9_pos[keep] ~ 0 + mCA_pos_5248_BA9_pos[keep])
abline(a = 0, b = coef(fit0), lty = 2, col = "red")
dev.off()
