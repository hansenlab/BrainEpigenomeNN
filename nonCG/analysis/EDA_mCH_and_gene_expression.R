# EDA of mCH and its association with gene expression
# 2018-02-02
# Peter Hickey

library(bsseq)
library(dplyr)
library(DelayedMatrixStats)
library(ggplot2)
library(limma)
library(scales)

### ============================================================================
### Load data
###

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$successful_permutations <= 0.05, ]
})
CH_DMRs <- unname(Reduce(c, list_of_candidate_CH_DMRs))
CH_DMRs$meanDiff <- CH_DMRs$`NA_pos` - CH_DMRs$`BA9_pos`

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

# ------------------------------------------------------------------------------
# mCG
#

CG_small_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.small.sorted.somatic.all"))
CG_small_BSseq <- CG_small_BSseq[, grep("pos", colnames(CG_small_BSseq))]
CG_small_BSseq$col <- CG_small_BSseq$Tissue_color

CG_large_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.large.sorted.somatic.all"))
CG_large_BSseq <- CG_large_BSseq[, grep("pos", colnames(CG_large_BSseq))]
CG_large_BSseq$col <- CG_large_BSseq$Tissue_color

list_of_CG_BSseq <- list("mCG (S)" = CG_small_BSseq,
                         "mCG (L)" = CG_large_BSseq)

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
                                 ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq$Tissue <- gsub("NA", "NAcc", BSseq$Tissue)
  BSseq[, grep("pos", colnames(BSseq))]
})


# ------------------------------------------------------------------------------
# mC
#

list_of_BSseq <- c(list_of_CH_BSseq, list_of_CG_BSseq)

### ============================================================================
### Update `rna_atac_meth` data frame with mCH data
###

# NOTE: When necessary, using mCA (+) to represent all mCH

# ------------------------------------------------------------------------------
# CH-DMR in promoter?
#

rna_atac_meth <- mutate(rna_atac_meth,
                          pCHDMR = overlapsAny(promoters_by_gene[gene],
                                               CH_DMRs,
                                               ignore.strand = TRUE))

# ------------------------------------------------------------------------------
# What is the max(|meanDiff|) in those pCH-DMRs
#

ol <- findOverlaps(promoters_by_gene[rna_atac_meth$gene], CH_DMRs,
                   ignore.strand = TRUE)
md <- tapply(X = CH_DMRs$meanDiff[subjectHits(ol)],
       INDEX = queryHits(ol),
       FUN = function(x) x[which.max(abs(x))])

rna_atac_meth$pCHDMRMeanDiff <- NA_real_
rna_atac_meth$pCHDMRMeanDiff[as.integer(names(md))] <- md

# ------------------------------------------------------------------------------
# CH-DMR in gene body?
#

rna_atac_meth <- mutate(rna_atac_meth,
                        gbCHDMR = overlapsAny(gencode_features$genes[gene],
                                              CH_DMRs,
                                              ignore.strand = TRUE))

# ------------------------------------------------------------------------------
# What is the max(|meanDiff|) in those gbCH-DMRs
#

ol <- findOverlaps(gencode_features$genes[rna_atac_meth$gene], CH_DMRs,
                   ignore.strand = TRUE)
md <- tapply(X = CH_DMRs$meanDiff[subjectHits(ol)],
             INDEX = queryHits(ol),
             FUN = function(x) x[which.max(abs(x))])

rna_atac_meth$gbCHDMRMeanDiff <- NA_real_
rna_atac_meth$gbCHDMRMeanDiff[as.integer(names(md))] <- md


# ------------------------------------------------------------------------------
# What percentage of gene body is covered by a CH-DMR
#


ol <- findOverlaps(gencode_features$genes[rna_atac_meth$gene], reduce(CH_DMRs))
w <- width(pintersect(gencode_features$genes[rna_atac_meth$gene][queryHits(ol)],
                      reduce(CH_DMRs)[subjectHits(ol)],
                      ignore.strand = TRUE))
tw <- tapply(X = w,
             INDEX = queryHits(ol),
             FUN = sum)
rna_atac_meth$percCHDMR <- 0
rna_atac_meth$percCHDMR[as.integer(names(tw))] <- 100 *
  tw / width(gencode_features$genes[rna_atac_meth$gene[as.integer(names(tw))]])

# ------------------------------------------------------------------------------
# mean(mCH) over gene body
#

m <- getMeth(BSseq = list_of_CH_BSseq[["mCA (+)"]],
             regions = unstrand(gencode_features$genes[rna_atac_meth$gene]),
             type = "smooth",
             what = "perRegion")
ave_m <- sapply(c("BA9", "BA24", "NA", "HC"), function(tissue) {
  rowMeans2(m, cols = grep(tissue, colnames(m)))
})
colnames(ave_m) <- paste0("mCA_pos_", colnames(ave_m))
rna_atac_meth <- bind_cols(rna_atac_meth, as_data_frame(ave_m))
rna_atac_meth <- mutate(rna_atac_meth,
                        mCHMeanDiff = mCA_pos_NA - mCA_pos_BA9)

# ------------------------------------------------------------------------------
# Save data fram
#

saveRDS(rna_atac_meth, "../tmp/rna_atac_meth.rds")

### ============================================================================
### Exploratory plots and summaries
###

# ------------------------------------------------------------------------------
# logFC vs. gene body meanDiff(mCH) stratified by DE-status
#

x <- select(rna_atac_meth, gene, DE, expLogFC, mCHMeanDiff)

group_by(x, DE) %>%
  summarise(cor = cor(mCHMeanDiff, y = expLogFC, use = "complete.obs"))

xlim <- range(x$mCHMeanDiff, na.rm = TRUE)
ylim <- range(x$expLogFC, na.rm = TRUE)
par(mfrow = c(1, 3))
plot(expLogFC ~ mCHMeanDiff,
     data = x,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(x$mCHMeanDiff, x$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x$expLogFC ~ 0 + x$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
xDE <- filter(x, DE)
plot(expLogFC ~ mCHMeanDiff,
     data = xDE,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(xDE$mCHMeanDiff, xDE$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(xDE$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(xDE))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(xDE$expLogFC ~ 0 + xDE$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)
xNDE <- filter(x, !DE)
plot(expLogFC ~ mCHMeanDiff,
     data = xNDE,
     main = paste0("NAcc_pos vs BA9_pos: cor = ",
                   round(cor(xNDE$mCHMeanDiff, xNDE$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(xNDE$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(xNDE))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(xNDE$expLogFC ~ 0 + xNDE$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# ------------------------------------------------------------------------------
# logFC vs. gene body meanDiff(mCH) stratified by DE-status and pCH-DMR status
#

ggplot(aes(x = mCHMeanDiff, y = expLogFC), data = rna_atac_meth) +
  geom_point() +
  facet_grid(DE ~ pCHDMR, labeller = label_both)

fisher.test(table(rna_atac_meth$DE, rna_atac_meth$pCHDMR))

group_by(rna_atac_meth, DE, pCHDMR) %>%
  summarise(cor = cor(mCHMeanDiff, expLogFC, use = "complete.obs"))

x <- select(rna_atac_meth, gene, DE, pCHDMR, expLogFC, mCHMeanDiff)

xlim <- range(x$mCHMeanDiff, na.rm = TRUE)
ylim <- range(x$expLogFC, na.rm = TRUE)
par(mfrow = c(2, 2))

x1 <- x[!x$DE & !x$pCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x1,
     main = paste0("!DE & !pCH-DMR: cor = ",
                   round(cor(x1$mCHMeanDiff, x1$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x1$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x1))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x1$expLogFC ~ 0 + x1$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x2 <- x[!x$DE & x$pCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x2,
     main = paste0("!DE & pCH-DMR: cor = ",
                   round(cor(x2$mCHMeanDiff, x2$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x2$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x2))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x2$expLogFC ~ 0 + x2$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x3 <- x[x$DE & !x$pCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x3,
     main = paste0("DE & !pCH-DMR: cor = ",
                   round(cor(x3$mCHMeanDiff, x3$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x3$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x3))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x3$expLogFC ~ 0 + x3$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x4 <- x[x$DE & x$pCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x4,
     main = paste0("DE & pCH-DMR: cor = ",
                   round(cor(x4$mCHMeanDiff, x4$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x4$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x4))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x4$expLogFC ~ 0 + x4$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# ------------------------------------------------------------------------------
# Percent covered by CH-DMR stratified by DE-status and pCH-DMR status
#

ggplot(rna_atac_meth, aes(x = percCHDMR)) +
  geom_density(aes(col = interaction(DE, pCHDMR)))

ggplot(filter(rna_atac_meth, percCHDMR > 0), aes(x = percCHDMR)) +
  geom_density(aes(col = interaction(DE, pCHDMR)))

# ------------------------------------------------------------------------------
# logFC vs. gene body meanDiff(mCH) stratified by DE-status and gbCH-DMR status
#

ggplot(aes(x = mCHMeanDiff, y = expLogFC), data = rna_atac_meth) +
  geom_point() +
  facet_grid(DE ~ gbCHDMR, labeller = label_both)

fisher.test(table(rna_atac_meth$DE, rna_atac_meth$gbCHDMR))

group_by(rna_atac_meth, DE, gbCHDMR) %>%
  summarise(cor = cor(mCHMeanDiff, expLogFC, use = "complete.obs"))

x <- select(rna_atac_meth, gene, DE, gbCHDMR, expLogFC, mCHMeanDiff)

xlim <- range(x$mCHMeanDiff, na.rm = TRUE)
ylim <- range(x$expLogFC, na.rm = TRUE)
par(mfrow = c(2, 2))

x1 <- x[!x$DE & !x$gbCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x1,
     main = paste0("!DE & !gbCH-DMR: cor = ",
                   round(cor(x1$mCHMeanDiff, x1$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x1$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x1))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x1$expLogFC ~ 0 + x1$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x2 <- x[!x$DE & x$gbCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x2,
     main = paste0("!DE & gbCH-DMR: cor = ",
                   round(cor(x2$mCHMeanDiff, x2$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x2$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x2))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x2$expLogFC ~ 0 + x2$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x3 <- x[x$DE & !x$gbCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x3,
     main = paste0("DE & !gbCH-DMR: cor = ",
                   round(cor(x3$mCHMeanDiff, x3$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x3$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x3))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x3$expLogFC ~ 0 + x3$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

x4 <- x[x$DE & x$gbCHDMR, ]
plot(expLogFC ~ mCHMeanDiff,
     data = x4,
     main = paste0("DE & gbCH-DMR: cor = ",
                   round(cor(x4$mCHMeanDiff, x4$expLogFC,
                             use = "complete.obs"), 2)),
     xlab = "Gene body meanDiff(mCA (+))",
     ylab = "log-FC Gene Expression",
     pch = 16,
     cex = 0.8,
     col = ifelse(x4$DE, "orange", "gray30"),
     sub = paste0("n = ", sum(complete.cases(x4))),
     xlim = xlim,
     ylim = ylim)
abline(h = 0, v = 0, col = "gray10", lty = 2, lwd = 2)
u <- lm(x4$expLogFC ~ 0 + x4$mCHMeanDiff)
coef(summary(u))
abline(u, col = "red", lwd = 2)

# ------------------------------------------------------------------------------
# MA-plots of gene body methylation
#

x <- as.matrix(
  select(rna_atac_meth, mCA_pos_BA24, mCA_pos_BA9, mCA_pos_HC, mCA_pos_NA))

xlim <- c(0, 0.6)
ylim <- c(-0.5, 0.5)
par(mfrow = c(2, 3))
lapply(list(c(1, 2), c(1, 3), c(2, 3), c(1, 4), c(2, 4), c(3, 4)),
       function(columns) {
         mdplot(x,
                columns = columns,
                status = ifelse(rna_atac_meth$DE, "DE", "Non-DE"),
                hl.col = alpha("orange", 0.3),
                bg.col = alpha("black", 0.3),
                xlim = xlim,
                ylim = ylim)
       })

# Same plot, but stratified be DE rather than coloured
par(mfcol = c(4, 3))
lapply(list(c(1, 2), c(1, 3), c(2, 3), c(1, 4), c(2, 4), c(3, 4)),
       function(columns) {
         mdplot(x[rna_atac_meth$DE, ],
                columns = columns,
                bg.col = alpha("orange", 0.8),
                xlim = xlim,
                ylim = ylim)
         mdplot(x[!rna_atac_meth$DE, ],
                columns = columns,
                bg.col = alpha("black", 0.3),
                xlim = xlim,
                ylim = ylim)
       })

### ----------------------------------------------------------------------------
### Does it matter where in the gene body the CH-DMR (more generally, DiffEpi) is?
###

features <- list(`CG-DMRs` = dmrs_NAvsBA9pos,
                 `CH-DMRs` = reduce(CH_DMRs),
                 `DARs` = dars_pos,
                 `bigDARs` = dars_pos[abs(dars_pos$logFC) > 1])

hitsAlongGeneBody <- function(feature, g, ss = 0.01) {
  stopifnot(ss > 0, ss < 1)

  p <- seq(ss, 1, ss)
  w <- ss * width(g)
  ol_tbl <- sapply(p, function(pp) {
    on_plus <- which(strand(g) == "+" | strand(g) == "*")
    on_plus_TSS <- start(g)[on_plus]
    end(g)[on_plus] <- on_plus_TSS + pp * width(g)[on_plus] - 1L
    on_minus <- which(strand(g) == "-")
    on_minus_TSS <- end(g)[on_minus]
    start(g)[on_minus] <- on_minus_TSS - pp * width(g)[on_minus] + 1L
    g <- resize(g, w, fix = "end")
    overlapsAny(g, feature)
  })
  gene_body_perc <- 100 * colMeans2(ol_tbl)
  return(list(x = p,
              y = gene_body_perc))
}

# ------------------------------------------------------------------------------
# Gene body
#

# All genes
g <- gencode_features$genes[rna_atac_meth$gene]
vals <- lapply(features,
               hitsAlongGeneBody,
               g = g,
               ss = 0.01)

pdf("../figures/location_of_DiffEpi_along_gene_body.all_genes.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- rep(mean(val$y), length(x))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()

# All genes with CGI promoters
g <- gencode_features$genes[intersect(rna_atac_meth$gene,
                                      unlist(cgi_promoters$GENEID))]
vals <- lapply(features,
               hitsAlongGeneBody,
               g = g,
               ss = 0.01)

pdf("../figures/location_of_DiffEpi_along_gene_body.genes_with_CGI_promoter.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- rep(mean(val$y), length(x))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()

# All genes without CGI promoters
g <- gencode_features$genes[setdiff(rna_atac_meth$gene,
                                    unlist(cgi_promoters$GENEID))]
vals <- lapply(features,
               hitsAlongGeneBody,
               g = g,
               ss = 0.01)

pdf("../figures/location_of_DiffEpi_along_gene_body.genes_without_CGI_promoter.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- rep(mean(val$y), length(x))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()

# ------------------------------------------------------------------------------
# Gene body with (scaled) upstream and downstream
#

# Number of gene body lengths to go upstream/downstream
ngb <- 2

# All genes
g <- gencode_features$genes[rna_atac_meth$gene]
w <- width(g)
vals <- lapply(features,
               hitsAlongGeneBody,
               g = resize(resize(g, w + ngb * w, fix = "start"),
                          w + 2 * ngb * w, fix = "end"),
               ss = 0.01 / (2 * ngb + 1))
pdf("../figures/location_of_DiffEpi_along_gene_body.all_genes.upstream_and_downstream.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  inc <- 2 * ngb + 1
  first_cut <- length(x) * ngb / inc
  second_cut <- length(x) * (ngb + 1) / inc
  x[seq(1, first_cut)] <- -rev(seq(0, ngb * 100 - 1, length.out = first_cut))
  x[seq(first_cut, second_cut)] <- seq(0, 100, length.out = second_cut - first_cut + 1)
  x[seq(second_cut, length(x))] <- seq(100, (ngb + 1) * 100, length.out = length(x) - second_cut + 1)
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- do.call(c, tapply(val$y,
                       c(rep(1, first_cut),
                         rep(2, second_cut - first_cut + 1),
                         rep(3, length(x) - second_cut - 1)),
                       function(x) rep(mean(x), length(x)),
                       simplify = FALSE))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()

# All genes with CGI promoters
g <- gencode_features$genes[intersect(rna_atac_meth$gene,
                                      unlist(cgi_promoters$GENEID))]
w <- width(g)
vals <- lapply(features,
               hitsAlongGeneBody,
               g = resize(resize(g, w + ngb * w, fix = "start"),
                          w + 2 * ngb * w, fix = "end"),
               ss = 0.01 / (2 * ngb + 1))

pdf("../figures/location_of_DiffEpi_along_genes.genes_with_CGI_promoter.upstream_and_downstream.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  inc <- 2 * ngb + 1
  first_cut <- length(x) * ngb / inc
  second_cut <- length(x) * (ngb + 1) / inc
  x[seq(1, first_cut)] <- -rev(seq(0, ngb * 100 - 1, length.out = first_cut))
  x[seq(first_cut, second_cut)] <- seq(0, 100, length.out = second_cut - first_cut + 1)
  x[seq(second_cut, length(x))] <- seq(100, (ngb + 1) * 100, length.out = length(x) - second_cut + 1)
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- do.call(c, tapply(val$y,
                         c(rep(1, first_cut),
                           rep(2, second_cut - first_cut + 1),
                           rep(3, length(x) - second_cut - 1)),
                         function(x) rep(mean(x), length(x)),
                         simplify = FALSE))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()

# All genes without CGI promoters
g <- gencode_features$genes[setdiff(rna_atac_meth$gene,
                                    unlist(cgi_promoters$GENEID))]
w <- width(g)
vals <- lapply(features,
               hitsAlongGeneBody,
               g = resize(resize(g, w + ngb * w, fix = "start"),
                          w + 2 * ngb * w, fix = "end"),
               ss = 0.01 / (2 * ngb + 1))

pdf("../figures/location_of_DiffEpi_along_genes.genes_without_CGI_promoter.upstream_and_downstream.pdf",
    height = 7,
    width = 7)
par(mfrow = c(2, 2))
lapply(seq_along(vals), function(k) {
  val <- vals[[k]]
  x <- 100 * val$x
  inc <- 2 * ngb + 1
  first_cut <- length(x) * ngb / inc
  second_cut <- length(x) * (ngb + 1) / inc
  x[seq(1, first_cut)] <- -rev(seq(0, ngb * 100 - 1, length.out = first_cut))
  x[seq(first_cut, second_cut)] <- seq(0, 100, length.out = second_cut - first_cut + 1)
  x[seq(second_cut, length(x))] <- seq(100, (ngb + 1) * 100, length.out = length(x) - second_cut + 1)
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = paste0("% genes with ", names(vals)[k]),
       main = names(vals)[k],
       type = "s",
       ylim = c(0, 1.2 * max(val$y)))
  m <- do.call(c, tapply(val$y,
                         c(rep(1, first_cut),
                           rep(2, second_cut - first_cut + 1),
                           rep(3, length(x) - second_cut - 1)),
                         function(x) rep(mean(x), length(x)),
                         simplify = FALSE))
  lines(x = x, y = m, col = "red", lwd = 2)
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
})
dev.off()
