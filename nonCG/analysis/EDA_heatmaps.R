# EDA of mCG, mCH, and ATAC-seq data around CG-DMRs and CH-DMRs
# Peter Hickey
# 2017-09-15

library(bsseq)
library(EnrichedHeatmap)
library(matrixStats)
library(viridis)
library(circlize)
library(scales)

RefSeq_exons3 <- readRDS("../objects/RefSeq_exons3.rds")
source("EDA_functions.R")

extdir <- "../extdata"

options("mc.cores" = 18)

### ============================================================================
### Load data
###

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

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../../FlowSortingProject/Objects/All_BLOCK_POS_DMRs_fwer50.rda")
CG_DMRs <- dmrs_pos
CG_blocks <- makeGRangesFromDataFrame(sig_block_dmrs)

# ------------------------------------------------------------------------------
# mCH
#

strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")
CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
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

list_of_candidate_CH_DMRs <- readRDS("../objects/list_of_candidate_CH_DMRs.rds")
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$n <= 0.05]
})

### ============================================================================
### Average data in NeuN+ samples
###

# ------------------------------------------------------------------------------
# Regions
#

regions <- c("CG DMRs" = CG_DMRs, list_of_CH_DMRs)
regions <- lapply(regions, function(x) {
  resize(x, width(x) + 15000 * 2, fix = "center")
})
reduced_regions <- sort(reduce(unlist(GRangesList(lapply(regions, granges)))))

# ------------------------------------------------------------------------------
# GRanges of mean mCG and mCH
#

list_of_mean_mC_GR <-
  mclapply(c(list_of_CG_BSseq, list_of_CH_BSseq), function(BSseq) {
    message(nrow(BSseq))
    j <- grep("pos", colnames(BSseq))
    se <- SummarizedExperiment(
      assays = as.matrix(getMeth(BSseq[, j])),
      rowRanges = rowRanges(BSseq),
      colData = colData(BSseq[, j]))
    se <- subsetByOverlaps(se, reduced_regions)
    meth <- sapply(c("BA9", "BA24", "HC", "NA"), function(j) {
      rowMeans2(assay(se, withDimnames = FALSE), cols = grep(j, colnames(se)),
                na.rm = TRUE)
    })
    GR <- granges(se)
    mcols(GR) <- as(meth, "DataFrame")
    GR
  }, mc.cores = options("mc.cores"))
saveRDS(list_of_mean_mC_GR, "../objects/list_of_mean_mC_GR.rds")

### ============================================================================
### Colours
###

ATAC_col <- viridis(200)
mCG_col <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
mCH_col <- plasma(200)

### ============================================================================
### CG-DMRs
###

CG_DMRs_target <- resize(CG_DMRs,
                         width = 1,
                         fix = "center")

list_of_CG_DMRs_mat <- mclapply(list_of_mean_mC_GR, function(mC_GR) {
  cns <- colnames(mcols(mC_GR))
  setNames(mclapply(cns, function(cn) {
    normalizeToMatrix(signal = mC_GR,
                      target = CG_DMRs_target,
                      extend = 1500,
                      w = 30,
                      value_column = cn,
                      mean_mode = "absolute",
                      empty_value = NA,
                      smooth = TRUE)
  }, mc.cores = 6),
  cns)
}, mc.cores = 3)
saveRDS(list_of_CG_DMRs_mat, "../objects/list_of_CG_DMRs_mat.rds")

pdf(file = "../figures/CG-DMR_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCG (S)"]][["NA"]],
                col = mCG_col,
                name = "mCG (S)",
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NAcc_pos",
                width = 1,
                axis_name = c("-1500", "CG-DMR", "1500"),
                use_raster = TRUE,
                raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCG (S)"]][["BA9"]],
                  col = mCG_col,
                  name = "mCG (S)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "BA9_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCA (+)"]][["NA"]],
                  col = mCH_col,
                  name = "mCA (+)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "NAcc_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCA (+)"]][["BA9"]],
                  col = mCH_col,
                  name = "mCA (+)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "BA9_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCA (-)"]][["NA"]],
                  col = mCH_col,
                  name = "mCA (-)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "NAcc_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCA (-)"]][["BA9"]],
                  col = mCH_col,
                  name = "mCA (-)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "BA9_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCT (+)"]][["NA"]],
                  col = mCH_col,
                  name = "mCT (+)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "NAcc_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCT (+)"]][["BA9"]],
                  col = mCH_col,
                  name = "mCT (+)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "BA9_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCT (-)"]][["NA"]],
                  col = mCH_col,
                  name = "mCT (-)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "NAcc_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG") +
  EnrichedHeatmap(mat = list_of_CG_DMRs_mat[["mCT (-)"]][["BA9"]],
                  col = mCH_col,
                  name = "mCT (-)",
                  pos_line = FALSE,
                  column_title_gp = gpar(fontsize = 15),
                  axis_name_gp = gpar(fontsize = 6),
                  column_title = "BA9_pos",
                  width = 1,
                  axis_name = c("-1500", "CG-DMR", "1500"),
                  use_raster = TRUE,
                  raster_device = "CairoPNG")
dev.off()

### ============================================================================
### CH-DMRs
###

# TODO: Save `split` variables for CH-DMRs (CG-DMR overlap and k-means with k = 6)

list_of_CH_DMRs_mat <- lapply(list_of_CH_DMRs, function(CH_DMRs) {
  target <- resize(CH_DMRs,
                   width = 1,
                   fix = "center")
  mclapply(list_of_mean_mC_GR, function(mC_GR) {
    cns <- colnames(mcols(mC_GR))
    setNames(mclapply(cns, function(cn) {
      normalizeToMatrix(signal = mC_GR,
                        target = target,
                        extend = 15000,
                        w = 5,
                        value_column = cn,
                        mean_mode = "absolute",
                        empty_value = NA,
                        smooth = TRUE)
    }, mc.cores = 6),
    cns)
  }, mc.cores = 3)
})
saveRDS(list_of_CH_DMRs_mat, "../objects/list_of_CH_DMRs_mat.rds")

list_of_CH_DMRs_centered_mat <- lapply(list_of_CH_DMRs, function(CH_DMRs) {
  target <- CH_DMRs
  mclapply(list_of_mean_mC_GR, function(mC_GR) {
    cns <- colnames(mcols(mC_GR))
    setNames(mclapply(cns, function(cn) {
      normalizeToMatrix(signal = mC_GR,
                        target = target,
                        extend = 15000,
                        w = 30,
                        value_column = cn,
                        mean_mode = "absolute",
                        empty_value = NA,
                        smooth = TRUE)
    }, mc.cores = 6),
    cns)
  }, mc.cores = 3)
})
saveRDS(list_of_CH_DMRs_centered_mat,
        "../objects/list_of_CH_DMRs_centered_mat.rds")

# --------------------------------------------------------------------------------
# pos CA DMRs
#

pdf(file = "../figures/pos_CA-DMR_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCA (+)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  split = split,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  split = split,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/pos_CA-DMRs.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = resize(list_of_CH_DMRs[["mCA (+)"]][split == i][1:20], 1,
                     fix = "center"),
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

pdf(file = "../figures/pos_CA-DMR_centered_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCA (+)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  split = split,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["NA"]],
  col = mCH_col,
  name = "mCA (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  split = split,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    split = split,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/pos_CA-DMRs_centered.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = list_of_CH_DMRs[["mCA (+)"]][split == i][1:20],
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

# --------------------------------------------------------------------------------
# neg CA DMRs
#

pdf(file = "../figures/neg_CA-DMR_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCA (-)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/neg_CA-DMRs.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = resize(list_of_CH_DMRs[["mCA (-)"]][split == i][1:20], 1,
                     fix = "center"),
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

pdf(file = "../figures/neg_CA-DMR_centered_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCA (-)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["NA"]],
  col = mCH_col,
  name = "mCA (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCA (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/neg_CA-DMRs_centered.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = list_of_CH_DMRs[["mCA (-)"]][split == i][1:20],
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

# --------------------------------------------------------------------------------
# pos CT DMRs
#

pdf(file = "../figures/pos_CT-DMR_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCT (+)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf(file = "../figures/pos_CT-DMR_centered_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCT (+)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["NA"]],
  col = mCH_col,
  name = "mCT (+)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["NA"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (+)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/pos_CT-DMRs_centered.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = list_of_CH_DMRs[["mCT (+)"]][split == i][1:20],
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

# --------------------------------------------------------------------------------
# neg CT DMRs
#

pdf(file = "../figures/neg_CT-DMR_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCT (-)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "CH-DMR", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "CH-DMR", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/neg_CT-DMRs.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = resize(list_of_CH_DMRs[["mCT (-)"]][split == i][1:20], 1,
                     fix = "center"),
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

pdf(file = "../figures/neg_CT-DMR_centered_heatmap.mC.pdf",
    width = 14,
    height = 14)
EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    use_raster = TRUE,
    raster_device = "CairoPNG")

split <- ifelse(overlapsAny(list_of_CH_DMRs[["mCT (-)"]], CG_DMRs),
                "CG-DMR", "Not")

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")

set.seed(666)
split <- kmeans(list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
                centers = 4)$cluster

EnrichedHeatmap(
  mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["NA"]],
  col = mCH_col,
  name = "mCT (-)",
  pos_line = FALSE,
  column_title_gp = gpar(fontsize = 15),
  axis_name_gp = gpar(fontsize = 6),
  column_title = "NAcc_pos",
  width = 1,
  axis_name = c("-15000", "start", "end", "15000"),
  split = split,
  use_raster = TRUE,
  raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["NA"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCT (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCT (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["NA"]],
    col = mCH_col,
    name = "mCA (-)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (-)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["NA"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCA (+)"]][["BA9"]],
    col = mCH_col,
    name = "mCA (+)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["NA"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "NAcc_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG") +
  EnrichedHeatmap(
    mat = list_of_CH_DMRs_centered_mat[["mCT (-)"]][["mCG (S)"]][["BA9"]],
    col = mCG_col,
    name = "mCG (S)",
    pos_line = FALSE,
    column_title_gp = gpar(fontsize = 15),
    axis_name_gp = gpar(fontsize = 6),
    column_title = "BA9_pos",
    width = 1,
    axis_name = c("-15000", "start", "end", "15000"),
    split = split,
    use_raster = TRUE,
    raster_device = "CairoPNG")
dev.off()

pdf("../figures/neg_CT-DMRs_centered.clusters.pdf")
for (i in 1:4) {
  plotManyRegions2(
    list_of_CG_BSseq = list_of_CG_BSseq,
    list_of_CH_BSseq = list_of_CH_BSseq,
    regions = list_of_CH_DMRs[["mCT (-)"]][split == i][1:20],
    extend = 15000,
    list_of_CG_addRegions = list("mCG (S)" = CG_DMRs, "mCG (L)" = CG_blocks),
    list_of_CH_addRegions = list_of_CH_DMRs,
    geneTrack = RefSeq_exons3,
    addNames = TRUE,
    addTicks = FALSE,
    verbose = TRUE)
}
dev.off()

### ============================================================================
### TODO
###

# TODO: ATAC
# TODO: Threshold mCH and ATAC based on quantiles (try 99.5%)
# TODO: Put mCA and mCT on common scale?
