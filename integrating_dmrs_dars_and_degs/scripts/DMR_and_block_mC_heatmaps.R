# Heatmap of average mC in CG-DMRs, CG-blocks, and CH-DMRs
# Peter Hickey
# 2017-12-11

library(SummarizedExperiment)
library(dplyr)
library(scales)
library(ComplexHeatmap)
library(matrixStats)
library(png)

list_of_SEs <- readRDS("../objects/list_of_SEs.DMRs_and_blocks.rds")
list_of_sample_level_SEs <- lapply(list_of_SEs, "[[", "sample_level")
list_of_condition_level_SEs <- lapply(list_of_SEs, "[[", "condition_level")

# Only retain those with FWER <= 0.05
list_of_sample_level_SEs <- lapply(list_of_sample_level_SEs, function(se) {
  if (!is.null(rowRanges(se)$fwer)) {
    return(se[rowRanges(se)$fwer <= 50])
  } else {
    se
  }
})
list_of_condition_level_SEs <-
  lapply(list_of_condition_level_SEs, function(se) {
    if (!is.null(rowRanges(se)$fwer)) {
      return(se[rowRanges(se)$fwer <= 50])
    } else {
      se
    }
  })

z_sample_level <- do.call(
  rbind, lapply(list_of_sample_level_SEs, function(se) {
    # Z-score
    a <- do.call(cbind, lapply(assays(se), function(a) {
      (a - rowMeans2(a)) / rowSds(a)
    }))
    colnames(a) <- paste0(colnames(a), ".",
                          rep(assayNames(se), each = ncol(se)))
    rownames(a) <- NULL
    a
  }))

z_condition_level <-
  do.call(rbind, lapply(list_of_condition_level_SEs, function(se) {
    # Z-score
    a <- do.call(cbind, lapply(assays(se), function(a) {
      (a - rowMeans2(a)) / rowSds(a)
    }))
    colnames(a) <- paste0(colnames(a), ".",
                          rep(assayNames(se), each = ncol(se)))
    rownames(a) <- NULL
    a
  }))

se_sample_level <- SummarizedExperiment(
  assays = z_sample_level,
  rowRanges = unlist(GRangesList(lapply(list_of_sample_level_SEs, function(x) {
    granges(rowRanges(x))
  })), use.names = FALSE))
colData(se_sample_level) <- DataFrame(
  Donor = substr(colnames(se_sample_level), 1, 4),
  Strand = ifelse(grepl("+", colnames(se_sample_level), fixed = TRUE),
                  "+",
                  ifelse(grepl("-", colnames(se_sample_level), fixed = TRUE),
                         "-",
                         "*")),
  Context = ifelse(grepl("CA", colnames(se_sample_level)),
                   "CA",
                   ifelse(grepl("CT", colnames(se_sample_level)),
                          "CT",
                          "CG")),
  Tissue = sapply(strsplit(colnames(se_sample_level), "_", fixed = TRUE),
                  "[[", 2))

se_condition_level <- SummarizedExperiment(
  assays = z_condition_level,
  rowRanges = unlist(GRangesList(
    lapply(list_of_condition_level_SEs, function(x) {
      granges(rowRanges(x))
    })), use.names = FALSE))
colData(se_condition_level) <- DataFrame(
  Strand = ifelse(grepl("+", colnames(se_condition_level), fixed = TRUE),
                  "+",
                  ifelse(grepl("-", colnames(se_condition_level),
                               fixed = TRUE),
                         "-",
                         "*")),
  Context = ifelse(grepl("CA", colnames(se_condition_level)),
                   "CA",
                   ifelse(grepl("CT", colnames(se_condition_level)),
                          "CT",
                          "CG")),
  Tissue = sapply(strsplit(colnames(se_condition_level), ".", fixed = TRUE),
                  "[[", 1))


# Plot
myPlot <- function(feature, se, km = 1L, feature2 = NULL) {
  i <- grep(feature, rownames(se), fixed = TRUE)
  # i <- sort(sample(i, min(1000, length(i))))
  j <- seq_len(ncol(se))
  se_ <- se[i, j]

  if (is.null(se$Donor)) {
    top_anno <- HeatmapAnnotation(
      df = data.frame(
        Strand = se$Strand,
        Context = se$Context,
        Tissue = se$Tissue),
      col = list(Strand = c("+" = "black",
                            "-" = grey(0.9),
                            "*" = grey(0.5)),
                 Context = c("CA" = "#66a61e",
                             "CT" = "#e6ab02",
                             "CG" = "darkviolet"),
                 Tissue = c("BA24" = "deeppink",
                            "BA9" = "deepskyblue",
                            "HC" = "darkgrey",
                            "NA" = "chocolate1")))
    col <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
    row_titles <- c("mCG (S)" = "CG-DMRs",
                    "mCG (L)" = "CG-blocks",
                    "mCA (+)" = "CA(+)-DMRs",
                    "mCA (-)" = "CA(-)-DMRs",
                    "mCT (+)" = "CT(+)-DMRs",
                    "mCT (-)" = "CT(-)-DMRs")
  } else {
    top_anno <- HeatmapAnnotation(
      df = data.frame(
        Strand = se$Strand,
        Context = se$Context,
        Donor = se$Donor,
        Tissue = se$Tissue),
      col = list(Strand = c("+" = "black",
                            "-" = grey(0.9),
                            "*" = grey(0.5)),
                 Context = c("CA" = "#66a61e",
                             "CT" = "#e6ab02",
                             "CG" = "darkviolet"),
                 Donor = c("5248" = "aquamarine",
                           "5284" = "chocolate1",
                           "5552" = "firebrick3",
                           "5569" = "deepskyblue",
                           "5570" = "purple",
                           "5628" = "darkgrey"),
                 Tissue = c("BA24" = "deeppink",
                            "BA9" = "deepskyblue",
                            "HC" = "darkgrey",
                            "NA" = "chocolate1")))
    col <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
    row_titles <- c("mCG (S)" = "CG-DMRs",
                    "mCG (L)" = "CG-blocks",
                    "mCA (+)" = "CA(+)-DMRs",
                    "mCA (-)" = "CA(-)-DMRs",
                    "mCT (+)" = "CT(+)-DMRs",
                    "mCT (-)" = "CT(-)-DMRs")
  }

  col <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
  row_titles <- c("mCG (S)" = "CG-DMRs",
                  "mCG (L)" = "CG-blocks",
                  "mCA (+)" = "CA(+)-DMRs",
                  "mCA (-)" = "CA(-)-DMRs",
                  "mCT (+)" = "CT(+)-DMRs",
                  "mCT (-)" = "CT(-)-DMRs")
  row_title <- paste0(row_titles[feature], " (n = ", length(i), ")")

  if (!is.null(feature2)) {
    split <- overlapsAny(se_, feature2)
  } else {
    split <- NULL
  }

  set.seed(666)
  # TODO: `raster` (can't install 'png' package on my laptop)
  Heatmap(matrix = assay(se_),
          col = col,
          name = "mC Z-score",
          na_col = "black",
          row_title = row_title,
          column_title = "Samples",
          show_row_dend = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          top_annotation = top_anno,
          km = km,
          split = split,
          use_raster = TRUE)
}

pdf("../figures/mC.heatmap.condition_level.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot,
       se = se_condition_level,
       feature2 = list_of_SEs[["mCG (S)"]][["sample_level"]])
dev.off()

pdf("../figures/mC.heatmap.sample_level.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, se = se_sample_level[, se_sample_level$Context != "CG"])
dev.off()

pdf("../figures/mC.heatmap.condition_level.overlaps_CG-DMR.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot,
       se = subsetByOverlaps(se_condition_level, list_of_SEs[["mCG (S)"]][["sample_level"]]),
       feature2 = list_of_SEs[["mCG (S)"]][["sample_level"]])
dev.off()

pdf("../figures/mC.heatmap.sample_level.overlaps_CG-DMR.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot,
       se = subsetByOverlaps(se_sample_level, list_of_SEs[["mCG (S)"]][["sample_level"]]))
dev.off()


# TODO: Get coordinates of regions that 'flip' and plotManyRegions2(). Check
#       whether it's a big change in mC or a small artefact
