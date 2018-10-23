# Heatmap of average mCH in CH-DMRs
# Peter Hickey
# 2017-12-05

library(SummarizedExperiment)
library(dplyr)
library(scales)
library(gplots)
library(matrixStats)

# TODO: Make versions at condition-level
list_of_SEs <- readRDS("../objects/list_of_SEs.DMRs_and_blocks.rds")
list_of_SEs <- lapply(list_of_SEs, "[[", "sample_level")

# TODO: CG-DMRs use FWER <= 0.05 rather than FWER < 0.05
# Only retain those with FWER < 0.05
list_of_SEs <- lapply(list_of_SEs, function(se) {
  if (!is.null(rowRanges(se)$fwer)) {
    return(se[rowRanges(se)$fwer <= 50])
  } else {
    se
  }
})

a <- do.call(rbind, lapply(list_of_SEs, function(se) {
  # Z-score
  a <- do.call(cbind, lapply(assays(se), function(a) {
    (a - rowMeans2(a)) / rowSds(a)
  }))
  # a <- do.call(cbind, assays(se))
  colnames(a) <- paste0(colnames(a), ".", rep(assayNames(se), each = ncol(se)))
  rownames(a) <- NULL
  a
}))

se <- SummarizedExperiment(
  assays = a,
  rowRanges = unlist(GRangesList(lapply(list_of_SEs, function(x) {
    granges(rowRanges(x))
  })), use.names = FALSE))


# Plot
col <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

myPlot <- function(feature, plot = 1, mCH_only = FALSE) {
  i <- grep(feature, rownames(se), fixed = TRUE)
  # i <- seq_len(nrow(se))
  i <- sort(sample(i, min(1000, length(i))))
  # j <- grep("mCA", colnames(se), invert = FALSE)
  # j <- grep("mCG", colnames(se), invert = TRUE)
  if (mCH_only) {
    j <- grep("mCG", colnames(se), invert = TRUE)
  } else {
    j <- seq_len(ncol(se))
  }
  se_ <- se[i, j]

  strand_cols <- case_when(grepl("mCA \\(\\+", colnames(se_)) ~
                             "#66a61e",
                           grepl("mCA \\(\\-", colnames(se_)) ~
                             alpha("#66a61e", 0.5),
                           grepl("mCT \\(\\+", colnames(se_)) ~
                             "#e6ab02",
                           grepl("mCT \\(\\-", colnames(se_)) ~
                             alpha("#e6ab02", 0.5),
                           grepl("mCG \\(S\\)", colnames(se_)) ~ "black",
                           grepl("mCG \\(L\\)", colnames(se_)) ~ "grey")
  tissue_cols <- case_when(grepl("BA24", colnames(se_)) ~ "deeppink",
                           grepl("BA9", colnames(se_)) ~ "deepskyblue",
                           grepl("HC", colnames(se_)) ~ "darkgrey",
                           grepl("NA", colnames(se_)) ~ "chocolate1")
  strand_cols2 <- case_when(grepl("mCA \\(\\+", rownames(se_)) ~
                              "#66a61e",
                            grepl("mCA \\(\\-", rownames(se_)) ~
                              alpha("#66a61e", 0.5),
                            grepl("mCT \\(\\+", rownames(se_)) ~
                              "#e6ab02",
                            grepl("mCT \\(\\-", rownames(se_)) ~
                              alpha("#e6ab02", 0.5),
                            grepl("mCG \\(S\\)", rownames(se_)) ~ "black",
                            grepl("mCG \\(L\\)", rownames(se_)) ~ "grey")

  # NOTE: Setting seed to ensure reproducibility of clusters
  set.seed(666)
  if (plot == 1) {
    heatmap.2(x = assay(se_),
              Rowv = FALSE,
              dendrogram = "column",
              trace = "none",
              labRow = NA,
              margins = c(10, 5),
              ColSideColors = tissue_cols,
              RowSideColors = strand_cols2,
              col = col,
              main = feature)
  } else if (plot == 2) {
    heatmap.2(x = assay(se_),
              Rowv = TRUE,
              dendrogram = "column",
              trace = "none",
              labRow = NA,
              margins = c(10, 5),
              ColSideColors = tissue_cols,
              RowSideColors = strand_cols2,
              col = col,
              main = feature)
  } else if (plot == 3) {
    heatmap.2(x = assay(se_),
              Rowv = TRUE,
              dendrogram = "column",
              trace = "none",
              labRow = NA,
              margins = c(10, 5),
              ColSideColors = strand_cols,
              RowSideColors = strand_cols2,
              col = col,
              main = feature)
  }
}

pdf("../figures/mC.heatmap.plot1.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot)
dev.off()
pdf("../figures/mC.heatmap.plot1.mCH_only.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, mCH_only = TRUE)
dev.off()

pdf("../figures/mC.heatmap.plot2.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, plot = 2)
dev.off()
pdf("../figures/mC.heatmap.plot2.mCH_only.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, plot = 2, mCH_only = TRUE)
dev.off()

pdf("../figures/mC.heatmap.plot3.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, plot = 3)
dev.off()
pdf("../figures/mC.heatmap.plot3.mCH_only.pdf",
    height = 7,
    width = 7)
lapply(c("mCG (S)", "mCG (L)", "mCA (+)", "mCA (-)", "mCT (+)", "mCT (-)"),
       myPlot, plot = 3, mCH_only = TRUE)
dev.off()

# TODO: Handle overlapping CH-DMRs somehow
# TODO: Get coordinates of regions that 'flip' and plotManyRegions2(). Check
#       whether it's a big change in mC or a small artefact
