# Create a dataset matrix
# Peter Hickey
# 2018-03-19

# ------------------------------------------------------------------------------
# Setup
#

library(SummarizedExperiment)
library(HDF5Array)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

extdir <- "../extdata"

# ------------------------------------------------------------------------------
# Load data
#

CG_BSseq_unsorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))
CG_BSseq_small_smooth <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
ATAC_CD <-
  readRDS("../../ATAC-seq/objects/colData-flow-sorted-brain-atac-seq.rds")
RNA_CD <- readRDS(file.path(extdir, "flow-sorted-brain-rna-seq", "objects",
                            "colData-flow-sorted-brain-rna-seq.rds"))

# ------------------------------------------------------------------------------
# Create data frames
#

datasets <- bind_rows(
  as.data.frame(
    colData(CG_BSseq_unsorted)) %>%
    select(Individual, Tissue) %>%
    mutate(
      Assay = "WGBS",
      NeuN = "Unsorted") %>%
    filter(Tissue != "caudate"),
  as.data.frame(
    colData(CG_BSseq_small_smooth)) %>%
    select(Individual, Tissue, NeuN) %>%
    mutate(Assay = "WGBS"),
  as.data.frame(ATAC_CD) %>%
    filter(REPLICATE == "rep1") %>%
    select(DONOR, TISSUE, NEUN) %>%
    rename(Individual = DONOR, Tissue = TISSUE, NeuN = NEUN) %>%
    mutate(Assay = "ATAC-seq"),
  as.data.frame(RNA_CD) %>%
    select(DONOR, TISSUE, NEUN) %>%
    rename(Individual = DONOR, Tissue = TISSUE, NeuN = NEUN) %>%
    mutate(Assay = "RNA-seq")) %>%
  mutate(
    Tissue = gsub("^NA$", "NAcc", Tissue),
    Tissue = factor(Tissue, levels = c("BA9", "BA24", "HC", "NAcc")),
    NeuN = case_when(
      NeuN == "pos" ~ "NeuN+",
      NeuN == "neg" ~ "NeuN-",
      NeuN == "Unsorted" ~ "Unsorted"),
    NeuN = factor(NeuN, levels = c("NeuN+", "NeuN-", "Unsorted"))) %>%
  as_data_frame()

tissue_colours <- as.data.frame(
  colData(CG_BSseq_small_smooth)) %>%
  select(Tissue, Tissue_color) %>%
  distinct() %>%
  rename(Colour = Tissue_color) %>%
  mutate(
    Tissue = gsub("^NA$", "NAcc", Tissue),
    Tissue = factor(Tissue, levels = c("BA9", "BA24", "HC", "NAcc")))

neun_colours <- data_frame(
  NeuN = factor(c("NeuN+", "NeuN-", "Unsorted")),
  Colour = c("darkgreen", "purple", "black"))

# ------------------------------------------------------------------------------
# Plot
#

individual_wide <- datasets %>%
  mutate(Assayed = TRUE) %>%
  complete(Individual, Assay, Tissue, NeuN, fill = list(Assayed = FALSE)) %>%
  group_by(Assay, Tissue, NeuN) %>%
  filter(any(Assayed)) %>%
  ungroup() %>%
  mutate(Assayed = as.integer(Assayed)) %>%
  spread(-Assayed, Assayed)

individual_wide_matrix <- select(individual_wide, -Assay, -Tissue, -NeuN) %>%
  as.matrix()
rownames(individual_wide_matrix) <- select(
  individual_wide, Assay, Tissue, NeuN) %>%
  mutate(Rowname = paste0(Assay, ": ", Tissue, " (", NeuN, ")")) %>%
  pull(Rowname)

tissue_neun_wide <- datasets %>%
  select(Assay, Tissue, NeuN) %>%
  distinct() %>%
  mutate(Condition = paste0(Tissue, " (", NeuN, ")")) %>%
  select(-Tissue, -NeuN) %>%
  mutate(Assayed = TRUE) %>%
  complete(Condition, Assay, fill = list(Assayed = FALSE)) %>%
  mutate(Assayed = as.integer(Assayed)) %>%
  spread(-Assayed, Assayed)

tissue_neun_wide_matrix <- select(tissue_neun_wide, -Assay) %>%
  as.matrix() %>%
  t()
colnames(tissue_neun_wide_matrix) <- pull(tissue_neun_wide, Assay)

ra1 <- HeatmapAnnotation(
  df = as.data.frame(select(individual_wide, Assay)),
  which = "row")
ra2 <- HeatmapAnnotation(
  df = select(individual_wide, Tissue, NeuN) %>%
    rename(`Brain Region` = Tissue) %>%
    as.data.frame(),
  col = list(
    `Brain Region` = setNames(tissue_colours$Colour, tissue_colours$Tissue),
    NeuN = setNames(neun_colours$Colour, neun_colours$NeuN)),
  which = "row")
hm <- Heatmap(
  matrix = individual_wide_matrix,
  col = c("grey", "darkblue"),
  name = "Assayed",
  rect_gp = gpar(col = "white"),
  column_title = "Donor",
  row_title = "Assay: Brain Region (NeuN-status)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  heatmap_legend_param = list(labels = c("False", "True")))

pdf("../figures/individual-level_dataset_matrix.pdf")
hm + ra1 + ra2
dev.off()

ra <- HeatmapAnnotation(
  df = data_frame(
    `Brain Region` =
      sapply(strsplit(rownames(tissue_neun_wide_matrix), " "), "[[", 1),
    NeuN = gsub(
      "\\)",
      "",
      gsub(
        "\\(",
        "",
        sapply(strsplit(rownames(tissue_neun_wide_matrix), " "), "[[", 2)))) %>%
    as.data.frame(),
  col = list(
    `Brain Region` = setNames(tissue_colours$Colour, tissue_colours$Tissue),
    NeuN = setNames(neun_colours$Colour, neun_colours$NeuN)),
  which = "row")
hm <- Heatmap(
  matrix = tissue_neun_wide_matrix,
  col = c("grey", "darkblue"),
  name = "Assayed",
  rect_gp = gpar(col = "white"),
  column_title = "Assay",
  row_title = "Brain Region (NeuN-status)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = c(3, 1, 2),
  row_names_side = "left",
  column_names_side = "top",
  heatmap_legend_param = list(labels = c("False", "True")))

pdf("../figures/condition-level_dataset_matrix.pdf")
hm + ra
dev.off()
