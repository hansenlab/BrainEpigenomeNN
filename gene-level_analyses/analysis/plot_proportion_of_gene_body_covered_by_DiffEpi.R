# Plot proportion of gene body covered by DiffEpi
# Peter Hickey
# 2018-03-08

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(broom)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "number_of_bases_of_gene_body_covered_by_DiffEpi.rds"))

### ----------------------------------------------------------------------------
### Data wrangling
###

sx <- inner_join(
  x,
  select(genes_df, gene, width, logFC, DE, gene_type)) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)", "CH-DMRs"))) %>%
  # NOTE: Not using bigDARs, CA-DMRs (opposite strand), or CH-DMRs
  filter(!DiffEpi %in%
           c("bigDARs", "CA-DMRs (opposite strand)", "CH-DMRs")) %>%
  # NOTE: Aggregate CA-DMRs (same strand)
  group_by(gene, DiffEpi, width, logFC, DE, gene_type) %>%
  summarise(n_bases_covered = min(sum(n_bases_covered), width)) %>%
  mutate(proportion_covered = n_bases_covered / width)

# NOTE: Colours are from Dark2 Colour Brewer palette
#       (http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=8) manually
#       chosen to match those in LDSC figures (CA-DMRs coloured as CH-DMRs)
cols <- c("CG-DMRs" = "#1b9e77",
          "CA-DMRs" = "#d95f02",
          "CA-DMRs (opposite strand)" = "#d95f02",
          "DARs" = "#e7298a",
          "CG-blocks" = "dodgerBlue")

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# Proportion of gene body covered by DiffEpi vs. |logFC| for genes with > 0 %
# of gene body covered by DiffEpi
#

g <- ggplot(sx, aes(x = proportion_covered, y = abs(logFC))) +
  facet_grid(DiffEpi ~ gene_type, margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              se = TRUE)
save_plot(
  "../figures/proportion_of_gene_body_covered_by_DiffEpi_vs_abs_logFC.pdf",
  g,
  base_height = 10)

# ------------------------------------------------------------------------------
# Proportion of genes differentially expressed by % of gene covered by DiffEpi
#

# NOTE: Adapted from
#       http://www.cookbook-r.com/Statistical_analysis/Logistic_regression/

g <- ggplot(sx, aes(x = proportion_covered, y = as.integer(DE))) +
  facet_grid(DiffEpi ~ gene_type,
             margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              method.args = list(family = "binomial"),
              se = TRUE,
              lty = 2) +
  ylab("Proportion of genes DE") +
  ylim(c(0, 1))
save_plot(
  "../figures/proportion_of_gene_covered_by_DiffEpi_vs_proportion_of_genes_DE.GLM_and_GAM.pdf",
  g,
  base_height = 10)

g <- ggplot(sx,
            aes(x = proportion_covered, y = as.integer(DE), col = DiffEpi,
                fill = DiffEpi)) +
  facet_grid( ~ gene_type,
              margins = "gene_type") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              method.args = list(family = "binomial"),
              se = TRUE,
              lty = 2) +
  ylab("Proportion of genes DE") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ylim(c(0, 1))
save_plot(
  "../figures/proportion_of_gene_covered_by_DiffEpi_vs_proportion_of_genes_DE.GAM.pdf",
  g,
  base_height = 10)

g <- ggplot(filter(sx, gene_type == "PC"),
            aes(x = proportion_covered, y = as.integer(DE), col = DiffEpi,
                fill = DiffEpi)) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              method.args = list(family = "binomial"),
              se = TRUE,
              lty = 2) +
  ylab("Proportion of genes DE") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ylim(c(0, 1))
save_plot(
  "../figures/proportion_of_gene_covered_by_DiffEpi_vs_proportion_of_genes_DE.GAM.PC_only.pdf",
  g,
  base_height = 10)
