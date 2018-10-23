# Plot presence of DiffEpi around genes using equal number of bins along gene
# body +/- 2 gene body equivalents (GBE)
# Peter Hickey
# 2018-03-08

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)
library(broom)
library(tidyr)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "presence_of_DiffEpi_around_scaled_gene.rds"))

### ----------------------------------------------------------------------------
### Plots
###

g <- x %>%
  inner_join(select(genes_df, gene, DE, CGI_promoter, gene_type)) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  group_by(DiffEpi, DE, CGI_promoter, bin, gene_type) %>%
  summarise(proportion = sum(presence_of_DiffEpi) /
              length(presence_of_DiffEpi)) %>%
  ungroup() %>%
  # NOTE: Not using bigDARs or CA-DMRs (opposite strand)
  filter(!DiffEpi %in% c("bigDARs", "CA-DMRs (opposite strand)")) %>%
  ggplot(aes(x = bin, y = proportion, col = DE)) +
  # facet_grid(DiffEpi ~ gene_type + CGI_promoter,
  #            scales = "free",
  #            labeller = labeller(CGI_promoter = label_both),
  #            margins = "gene_type") +
  facet_grid(DiffEpi ~ gene_type + CGI_promoter,
             scales = "free",
             labeller = labeller(CGI_promoter = label_both)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 10, lty = 2) +
  geom_vline(xintercept = 15, lty = 2) +
  geom_step(direction = "vh") +
  ylab("Proportion of genes") +
  scale_colour_manual(values = c("FALSE" = "black",
                                 "TRUE" = "orange"))
save_plot(
  "../figures/presence_of_DiffEpi_around_scaled_gene.pdf",
  g,
  base_height = 10)

g <- x %>%
  inner_join(select(genes_df, gene, DE, gene_type)) %>%
  group_by(DiffEpi, DE, bin, gene_type) %>%
  summarise(proportion = sum(presence_of_DiffEpi) /
              length(presence_of_DiffEpi)) %>%
  ungroup() %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  # NOTE: Not using bigDARs or CA-DMRs (opposite strand)
  filter(!DiffEpi %in% c("bigDARs", "CA-DMRs (opposite strand)")) %>%
  ggplot(aes(x = bin, y = proportion, col = DE)) +
  facet_grid(DiffEpi ~ gene_type,
             scales = "free",
             margins = "gene_type") +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_vline(xintercept = 300, lty = 2) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  ylab("Proportion of genes") +
  scale_colour_manual(values = c("FALSE" = "black",
                                 "TRUE" = "orange"))
save_plot(
  "../figures/presence_of_DiffEpi_around_scaled_gene.no_stratification_by_CGI-promoter.pdf",
  g,
  base_height = 10)

### ----------------------------------------------------------------------------
### Tables
###

x %>%
  filter(bin >= 200, bin <= 300) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  # NOTE: Not using bigDARs or CA-DMRs (opposite strand)
  filter(!DiffEpi %in% c("bigDARs", "CA-DMRs (opposite strand)")) %>%
  inner_join(select(genes_df, gene, DE, CGI_promoter, gene_type)) %>%
  group_by(DiffEpi, CGI_promoter, gene_type, DE, gene) %>%
  summarise(presence_of_DiffEpi = any(presence_of_DiffEpi)) %>%
  group_by(DiffEpi, CGI_promoter, gene_type, DE) %>%
  summarise(x = sum(presence_of_DiffEpi),
            n = length(presence_of_DiffEpi))
