# Plot average mC levels at gene features against gene expression data
# Peter Hickey
# 2018-03-07

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
    "mC_at_gene_features.rds"))

### ----------------------------------------------------------------------------
### Functions
###

tidyCor <- function(x, y) {
  val <- try(cor.test(x, y, na.action = "na.omit"), silent = TRUE)
  if (is(val, "try-error")) {
    return(data_frame(estimate = NA_real_))
  }
  tidy(val)
}

### ----------------------------------------------------------------------------
### Data wrangling
###

# Create the grouped data
grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC, DE)) %>%
  group_by(context, gene, gene_type, logFC, DE, feature) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, DE, feature, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = as.character(DE))

# Create the partially grouped data
gene_type_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC)) %>%
  group_by(context, gene, gene_type, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = "(all)",
         feature = "(all)") %>%
  select(gene_type, DE, feature, everything())
DE_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, DE, logFC)) %>%
  group_by(context, gene, DE, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(DE, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         feature = "(all)",
         DE = as.character(DE)) %>%
  select(gene_type, DE, feature, everything())
feature_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, logFC)) %>%
  group_by(context, gene, feature, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(feature, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         DE = "(all)") %>%
  select(gene_type, DE, feature, everything())
gene_type_and_DE_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, DE, logFC)) %>%
  group_by(context, gene, gene_type, DE, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, DE, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = as.character(DE),
         feature = "(all)") %>%
  select(gene_type, DE, feature, everything())
gene_type_and_feature_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC)) %>%
  group_by(context, gene, gene_type, feature, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, feature, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = "(all)") %>%
  select(gene_type, DE, feature, everything())
DE_and_feature_grouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, DE, logFC)) %>%
  group_by(context, gene, DE, feature, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(DE, feature, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         DE = as.character(DE)) %>%
  select(gene_type, DE, feature, everything())

# Create the ungrouped data
ungrouped_df <- filter(x, !is.na(mC)) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, logFC)) %>%
  group_by(context, gene, logFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         DE = "(all)",
         feature = "(all)") %>%
  select(gene_type, DE, feature, everything())

# bind grouped, partially grouped, and ungrouped together with factors to
# order the data
sx <- bind_rows(grouped_df,
                gene_type_grouped_df, DE_grouped_df, feature_grouped_df,
                gene_type_and_DE_grouped_df, gene_type_and_feature_grouped_df,
                DE_and_feature_grouped_df,
                ungrouped_df) %>%
  mutate(gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")),
         DE = factor(DE,
                     levels = c("FALSE", "TRUE", "(all)")),
         feature = factor(feature,
                          levels = c("promoters", "first_exon", "first_intron",
                                     "remaining_exons", "remaining_introns",
                                     "(all)")))

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# NAcc_pos
#

# NOTE: Omitting 1% of genes with largest RPKM
g <- x %>%
  filter(condition == "NAcc_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_NAcc_pos)) %>%
  group_by(condition, context, gene, feature, gene_type, rpkm_NAcc_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(rpkm_NAcc_pos < quantile(rpkm_NAcc_pos, .99, na.rm = TRUE)) %>%
  ggplot(aes(x = mC, y = rpkm_NAcc_pos, col = feature)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_smooth(method = "loess", se = TRUE) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/average_mC_at_gene_features_vs_rpkm.NAcc_pos.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos
#

# NOTE: Omitting 1% of genes with largest RPKM
g <- x %>%
  filter(condition == "BA9_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_BA9_pos)) %>%
  group_by(condition, context, gene, feature, gene_type, rpkm_BA9_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(rpkm_BA9_pos < quantile(rpkm_BA9_pos, .99, na.rm = TRUE)) %>%
  ggplot(aes(x = mC, y = rpkm_BA9_pos, col = feature)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_smooth(method = "loess", se = TRUE) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/average_mC_at_gene_features_vs_rpkm.BA9_pos.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos vs. BA9_pos
#

g <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "CpG",
           grepl("opposite", context) ~ "CpA (opposite strand)",
           TRUE ~ "CpA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC, DE)) %>%
  group_by(context, gene, gene_type, feature, logFC, DE) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  ggplot(aes(x = delta_mC, y = logFC, col = DE)) +
  facet_grid(gene_type + feature ~ context,
             scales = "free",
             margins = "gene_type") +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = delta_mC, y = logFC),
              method = "loess",
              col = "dodgerBlue",
              se = FALSE,
              lty = 2) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("black", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/average_mC_at_gene_features_vs_logFC.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos vs. BA9_pos correlation plot
#

g <- sx %>%
  ggplot(aes(x = feature, y = estimate, col = DE)) +
  facet_grid(gene_type ~ context) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  scale_y_continuous(limits = c(-1, 0)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "orange", `(all)` = "dodgerBlue")) +
  ylab("Pearson correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/cor_average_mC_at_gene_features_and_logFC.pdf",
          g,
          base_height = 10)
