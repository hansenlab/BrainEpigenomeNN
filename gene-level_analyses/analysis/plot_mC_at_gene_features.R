# Plot mC levels at gene features
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

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "mC_at_gene_features.rds"))

### ----------------------------------------------------------------------------
### Data wrangling
###

genes_df <- genes_df %>%
  mutate(qrpkm_NAcc_pos = cut(rpkm_NAcc_pos,
                              quantile(rpkm_NAcc_pos, 0:4 / 4, na.rm = TRUE),
                              ordered_result = TRUE,
                              labels = paste0("Q", 1:4)),
         qrpkm_BA9_pos = cut(rpkm_BA9_pos,
                             quantile(rpkm_BA9_pos, 0:4 / 4, na.rm = TRUE),
                             ordered_result = TRUE,
                             labels = paste0("Q", 1:4)),
         qLogFC = cut(logFC,
                      quantile(logFC, 0:4 / 4, na.rm = TRUE),
                      ordered_result = TRUE,
                      labels = paste0("Q", 1:4)))

x <- x %>%
  mutate(feature = factor(feature,
                          levels = c("promoters", "first_exon", "first_intron",
                                     "remaining_introns", "remaining_exons")))

# Create the grouped data
grouped_df <- x %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type))) %>%
  group_by(CGI_promoter, gene_type, condition, context, feature) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CGI_promoter = as.character(CGI_promoter))
NAcc_pos_rpkm_grouped_df <- x %>%
  filter(condition == "NAcc_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type,
                  qrpkm_NAcc_pos),
           !is.na(qrpkm_NAcc_pos))) %>%
  group_by(CGI_promoter, gene_type, condition, context, feature,
           qrpkm_NAcc_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CGI_promoter = as.character(CGI_promoter))
BA9_pos_rpkm_grouped_df <- x %>%
  filter(condition == "BA9_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type,
                  qrpkm_BA9_pos),
           !is.na(qrpkm_BA9_pos))) %>%
  group_by(CGI_promoter, gene_type, condition, context, feature,
           qrpkm_BA9_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CGI_promoter = as.character(CGI_promoter))
delta_mC_grouped_df <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos) %>%
  select(-NAcc_pos, -BA9_pos) %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type))) %>%
  group_by(CGI_promoter, gene_type, context, feature) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CGI_promoter = as.character(CGI_promoter))
logFC_delta_mC_grouped_df <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos) %>%
  select(-NAcc_pos, -BA9_pos) %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type,
                  qLogFC),
           !is.na(qLogFC))) %>%
  group_by(CGI_promoter, gene_type, context, feature, qLogFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CGI_promoter = as.character(CGI_promoter))

# Create the ungrouped data
ungrouped_df <- x %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene))) %>%
  group_by(condition, context, feature) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  mutate(gene_type = "(all)",
         CGI_promoter = "(all)") %>%
  select(CGI_promoter, gene_type, everything())
NAcc_pos_rpkm_ungrouped_df <- x %>%
  filter(condition == "NAcc_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, qrpkm_NAcc_pos),
           !is.na(qrpkm_NAcc_pos))) %>%
  group_by(condition, context, feature, qrpkm_NAcc_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  mutate(gene_type = "(all)",
         CGI_promoter = "(all)") %>%
  select(CGI_promoter, gene_type, everything())
BA9_pos_rpkm_ungrouped_df <- x %>%
  filter(condition == "BA9_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, qrpkm_BA9_pos),
           !is.na(qrpkm_BA9_pos))) %>%
  group_by(condition, context, feature, qrpkm_BA9_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  mutate(gene_type = "(all)",
         CGI_promoter = "(all)") %>%
  select(CGI_promoter, gene_type, everything())
delta_mC_ungrouped_df <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos) %>%
  select(-NAcc_pos, -BA9_pos) %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, CGI_promoter, gene_type))) %>%
  group_by(context, feature) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  mutate(gene_type = "(all)",
         CGI_promoter = "(all)") %>%
  select(CGI_promoter, gene_type, everything())
logFC_delta_mC_ungrouped_df <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos) %>%
  select(-NAcc_pos, -BA9_pos) %>%
  mutate(context = case_when(
    context == "CpG" ~ "CpG",
    grepl("opposite", context) ~ "CpA (opposite strand)",
    TRUE ~ "CpA (same strand)")) %>%
  inner_join(
    filter(select(genes_df, gene, qLogFC),
           !is.na(qLogFC))) %>%
  group_by(context, feature, qLogFC) %>%
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  mutate(gene_type = "(all)",
         CGI_promoter = "(all)") %>%
  select(CGI_promoter, gene_type, everything())

# bind grouped and ungrouped together with factors to order the data
sx <- bind_rows(grouped_df, ungrouped_df) %>%
  mutate(CGI_promoter = factor(CGI_promoter,
                               levels = c("FALSE", "TRUE", "(all)")),
         gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))
NAcc_pos_rpkm_sx <- bind_rows(
  NAcc_pos_rpkm_grouped_df, NAcc_pos_rpkm_ungrouped_df) %>%
  mutate(CGI_promoter = factor(CGI_promoter,
                               levels = c("FALSE", "TRUE", "(all)")),
         gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))
BA9_pos_rpkm_sx <- bind_rows(
  BA9_pos_rpkm_grouped_df, BA9_pos_rpkm_ungrouped_df) %>%
  mutate(CGI_promoter = factor(CGI_promoter,
                               levels = c("FALSE", "TRUE", "(all)")),
         gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))
delta_mC_sx <- bind_rows(
  delta_mC_grouped_df, delta_mC_ungrouped_df) %>%
  mutate(CGI_promoter = factor(CGI_promoter,
                               levels = c("FALSE", "TRUE", "(all)")),
         gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))
logFC_delta_mC_sx <- bind_rows(
  logFC_delta_mC_grouped_df, logFC_delta_mC_ungrouped_df) %>%
  mutate(CGI_promoter = factor(CGI_promoter,
                               levels = c("FALSE", "TRUE", "(all)")),
         gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# NAcc_pos
#

g <- sx %>%
  filter(condition == "NAcc_pos") %>%
  ggplot(aes(x = feature, y = mC, group = condition)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             scales = "free_y",
             labeller = labeller(CGI_promoter = label_both)) +
  ylim(0, NA) +
  ggtitle("NAcc_pos") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.NAcc_pos.pdf",
          g,
          base_height = 10)

g <- NAcc_pos_rpkm_sx %>%
  ggplot(aes(x = feature, y = mC, col = qrpkm_NAcc_pos,
             group = qrpkm_NAcc_pos)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             scales = "free_y",
             labeller = labeller(CGI_promoter = label_both)) +
  ylim(0, NA) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("NAcc_pos") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.NAcc_pos.rpkm.pdf",
          g,
          base_height = 7,
          base_width = 14)

# ------------------------------------------------------------------------------
# BA9_pos
#

g <- sx %>%
  filter(condition == "BA9_pos") %>%
  ggplot(aes(x = feature, y = mC, group = condition)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             scales = "free_y",
             labeller = labeller(CGI_promoter = label_both)) +
  ylim(0, NA) +
  ggtitle("BA9_pos") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.BA9_pos.pdf",
          g,
          base_height = 10)

g <- BA9_pos_rpkm_sx %>%
  ggplot(aes(x = feature, y = mC, col = qrpkm_BA9_pos,
             group = qrpkm_BA9_pos)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             scales = "free_y",
             labeller = labeller(CGI_promoter = label_both)) +
  ylim(0, NA) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("BA9_pos") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.BA9_pos.rpkm.pdf",
          g,
          base_height = 7,
          base_width = 14)

# ------------------------------------------------------------------------------
# NAcc_pos vs. BA9_pos
#

g <- delta_mC_sx %>%
  ggplot(aes(x = feature, y = delta_mC, group = context)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             labeller = labeller(CGI_promoter = label_both)) +
  ggtitle("NAcc_pos vs. BA9_pos") +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.NAcc_pos_vs_BA9_pos.pdf",
          g,
          base_height = 10)

g <- logFC_delta_mC_sx %>%
  ggplot(aes(x = feature, y = delta_mC, col = qLogFC, group = qLogFC)) +
  geom_point() +
  geom_path() +
  facet_grid(context ~ gene_type + CGI_promoter,
             labeller = labeller(CGI_promoter = label_both)) +
  ggtitle("NAcc_pos vs. BA9_pos") +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot("../figures/mC_at_gene_features.NAcc_pos_vs_BA9_pos.logFC.pdf",
          g,
          base_height = 7,
          base_width = 14)
