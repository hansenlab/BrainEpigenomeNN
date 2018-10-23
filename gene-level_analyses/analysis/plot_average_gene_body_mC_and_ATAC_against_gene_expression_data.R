# Plot average gene body mC levels and ATAC counts against gene expression data
# Peter Hickey
# 2018-03-23

### ----------------------------------------------------------------------------
### Setup
###

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(broom)
library(viridis)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
x <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "mC_around_scaled_gene.rds"))
se <- readRDS(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "ATAC_at_gene_level.rds"))

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

# Replace 0 in ATAC RPKM with NA, otherwise these screw up the log ratios
assay(se)[assay(se) == 0] <- NA_real_

# Create the grouped data
grouped_df <- filter(x, !is.na(mC), bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC, DE)) %>%
  group_by(context, gene, gene_type, logFC, DE) %>%
  # NOTE: Aggregate each of CpA (same strand) and CpA (opposite strand)
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, DE, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = as.character(DE)) %>%
  bind_rows(
    data_frame(gene = rownames(se),
           delta_ATAC = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
           context = "ATAC") %>%
  inner_join(genes_df) %>%
  group_by(gene_type, DE, context) %>%
  summarise(cor = list(tidyCor(delta_ATAC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = as.character(DE)))

# Create the partially grouped data
gene_type_grouped_df <- filter(x, !is.na(mC), bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, logFC, gene_type)) %>%
  group_by(context, gene, logFC, gene_type) %>%
  # NOTE: Aggregate each of CpA (same strand) and CpA (opposite strand)
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(gene_type, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(DE = "(all)") %>%
  select(gene_type, DE, everything()) %>%
  bind_rows(
    data_frame(
      gene = rownames(se),
      delta_ATAC = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
                        context = "ATAC") %>%
      inner_join(genes_df) %>%
      group_by(gene_type, context) %>%
      summarise(cor = list(tidyCor(delta_ATAC, logFC))) %>%
      unnest() %>%
      ungroup() %>%
      mutate(DE = "(all)"))
DE_grouped_df <- filter(x, !is.na(mC), bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, logFC, DE)) %>%
  group_by(context, gene, logFC, DE) %>%
  # NOTE: Aggregate each of CpA (same strand) and CpA (opposite strand)
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(DE, context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         DE = as.character(DE)) %>%
  select(gene_type, DE, everything()) %>%
  bind_rows(
    data_frame(
      gene = rownames(se),
      delta_ATAC = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
      context = "ATAC") %>%
      inner_join(genes_df) %>%
      group_by(DE, context) %>%
      summarise(cor = list(tidyCor(delta_ATAC, logFC))) %>%
      unnest() %>%
      ungroup() %>%
      mutate(gene_type = "(all)",
             DE = as.character(DE)))

# Create the ungrouped data
ungrouped_df <- filter(x, !is.na(mC), bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, logFC)) %>%
  group_by(context, gene, logFC) %>%
  # NOTE: Aggregate each of CpA (same strand) and CpA (opposite strand)
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  group_by(context) %>%
  summarise(cor = list(tidyCor(delta_mC, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)",
         DE = "(all)") %>%
  select(gene_type, DE, everything()) %>%
  bind_rows(
    data_frame(
      gene = rownames(se),
      delta_ATAC = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
      context = "ATAC") %>%
      inner_join(genes_df) %>%
      group_by(context) %>%
      summarise(cor = list(tidyCor(delta_ATAC, logFC))) %>%
      unnest() %>%
      ungroup() %>%
      mutate(gene_type = "(all)",
             DE = "(all)"))

# bind grouped, partially grouped, and ungrouped together with factors to
# order the data
sx <- bind_rows(grouped_df, gene_type_grouped_df, DE_grouped_df,
                ungrouped_df) %>%
  mutate(gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")),
         DE = factor(DE,
                     levels = c("FALSE", "TRUE", "(all)")),
         context = factor(context,
                          levels = c("mCA (opposite strand)",
                                     "mCA (same strand)", "mCG", "ATAC")))

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# Contour plots of epigenetic mark against ranked RPKM
#

g <- filter(x, bin >= 200, bin <= 300) %>%
  mutate(context = case_when(
    context == "CpG" ~ "mCG",
    grepl("opposite", context) ~ "mCA (opposite strand)",
    TRUE ~ "mCA (same strand)"),
    context = factor(context,
                     levels = c("mCA (opposite strand)",
                                "mCA (same strand)", "mCG", "ATAC"))) %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_NAcc_pos, rpkm_BA9_pos)) %>%
  mutate(rpkm =
           ifelse(condition == "NAcc_pos", rpkm_NAcc_pos, rpkm_BA9_pos)) %>%
  group_by(condition, context, gene, gene_type, rpkm) %>%
  summarise(val = weighted.mean(mC, nC, na.rm = TRUE))  %>%
  filter(gene_type == "PC") %>%
  group_by(condition, context) %>%
  mutate(ranked_exp = rank(rpkm)) %>%
  bind_rows(
    data_frame(gene = rep(rownames(se), 2),
               condition = rep(c("NAcc_pos", "BA9_pos"), each = nrow(se)),
               val = c(assay(se)[, "NAcc_pos"], assay(se)[, "BA9_pos"]),
               context = "ATAC") %>%
      inner_join(select(genes_df, gene, gene_type, rpkm_NAcc_pos, rpkm_BA9_pos)) %>%
      mutate(rpkm =
               ifelse(condition == "NAcc_pos", rpkm_NAcc_pos, rpkm_BA9_pos)) %>%
      filter(gene_type == "PC") %>%
      group_by(condition, context) %>%
      mutate(ranked_exp = rank(rpkm))
  ) %>%
  ggplot(aes(x = ranked_exp, y = val)) +
  facet_grid(context ~ condition,
             scales = "free_y") +
  geom_density_2d(aes(colour = ..level..)) +
  geom_smooth(se = FALSE, col = "red") +
  xlab("Genes ranked lowest-to-highest RPKM") +
  scale_colour_viridis()
save_plot("../figures/average_gene_body_mC_and_ATAC_vs_ranked_rpkm.PC_only.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos
#

# NOTE: Omitting 1% of genes with largest RPKM
g <- filter(x, bin >= 200, bin <= 300) %>%
  filter(condition == "NAcc_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "mCG",
    grepl("opposite", context) ~ "mCA (opposite strand)",
    TRUE ~ "mCA (same strand)"),
    context = factor(context,
                     levels = c("mCA (opposite strand)",
                                "mCA (same strand)", "mCG", "ATAC"))) %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_NAcc_pos)) %>%
  group_by(condition, context, gene, gene_type, rpkm_NAcc_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ggplot(aes(x = mC, y = rpkm_NAcc_pos)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  ylim(0, quantile(genes_df$rpkm_NAcc_pos, 0.99))
save_plot("../figures/average_gene_body_mC_vs_rpkm.NAcc_pos.pdf",
          g,
          base_height = 10)

# NOTE: Omitting 1% of genes with largest RPKM
g <- data_frame(gene = rownames(se),
                ATAC = assay(se)[, "NAcc_pos"],
                context = "ATAC") %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_NAcc_pos)) %>%
  ggplot(aes(x = ATAC, y = rpkm_NAcc_pos)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  ylim(0, quantile(genes_df$rpkm_NAcc_pos, 0.99))
save_plot("../figures/average_gene_body_ATAC_vs_rpkm.NAcc_pos.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# BA9_pos
#

# NOTE: Omitting 1% of genes with largest RPKM
g <- filter(x, bin >= 200, bin <= 300) %>%
  filter(condition == "BA9_pos") %>%
  mutate(context = case_when(
    context == "CpG" ~ "mCG",
    grepl("opposite", context) ~ "mCA (opposite strand)",
    TRUE ~ "mCA (same strand)"),
    context = factor(context,
                     levels = c("mCA (opposite strand)",
                                "mCA (same strand)", "mCG", "ATAC"))) %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_BA9_pos)) %>%
  group_by(condition, context, gene, gene_type, rpkm_BA9_pos) %>%
  summarise(mC = weighted.mean(mC, nC, na.rm = TRUE)) %>%
  ggplot(aes(x = mC, y = rpkm_BA9_pos)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  ylim(0, quantile(genes_df$rpkm_BA9_pos, 0.99))
save_plot("../figures/average_gene_body_mC_vs_rpkm.BA9_pos.pdf",
          g,
          base_height = 10)

# NOTE: Omitting 1% of genes with largest RPKM
g <- data_frame(gene = rownames(se),
                ATAC = assay(se)[, "BA9_pos"],
                context = "ATAC") %>%
  inner_join(select(genes_df, gene, gene_type, rpkm_BA9_pos)) %>%
  ggplot(aes(x = ATAC, y = rpkm_BA9_pos)) +
  facet_grid(gene_type ~ context,
             scales = "free_x",
             margins = "gene_type") +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  ylim(0, quantile(genes_df$rpkm_BA9_pos, 0.99))
save_plot("../figures/average_gene_body_ATAC_vs_rpkm.BA9_pos.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos vs. BA9_pos
#

g <- filter(x, bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC, DE)) %>%
  group_by(context, gene, gene_type, logFC, DE) %>%
  summarise(delta = mean(delta_mC, na.rm = TRUE)) %>%
  bind_rows(
    data_frame(
      gene = rownames(se),
      delta = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
      context = "ATAC") %>%
      inner_join(select(genes_df, gene, gene_type, logFC, DE))
  ) %>%
  ungroup() %>%
  mutate(context = factor(context,
                          levels = c("mCA (opposite strand)",
                                     "mCA (same strand)", "mCG", "ATAC"))) %>%
  ggplot(aes(x = delta, y = logFC, col = DE)) +
  facet_grid(gene_type ~ context,
             scales = "free",
             margins = "gene_type")  +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = delta, y = logFC),
              method = "lm",
              col = "dodgerBlue",
              se = FALSE,
              lty = 2) +
  geom_smooth(method = "lm",
              se = FALSE) +
  scale_color_manual(values = c("black", "orange"))
save_plot("../figures/average_gene_body_mC_and_ATAC_vs_logFC.pdf",
          g,
          base_height = 10)

g <- filter(x, bin >= 200, bin <= 300) %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos,
         context = case_when(
           context == "CpG" ~ "mCG",
           grepl("opposite", context) ~ "mCA (opposite strand)",
           TRUE ~ "mCA (same strand)")) %>%
  inner_join(select(genes_df, gene, gene_type, logFC, DE)) %>%
  group_by(context, gene, gene_type, logFC, DE) %>%
  summarise(delta = mean(delta_mC, na.rm = TRUE)) %>%
  bind_rows(
    data_frame(
      gene = rownames(se),
      delta = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]),
      context = "ATAC") %>%
      inner_join(select(genes_df, gene, gene_type, logFC, DE))
  ) %>%
  ungroup() %>%
  mutate(context = factor(context,
                          levels = c("mCA (opposite strand)",
                                     "mCA (same strand)", "mCG", "ATAC"))) %>%
  filter(gene_type == "PC") %>%
  ggplot(aes(x = delta, y = logFC, col = DE)) +
  facet_grid(gene_type ~ context,
             scales = "free")  +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = delta, y = logFC),
              method = "lm",
              col = "dodgerBlue",
              se = FALSE,
              lty = 2) +
  geom_smooth(method = "lm",
              se = FALSE) +
  scale_color_manual(values = c("black", "orange"))
save_plot("../figures/average_gene_body_mC_and_ATAC_vs_logFC.PC_only.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# NAcc_pos vs. BA9_pos correlation plot
#

g <- sx %>%
  mutate(feature = "gene_body") %>%
  ggplot(aes(x = feature, y = abs(estimate), col = DE)) +
  facet_grid(gene_type ~ context) +
  geom_point() +
  geom_errorbar(aes(ymax = abs(conf.low), ymin = abs(conf.high))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c(`FALSE` = "black", `TRUE` = "orange", `(all)` = "dodgerBlue"))
save_plot("../figures/cor_average_gene_body_mC_and_ATAC_and_logFC.pdf",
          g,
          base_height = 10)

select(sx, gene_type, DE, context, estimate, conf.low, conf.high) %>%
  arrange(gene_type, DE, context) %>%
  as.data.frame()
# gene_type    DE               context    estimate    conf.low  conf.high
# 1     lncRNA FALSE mCA (opposite strand) -0.24240877 -0.26434535 -0.2202215
# 2     lncRNA FALSE     mCA (same strand) -0.23777919 -0.25976370 -0.2155483
# 3     lncRNA FALSE                   mCG -0.13167580 -0.15463706 -0.1085724
# 4     lncRNA FALSE                  ATAC  0.10658341  0.07697580  0.1360032
# 5     lncRNA  TRUE mCA (opposite strand) -0.55852827 -0.61363083 -0.4980119
# 6     lncRNA  TRUE     mCA (same strand) -0.54818068 -0.60419785 -0.4867719
# 7     lncRNA  TRUE                   mCG -0.39089735 -0.45966902 -0.3174655
# 8     lncRNA  TRUE                  ATAC  0.42262677  0.33596608  0.5021919
# 9     lncRNA (all) mCA (opposite strand) -0.33006991 -0.35003427 -0.3098058
# 10    lncRNA (all)     mCA (same strand) -0.32499809 -0.34503257 -0.3046676
# 11    lncRNA (all)                   mCG -0.18736542 -0.20905354 -0.1654931
# 12    lncRNA (all)                  ATAC  0.17160318  0.14368209  0.1992513
# 13        PC FALSE mCA (opposite strand) -0.26432121 -0.27962907 -0.2488790
# 14        PC FALSE     mCA (same strand) -0.25571540 -0.27109910 -0.2402011
# 15        PC FALSE                   mCG -0.16515830 -0.18118869 -0.1490402
# 16        PC FALSE                  ATAC  0.10126708  0.08411198  0.1183622
# 17        PC  TRUE mCA (opposite strand) -0.63605677 -0.65928949 -0.6116107
# 18        PC  TRUE     mCA (same strand) -0.63045274 -0.65396275 -0.6057262
# 19        PC  TRUE                   mCG -0.51665236 -0.54538351 -0.4867085
# 20        PC  TRUE                  ATAC  0.55004826  0.52070991  0.5780890
# 21        PC (all) mCA (opposite strand) -0.46997882 -0.48179661 -0.4579901
# 22        PC (all)     mCA (same strand) -0.45961306 -0.47157852 -0.4474784
# 23        PC (all)                   mCG -0.32985958 -0.34340241 -0.3161796
# 24        PC (all)                  ATAC  0.25236159  0.23736934  0.2672336
# 25     (all) FALSE mCA (opposite strand) -0.25366479 -0.26626058 -0.2409824
# 26     (all) FALSE     mCA (same strand) -0.24697661 -0.25961742 -0.2342512
# 27     (all) FALSE                   mCG -0.14787703 -0.16106094 -0.1346403
# 28     (all) FALSE                  ATAC  0.09846529  0.08361895  0.1132679
# 29     (all)  TRUE mCA (opposite strand) -0.61156149 -0.63367583 -0.5884485
# 30     (all)  TRUE     mCA (same strand) -0.60501941 -0.62741266 -0.5816263
# 31     (all)  TRUE                   mCG -0.47947319 -0.50680157 -0.4511821
# 32     (all)  TRUE                  ATAC  0.51822841  0.48981687  0.5455404
# 33     (all) (all) mCA (opposite strand) -0.41597772 -0.42638537 -0.4054599
# 34     (all) (all)     mCA (same strand) -0.40723569 -0.41773376 -0.3966289
# 35     (all) (all)                   mCG -0.26846575 -0.28016361 -0.2566882
# 36     (all) (all)                  ATAC  0.22324780  0.20996936  0.2364439
