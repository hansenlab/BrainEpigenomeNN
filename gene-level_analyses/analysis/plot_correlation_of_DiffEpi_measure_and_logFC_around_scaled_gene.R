# Plot correlation of DiffEpi measure with gene expression logFC around genes
# using equal number of bins along gene body +/- 2 gene body equivalents (GBE)
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
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
# NOTE: Not using bigDARs
x <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "DiffEpi_measure_around_scaled_gene.rds")) %>%
  filter(DiffEpi != "bigDARs")

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

# NOTE: Have to make the '(all)' facet/group levels by hand
# NOTE: Stratifying by 'overlaps_another_gene' did not reveal any striking
#       pattern (data not shown)

# Create the grouped data
grouped_df <- inner_join(
  filter(x, !is.na(DiffEpi_measure)),
  select(genes_df, gene, logFC, gene_type)) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  group_by(gene_type, DiffEpi, bin) %>%
  summarise(cor = list(tidyCor(DiffEpi_measure, logFC))) %>%
  unnest() %>%
  ungroup()

# Create the ungrouped data
ungrouped_df <- inner_join(filter(x, !is.na(DiffEpi_measure)),
                           select(genes_df, gene, logFC)) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  group_by(DiffEpi, bin) %>%
  summarise(cor = list(tidyCor(DiffEpi_measure, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)") %>%
  select(gene_type, everything())

# Bind grouped and ungrouped together with factors to order the data
sx <- bind_rows(grouped_df, ungrouped_df) %>%
  mutate(gene_type = factor(gene_type,
                            levels = c("lncRNA", "PC", "(all)")))

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
# Full scale
#

g <- ggplot(sx,
            aes(x = bin, y = estimate, col = DiffEpi)) +
  facet_grid(gene_type ~ DiffEpi) +
  geom_vline(xintercept = 2000, lty = 2) +
  geom_vline(xintercept = 2100, lty = 2) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = -0.5, lty = 3) +
  geom_hline(yintercept = 0.5, lty = 3) +
  geom_smooth(aes(x = bin, y = conf.low),
              method = "loess", span = 0.1, se = FALSE, size = 0.3) +
  geom_smooth(aes(x = bin, y = estimate),
              method = "loess", span = 0.1, se = FALSE) +
  geom_smooth(aes(x = bin, y = conf.high),
              method = "loess", span = 0.1, se = FALSE, size = 0.3) +
  ylim(-1, 1) +
  ylab("Correlation with gene expression") +
  scale_colour_manual(values = cols)
save_plot(
  "../figures/cor_DiffEpi_and_logFC_around_scaled_gene.faceted_with_CI.pdf",
  g,
  base_height = 10)

# NOTE: Exclude 'CA-DMRs (opposite strand)' because these are basically the same as
#       'CA-DMRs'
g <- ggplot(filter(sx, DiffEpi != "CA-DMRs (opposite strand)"),
            aes(x = bin, y = abs(estimate), col = DiffEpi)) +
  facet_grid(gene_type ~ .) +
  geom_vline(xintercept = 2000, lty = 2) +
  geom_vline(xintercept = 2100, lty = 2) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.5, lty = 3) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  ylim(0, 1) +
  scale_colour_manual(values = cols)
save_plot(
  "../figures/absolute_cor_DiffEpi_and_logFC_around_scaled_gene.overlaid_without_CI.pdf",
  g,
  base_height = 10)

# ------------------------------------------------------------------------------
# Zoom
#

# u = Number of bases upstream/downstream in
#     plot_correlation_of_DiffEpi_measure_upstream_and_downstream_of_gene.R
u <- 100000
w <- median(genes_df$width)
# n = number of bins to pad around scaled gene body. Derived by solving
#     w / (u + w) = w / ((2 * n + 1) * w) for n
n <- u / (2 * w)
# nbgtile = Number of tiles used across gene body
ngbtile <- 100

# NOTE: Exclude 'CA-DMRs (opposite strand)' because these are basically the same as
#       'CA-DMRs'
g <- ggplot(filter(sx, DiffEpi != "CA-DMRs (opposite strand)"),
            aes(x = bin, y = abs(estimate), col = DiffEpi)) +
  facet_grid(gene_type ~ .) +
  geom_vline(xintercept = 2000, lty = 2) +
  geom_vline(xintercept = 2100, lty = 2) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.5, lty = 3) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  ylim(0, 1) +
  scale_colour_manual(values = cols) +
  coord_cartesian(xlim = c(2000 - ngbtile * n, 2100 + ngbtile * n))
save_plot(
  "../figures/absolute_cor_DiffEpi_and_logFC_around_scaled_gene.overlaid_without_CI.zoom.pdf",
  g,
  base_height = 10)
