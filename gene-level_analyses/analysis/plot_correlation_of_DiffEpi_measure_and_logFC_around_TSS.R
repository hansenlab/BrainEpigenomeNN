# Plot correlation of DiffEpi measure with gene expression logFC around TSS of
# genes
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
  file.path(
    extdir, "objects", "DiffEpi_measure_around_TSS.rds")) %>%
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

# Create the grouped data
grouped_df <- inner_join(
  filter(x, !is.na(DiffEpi_measure)),
  select(
    genes_df, gene, logFC, gene_type)) %>%
  mutate(DiffEpi = factor(case_when(
    grepl("opposite", DiffEpi) ~ "CA-DMRs (opposite strand)",
    grepl("CA_", DiffEpi) ~ "CA-DMRs",
    TRUE ~ DiffEpi),
    levels = c("DARs", "bigDARs",
               "CG-DMRs", "CG-blocks",
               "CA-DMRs", "CA-DMRs (opposite strand)"))) %>%
  group_by(gene_type, DiffEpi, distance) %>%
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
  group_by(DiffEpi, distance) %>%
  summarise(cor = list(tidyCor(DiffEpi_measure, logFC))) %>%
  unnest() %>%
  ungroup() %>%
  mutate(gene_type = "(all)") %>%
  select(gene_type, everything())

# bind grouped and ungrouped together with factors to order the data
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

g <- ggplot(sx, aes(x = distance, y = estimate, col = DiffEpi)) +
  facet_grid(gene_type ~ DiffEpi) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = -0.5, lty = 3) +
  geom_hline(yintercept = 0.5, lty = 3) +
  geom_smooth(aes(x = distance, y = conf.low),
              method = "loess", span = 0.1, se = FALSE, size = 0.3) +
  geom_smooth(aes(x = distance, y = estimate),
              method = "loess", span = 0.1, se = FALSE) +
  geom_smooth(aes(x = distance, y = conf.high),
              method = "loess", span = 0.1, se = FALSE, size = 0.3) +
  ylim(-1, 1) +
  ylab("Correlation with gene expression") +
  scale_colour_manual(values = cols)
save_plot(
  "../figures/cor_DiffEpi_and_logFC_around_TSS.faceted_with_CI.pdf",
  g,
  base_height = 10)

# NOTE: Exclude 'CA-DMRs (opposite strand)' because these are basically the same as
#       'CA-DMRs'
g <- ggplot(filter(sx, DiffEpi != "CA-DMRs (opposite strand)"),
            aes(x = distance, y = abs(estimate), col = DiffEpi)) +
  facet_grid(gene_type ~ .) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE) +
  ylim(0, 1) +
  scale_colour_manual(values = cols)
save_plot(
  "../figures/absolute_cor_DiffEpi_and_logFC_around_TSS.overlaid_without_CI.pdf",
  g,
  base_height = 10)
