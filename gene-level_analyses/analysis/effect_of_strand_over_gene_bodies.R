# Investigate effect of strand on mCA and mCT over gene bodies
# Peter Hickey
# 2018-03-25

# TODO: Same plot using 1 kb or 10 kb bins along genome

library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

se <- readRDS(
  file.path(extdir, "flow-sorted-brain-gene-level_analyses",  "objects",
            "mC_at_gene_level.same_strand_and_opposite_strand.rds"))

### ----------------------------------------------------------------------------
### Data wrangling
###

x <- data_frame(
  Individual = rep(se$Individual, each = nrow(se)),
  mC_ss = as.vector(assay(se, "mCA (same strand)")),
  mC_os = as.vector(assay(se, "mCA (opposite strand)")),
  Tissue = rep(se$Tissue, each = nrow(se)))

x %>%
  group_by(Individual, Tissue) %>%
  summarise(cor = cor(mC_ss, mC_os, use = "complete.obs")) %>%
  arrange(cor)

### ----------------------------------------------------------------------------
### Plots
###

g <- ggplot(aes(x = mC_ss, y = mC_os), data = x) +
  facet_grid(Individual ~ Tissue) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_density_2d(aes(colour = ..level..)) +
  geom_smooth(se = FALSE, col = "red", method = "lm") +
save_plot("~/kraken.pdf", g, base_height = 10)

# TODO: Look at missingness patterns
g <- ggplot(aes(x = mC_ss, y = mC_os),
            data = x) +
  facet_grid(Individual ~ Tissue) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_miss_point()
save_plot("~/slipslop.pdf", g, base_height = 10)


g <- ggplot(aes(x = mC_ss - mC_os), data = x) +
  facet_grid(Individual ~ Tissue) +
  geom_density()
save_plot("~/blonko.pdf", g, base_height = 10)
