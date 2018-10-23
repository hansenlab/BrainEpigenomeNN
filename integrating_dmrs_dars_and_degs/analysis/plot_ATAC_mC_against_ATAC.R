# Plot mC against ATAC in bins along genome
# Peter Hickey
# 2018-03-13

### ----------------------------------------------------------------------------
### Setup
###

library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(viridis)

### ----------------------------------------------------------------------------
### Load data
###

x <- readRDS("../objects/mC_in_bins.rds")
se <- readRDS("../objects/ATAC_in_bins.rds")

### ----------------------------------------------------------------------------
### Data wrangling
###

assay(se) <- pmax(assay(se), 0.01)

sx <- x %>%
  select(-nC) %>%
  spread(condition, mC) %>%
  mutate(delta_mC = NAcc_pos - BA9_pos) %>%
  group_by(context, bin) %>%
  # NOTE: Aggregate each of CpA (same strand) and CpA (opposite strand)
  summarise(delta_mC = mean(delta_mC, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(
    data_frame(
      bin = as.integer(rownames(se)),
      delta_ATAC = log2(assay(se)[, "NAcc_pos"] / assay(se)[, "BA9_pos"]))) %>%
  mutate(delta_ATAC = ifelse(is.infinite(delta_ATAC), NA, delta_ATAC))

# The aimed for binsize should be the mode of the distribution
binsize <- names(which.max(table(x$width)))

### ----------------------------------------------------------------------------
### Plots
###

# ------------------------------------------------------------------------------
# Contour plots of ATAC against mC
#

# NOTE: x-axis is truncated to exclude bins with largest 0.5% ATAC
xmax <- quantile(assay(se), 0.995)
g <- inner_join(
  x,
  data_frame(
    condition = rep(colnames(se), nrow(se)),
    bin = rep(as.integer(rownames(se)), ncol(se)),
    ATAC = as.vector(assay(se)))) %>%
  ggplot(aes(x = ATAC, y = mC)) +
  facet_grid(context ~ condition) +
  geom_density_2d(aes(colour = ..level..),
                  bins = 100) +
  scale_colour_viridis() +
  scale_x_continuous(limits = c(0, xmax)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("ATAC-seq RPKM") +
  ylab("mC") +
  ggtitle(paste0("binsize = ", binsize, " bp"))
save_plot("../figures/mC_against_ATAC_in_bins.contours.pdf",
          g,
          base_height = 10)

g <- inner_join(
  x,
  data_frame(
    condition = rep(colnames(se), nrow(se)),
    bin = rep(as.integer(rownames(se)), ncol(se)),
    ATAC = as.vector(assay(se)))) %>%
  ggplot(aes(x = ATAC, y = mC)) +
  facet_grid(context ~ condition) +
  stat_density_2d(aes(fill = ..density..),
                  geom = "raster",
                  contour = FALSE) +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(0, xmax)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("ATAC-seq RPKM") +
  ylab("mC") +
  ggtitle(paste0("binsize = ", binsize, " bp"))
save_plot("../figures/mC_against_ATAC_in_bins.filled.pdf",
          g,
          base_height = 10)

# ------------------------------------------------------------------------------
# Contour plots of delta-ATAC against delta-mC
#

# TODO: Don't remove bins with infinite logFC
# NOTE: xlim and ylim chosen by looking at data
g <- ggplot(sx, aes(x = delta_ATAC, y = delta_mC)) +
  facet_grid(context ~ .,
             scales = "free") +
  stat_density_2d(aes(fill = ..density..),
                  geom = "raster",
                  contour = FALSE) +
  scale_fill_viridis() +
  geom_smooth(col = "red", method = "lm", se = FALSE) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  xlab("ATAC-seq logFC") +
  ylab("mC meanDiff") +
  ggtitle(paste0("binsize = ", binsize, " bp"))
save_plot("~/blonko.pdf",
          g,
          base_height = 10)

### ----------------------------------------------------------------------------
### Tables
###

# ------------------------------------------------------------------------------
# Number of bins with mC and ATAC measurements
#

inner_join(
  x,
  data_frame(
    condition = rep(colnames(se), nrow(se)),
    bin = rep(as.integer(rownames(se)), ncol(se)),
    ATAC = as.vector(assay(se)))) %>%
  group_by(condition, context) %>%
  count(is.na(mC), ATAC == 0) %>%
  as.data.frame()
#    condition context is.na(mC) ATAC == 0      n
# 1    BA9_pos mCA (-)     FALSE     FALSE 267407
# 2    BA9_pos mCA (-)     FALSE      TRUE    467
# 3    BA9_pos mCA (-)      TRUE     FALSE    326
# 4    BA9_pos mCA (-)      TRUE      TRUE  19913
# 5    BA9_pos mCA (+)     FALSE     FALSE 267407
# 6    BA9_pos mCA (+)     FALSE      TRUE    491
# 7    BA9_pos mCA (+)      TRUE     FALSE    326
# 8    BA9_pos mCA (+)      TRUE      TRUE  19889
# 9    BA9_pos     mCG     FALSE     FALSE 267243
# 10   BA9_pos     mCG     FALSE      TRUE    378
# 11   BA9_pos     mCG      TRUE     FALSE    490
# 12   BA9_pos     mCG      TRUE      TRUE  20002
# 13  NAcc_pos mCA (-)     FALSE     FALSE 267388
# 14  NAcc_pos mCA (-)     FALSE      TRUE    469
# 15  NAcc_pos mCA (-)      TRUE     FALSE    353
# 16  NAcc_pos mCA (-)      TRUE      TRUE  19903
# 17  NAcc_pos mCA (+)     FALSE     FALSE 267402
# 18  NAcc_pos mCA (+)     FALSE      TRUE    481
# 19  NAcc_pos mCA (+)      TRUE     FALSE    339
# 20  NAcc_pos mCA (+)      TRUE      TRUE  19891
# 21  NAcc_pos     mCG     FALSE     FALSE 267244
# 22  NAcc_pos     mCG     FALSE      TRUE    377
# 23  NAcc_pos     mCG      TRUE     FALSE    497
# 24  NAcc_pos     mCG      TRUE      TRUE  19995

# ------------------------------------------------------------------------------
# mCH desserts
#

filter(x, context == "mCA (+)") %>%
  pull(mC) %>%
  S4Vectors::Rle()

