# Explore variability of WGBS data
# Peter Hickey
# 2017-06-15

library(bsseq)
library(matrixStats)
library(scales)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

extdir <- "../extdata"

### ============================================================================
### mCG (small smooth)
###

# ==============================================================================
# Load data
#

mCG_small_sorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
# NOTE: Tidy up sample colData
colnames(mCG_small_sorted) <- gsub("NA", "NAcc", colnames(mCG_small_sorted))
mCG_small_sorted$Tissue <- gsub("NA", "NAcc", mCG_small_sorted$Tissue)
mCG_small_sorted$NeuN <- factor(mCG_small_sorted$NeuN,
                                levels = c("bulk", "neg", "pos"))
# Construct design matrices
sorted_design_matrices <-
  lapply(list("Tissue:NeuN" = ~ 0 + Tissue:NeuN,
              "Donor:NeuN" = ~ 0 + Individual:NeuN,
              "Donor:Tissue" = ~ 0 + Individual:Tissue), function(d) {
                mm <- model.matrix(d, colData(mCG_small_sorted))
                colnames(mm) <- gsub("Tissue|NeuN|Individual", "", colnames(mm))
                mm
              })

mCG_small_bulk <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))
# NOTE: Tidy up sample colData
colnames(mCG_small_bulk) <- gsub("NA", "NAcc", colnames(mCG_small_bulk))
mCG_small_bulk$Tissue <- gsub("NA", "NAcc", mCG_small_bulk$Tissue)
mCG_small_bulk$NeuN <- factor(rep("bulk", ncol(mCG_small_bulk)),
                              levels = c("bulk", "neg", "pos"))
# NOTE: Not interested in caudate so removing
mCG_small_bulk <- mCG_small_bulk[, !grepl("caudate", colnames(mCG_small_bulk))]

# Construct design matrices
bulk_design_matrices <-
  lapply(list("Tissue:NeuN" = ~ 0 + Tissue:NeuN,
              "Donor:NeuN" = ~ 0 + Individual:NeuN,
              "Donor:Tissue" = ~ 0 + Individual:Tissue), function(d) {
                mm <- model.matrix(d, colData(mCG_small_bulk))
                colnames(mm) <- gsub("Tissue|NeuN|Individual", "", colnames(mm))
                mm
              })

# ==============================================================================
# Compute stats
#

computeStats <- function(BSseq, design_matrices) {
  meth <- as.array(getMeth(BSseq))
  lapply(design_matrices, function(design) {
    n <- colSums(design)
    # NOTE: Need at least 2 samples to compute an sd
    design <- design[, n > 1]
    groups <- lapply(seq_len(ncol(design)), function(j) which(design[, j] > 0))
    names(groups) <- colnames(design)
    lapply(groups, function(group) {
      mean <- matrixStats::rowMeans2(meth, cols = group, na.rm = TRUE)
      sd <- matrixStats::rowSds(meth, cols = group, na.rm = TRUE)
      list(mean = mean, sd = sd)
    })
  })
}

stats <- list("sorted" = computeStats(mCG_small_sorted, sorted_design_matrices),
              "bulk" = computeStats(mCG_small_bulk[], bulk_design_matrices))
# NOTE: Re-organise the list for easier access
stats <- map(transpose(stats), flatten)

stats_df <- map2_df(stats, names(stats), function(x, y) {
  g1_name <- map_chr(strsplit(y, ":"), 1)
  g1 <- map_chr(flatten(map(names(x), ~ strsplit(., ":"))), 1)
  g2_name <- map_chr(strsplit(y, ":"), 2)
  g2 <- map_chr(flatten(map(names(x), ~ strsplit(., ":"))), 2)
  data_frame(!!g1_name := g1,
             !!g2_name := g2,
             stats = map(x, as_data_frame))
})
rm(stats)

# ==============================================================================
# Boxplots of sd vs. binned mean
#

# ------------------------------------------------------------------------------
# Averaging over Donor within Tissue:NeuN
#

x <- filter(stats_df,
            is.na(Donor)) %>%
  select(-Donor) %>%
  unnest() %>%
  mutate(binned_mean = cut(mean,
                           breaks = 0:10 / 10,
                           include.lowest = TRUE))
g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ NeuN) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.Donor.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.Donor.2.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over Tissue within Donor:NeuN
#

x <- filter(stats_df, is.na(Tissue)) %>%
  select(-Tissue) %>%
  unnest()  %>%
  mutate(binned_mean = cut(mean,
                           breaks = 0:10 / 10,
                           include.lowest = TRUE))

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ NeuN) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.Tissue.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.Tissue.2.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over NeuN within Donor:Tissue
#

x <- filter(stats_df, is.na(NeuN)) %>%
  select(-NeuN) %>%
  unnest()  %>%
  mutate(binned_mean = cut(mean,
                           breaks = 0:10 / 10,
                           include.lowest = TRUE))
g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.NeuN.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  ylim(0, 0.3)
ggsave("../figures/mCG_small.sd_vs_mean.NeuN.2.pdf",
       g,
       width = 14,
       height = 14)

# ==============================================================================
# Scatterplot of mean(sd) vs. mean(mean) where each point is 1% of CpGs
# (inspired by Boyle [2017])
#

# ------------------------------------------------------------------------------
# Averaging over Donor within Tissue:NeuN
#

x <- filter(stats_df,
            is.na(Donor)) %>%
  select(-Donor) %>%
  unnest() %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:100 / 100),
                           include.lowest = TRUE)) %>%
  group_by(Tissue, NeuN, binned_mean) %>%
  summarise(mean = mean(mean),
            sd = mean(sd))

g <- ggplot(data = x, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ NeuN) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1"))
ggsave("../figures/mCG_small.sd_vs_mean.Donor.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  facet_grid(Tissue ~ NeuN)  +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1"))
ggsave("../figures/mCG_small.sd_vs_mean.Donor.4.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over Tissue within Donor:NeuN
#

x <- filter(stats_df,
            is.na(Tissue)) %>%
  select(-Tissue) %>%
  unnest() %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:100 / 100),
                           include.lowest = TRUE)) %>%
  group_by(Donor, NeuN, binned_mean) %>%
  summarise(mean = mean(mean),
            sd = mean(sd))

g <- ggplot(data = x, aes(x = mean, y = sd, colour = NeuN)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen"))
ggsave("../figures/mCG_small.sd_vs_mean.Tissue.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x, aes(x = mean, y = sd, colour = NeuN)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(NeuN ~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen"))
ggsave("../figures/mCG_small.sd_vs_mean.Tissue.4.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over NeuN within Donor:Tissue
#

x <- filter(stats_df,
            is.na(NeuN)) %>%
  select(-NeuN) %>%
  unnest() %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:100 / 100),
                           include.lowest = TRUE)) %>%
  group_by(Donor, Tissue, binned_mean) %>%
  summarise(mean = mean(mean),
            sd = mean(sd))

g <- ggplot(data = x, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1"))
ggsave("../figures/mCG_small.sd_vs_mean.NeuN.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(Tissue ~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1"))
ggsave("../figures/mCG_small.sd_vs_mean.NeuN.4.pdf",
       g,
       width = 12,
       height = 12)

# ==============================================================================
# Density plots of mean and sd
#

# TODO

# ==============================================================================
# ECDF plots of mean and sd
#

# TODO

# ==============================================================================
# TODO: t-stat of NA vs BA9 in NeuN+ and NeuN-.
#       Compare distribution of t-stat in NeuN+ to NeuN-
#       Compare meanDiff in NeuN+ to NeuN-
#
