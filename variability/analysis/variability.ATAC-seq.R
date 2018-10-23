# Explore variability of ATAC-seq counts in peaks
# Peter Hickey
# 2017-07-03

library(SummarizedExperiment)
library(matrixStats)
library(scales)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

extdir <- "../extdata"

### ============================================================================
### ATAC-seq data
###

# ==============================================================================
# Load data
#

conditions <- c("BA9-neg", "BA9-pos", "NA-neg", "NA-pos")
names(conditions) <- conditions
ATAC_SE <-
  readRDS("../../ATAC-seq/objects/flow-sorted-brain-atac.union_narrowPeak_reduced.se.rds")
ATAC_SE <- ATAC_SE[, ATAC_SE$REPLICATE == "rep1"]
colnames(ATAC_SE) <- gsub("NA", "NAcc", colnames(ATAC_SE))
ATAC_SE$TISSUE <- ifelse(ATAC_SE$TISSUE == "NA", "NAcc", ATAC_SE$TISSUE)

# Construct design matrices
design_matrices <-
  lapply(list("Tissue:NeuN" = ~ 0 + TISSUE:NEUN,
              "Donor:NeuN" = ~ 0 + DONOR:NEUN,
              "Donor:Tissue" = ~ 0 + DONOR:TISSUE), function(d) {
                mm <- model.matrix(d, colData(ATAC_SE))
                colnames(mm) <- gsub("TISSUE|NEUN|DONOR", "", colnames(mm))
                mm
              })

# ==============================================================================
# Compute stats
#

computeStats <- function(ATAC_SE, design_matrices, log = FALSE) {
  atac <- cpm(assay(ATAC_SE), log = log)
  lapply(design_matrices, function(design) {
    n <- colSums(design)
    # NOTE: Need at least 2 samples to compute an sd
    design <- design[, n > 1]
    groups <- lapply(seq_len(ncol(design)), function(j) which(design[, j] > 0))
    names(groups) <- colnames(design)
    lapply(groups, function(group) {
      mean <- matrixStats::rowMeans2(atac, cols = group, na.rm = TRUE)
      sd <- matrixStats::rowSds(atac, cols = group, na.rm = TRUE)
      list(mean = mean, sd = sd)
    })
  })
}

stats <- list("sorted" = computeStats(ATAC_SE, design_matrices, log = FALSE))
# NOTE: Re-organise the list for easier access
stats <- map(transpose(stats), flatten)

stats_log <- list("sorted" = computeStats(ATAC_SE, design_matrices, log = TRUE))
# NOTE: Re-organise the list for easier access
stats_log <- map(transpose(stats_log), flatten)

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

stats_log_df <- map2_df(stats_log, names(stats_log), function(x, y) {
  g1_name <- map_chr(strsplit(y, ":"), 1)
  g1 <- map_chr(flatten(map(names(x), ~ strsplit(., ":"))), 1)
  g2_name <- map_chr(strsplit(y, ":"), 2)
  g2 <- map_chr(flatten(map(names(x), ~ strsplit(., ":"))), 2)
  data_frame(!!g1_name := g1,
             !!g2_name := g2,
             stats = map(x, as_data_frame))
})
rm(stats_log)

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
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))

x_log <- filter(stats_log_df,
            is.na(Donor)) %>%
  select(-Donor) %>%
  unnest() %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ NeuN) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Donor.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ NeuN) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Donor.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Donor.2.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Donor.2.pdf",
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
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))

x_log <- filter(stats_log_df, is.na(Tissue)) %>%
  select(-Tissue) %>%
  unnest()  %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ NeuN)  +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Tissue.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ NeuN)  +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Tissue.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Tissue.2.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = NeuN),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Tissue.2.pdf",
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
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))
x_log <- filter(stats_log_df, is.na(NeuN)) %>%
  select(-NeuN) %>%
  unnest()  %>%
  mutate(binned_mean = cut(mean,
                           breaks = quantile(mean, 0:10 / 10),
                           labels = paste0("Q", 1:10),
                           include.lowest = TRUE))

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.NeuN.1.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Tissue),
               outlier.shape = NA) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.NeuN.1.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  ylim(0, 2) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.NeuN.2.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log) +
  geom_boxplot(aes(x = binned_mean, y = sd, colour = Donor),
               outlier.shape = NA) +
  facet_wrap(~ Tissue) +
  ylim(0, 2) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.NeuN.2.pdf",
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

x_log <- filter(stats_log_df,
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
                                 "NAcc" = "chocolate1")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Donor.3.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ NeuN) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Donor.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  facet_grid(Tissue ~ NeuN)  +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Donor.4.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  facet_grid(Tissue ~ NeuN)  +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Donor.4.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over Tissue within Donor:NeuN
#

# NOTE: This is pretty useless here because only 2 samples per combination of
#       Donor:NeuN

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

x_log <- filter(stats_log_df,
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
                                 "pos" = "darkgreen")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Tissue.3.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = NeuN)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Tissue.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x, aes(x = mean, y = sd, colour = NeuN)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(NeuN ~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.Tissue.4.pdf",
       g,
       width = 14,
       height = 14)
g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = NeuN)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("bulk" = "black",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.Tissue.3.pdf",
       g,
       width = 14,
       height = 14)

# ------------------------------------------------------------------------------
# Averaging over NeuN within Donor:Tissue
#

# NOTE: This is pretty useless here because only 2 samples per combination of
#       Donor:Tissue

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

x_log <- filter(stats_log_df,
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
                                 "NAcc" = "chocolate1")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.NeuN.3.pdf",
       g,
       width = 14,
       height = 14)

g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.NeuN.3.pdf",
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
                                 "NAcc" = "chocolate1")) +
  ggtitle("cpm")
ggsave("../figures/ATAC_cpm.sd_vs_mean.NeuN.4.pdf",
       g,
       width = 12,
       height = 12)

g <- ggplot(data = x_log, aes(x = mean, y = sd, colour = Tissue)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(Tissue ~ Donor) +
  scale_colour_manual(values = c("BA24" = "deeppink",
                                 "BA9" = "deepskyblue",
                                 "HC" = "darkgrey",
                                 "NAcc" = "chocolate1")) +
  ggtitle("log2cpm")
ggsave("../figures/ATAC_logcpm.sd_vs_mean.NeuN.4.pdf",
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
