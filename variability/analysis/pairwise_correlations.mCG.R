# Look at pairwise correlation of mCG between samples
# Peter Hickey
# 2017-06-19

### ============================================================================
### Setup
###

library(bsseq)
library(dplyr)
library(gtools)
library(ggplot2)
library(tidyr)
library(purrr)
library(coop)
library(gplots)
library(RColorBrewer)

extdir <- "../extdata"

mCG_small_sorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
colnames(mCG_small_sorted) <- gsub("NA", "NAcc", colnames(mCG_small_sorted))
mCG_small_sorted$Tissue <- gsub("NA", "NAcc", mCG_small_sorted$Tissue)
mCG_small_bulk <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))
# NOTE: Not interested in caudate so removing
mCG_small_bulk <- mCG_small_bulk[, !grepl("caudate", colnames(mCG_small_bulk))]
colnames(mCG_small_bulk) <- gsub("NA", "NAcc", colnames(mCG_small_bulk))
mCG_small_bulk$Tissue <- gsub("NA", "NAcc", mCG_small_bulk$Tissue)
# NOTE: Add bulk as NeuN type
colnames(mCG_small_bulk) <- paste0(colnames(mCG_small_bulk), "_bulk")
mCG_small_bulk$NeuN <- "bulk"

# ------------------------------------------------------------------------------
# Construct SummarizedExperiment with methylation data for common CpGs
#

ol <- findOverlaps(mCG_small_sorted, mCG_small_bulk)
meth_sorted <- as.matrix(getMeth(mCG_small_sorted))[queryHits(ol), ]
meth_bulk <- as.matrix(getMeth(mCG_small_bulk))[subjectHits(ol), ]
combined <- SummarizedExperiment(
  assays = SimpleList(meth = cbind(meth_sorted, meth_bulk)),
  rowRanges = rowRanges(mCG_small_sorted)[queryHits(ol)],
  colData = combine(as.data.frame(colData(mCG_small_sorted)),
                    as.data.frame(colData(mCG_small_bulk))))
rm(meth_sorted, meth_bulk)
meth <- assay(combined, "meth")

# ------------------------------------------------------------------------------
# Construct SummarizedExperiment with methylation data for common CpGs with
# average methylation in 1 kb bins
#

bins <- tileGenome(seqlengths(combined),
                   tilewidth = 1000,
                   cut.last.tile.in.chrom = TRUE)
ol <- findOverlaps(bins, combined)


binned_meth <- do.call(rbind, lapply(as.list(ol), function(i) {
  colMeans2(assay(combined, "meth"), rows = i)
}))
binned_meth <- binned_meth[!rowAlls(binned_meth, value = NaN), ]
binned_meth_row_centered <- binned_meth - rowMeans2(binned_meth)
binned_meth_row_centered_cp <- crossprod(binned_meth_row_centered)
binned_meth_row_centered_cp_svd <- svd(binned_meth_row_centered_cp)

### ============================================================================
### Compute correlations
###

# ------------------------------------------------------------------------------
# Using all common CpGs
#

pcors <- pcor(meth)
x <- pcors
x[upper.tri(x)] <- NA_real_
pcors_df <- data_frame(cor = as.vector(x),
                       S1 = rep(rownames(x), times = ncol(x)),
                       S2 = rep(colnames(x), each = nrow(x))) %>%
  # NOTE: Remove diagonal and upper triangle
  filter(S1 != S2,
         !is.na(cor)) %>%
  extract(data = .,
          col = S1,
          into = c("Donor1", "Tissue1", "NeuN1"),
          regex = "([[:alnum:]]+)\\_([[:alnum:]]+)\\_([[:alnum:]]+)") %>%
  extract(data = .,
          col = S2,
          into = c("Donor2", "Tissue2", "NeuN2"),
          regex = "([[:alnum:]]+)\\_([[:alnum:]]+)\\_([[:alnum:]]+)") %>%
  mutate(NeuN1 = ifelse(is.na(NeuN1), "bulk", NeuN1),
          NeuN2 = ifelse(is.na(NeuN2), "bulk", NeuN2))

# ------------------------------------------------------------------------------
# Using average methylation in 1kb bins
#

pcors_binned <- pcor(binned_meth)
x <- pcors_binned
x[upper.tri(x)] <- NA_real_
pcors_binned_df <- data_frame(cor = as.vector(x),
                              S1 = rep(rownames(x), times = ncol(x)),
                              S2 = rep(colnames(x), each = nrow(x))) %>%
  # NOTE: Remove diagonal and upper triangle
  filter(S1 != S2,
         !is.na(cor)) %>%
  extract(data = .,
          col = S1,
          into = c("Donor1", "Tissue1", "NeuN1"),
          regex = "([[:alnum:]]+)\\_([[:alnum:]]+)\\_([[:alnum:]]+)") %>%
  extract(data = .,
          col = S2,
          into = c("Donor2", "Tissue2", "NeuN2"),
          regex = "([[:alnum:]]+)\\_([[:alnum:]]+)\\_([[:alnum:]]+)") %>%
  mutate(NeuN1 = ifelse(is.na(NeuN1), "bulk", NeuN1),
         NeuN2 = ifelse(is.na(NeuN2), "bulk", NeuN2))

### ============================================================================
### Plot correlations
###

# ------------------------------------------------------------------------------
# Heatmap of correlation matrix
#

labels <- factor(gsub("^[0-9]+\\_", "", rownames(pcors)),
                 levels = c("BA9_neg", "BA9_pos", "BA9_bulk",
                            "BA24_neg", "BA24_pos", "BA24_bulk",
                            "HC_neg", "HC_pos", "HC_bulk",
                            "NAcc_neg", "NAcc_pos", "NAcc_bulk"))
row_side_colors <- rep(c("deeppink", "deepskyblue", "darkgrey", "chocolate1"),
                       each = 3)
names(row_side_colors) <- levels(labels)
col_side_colors <- rep(c("purple", "darkgreen", "black"), times = 4)
names(col_side_colors) <- levels(labels)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf("../figures/common_CpGs.correlation_heatmap.pdf")
heatmap.2(pcors,
          trace = "none",
          ColSideColors = col_side_colors[labels],
          RowSideColors = row_side_colors[labels],
          col = hmcol,
          cexRow = 0.1 + 1 / log10(ncol(pcors)),
          cexCol = 0.1 + 1 / log10(ncol(pcors)),
          main = "Common CpGs")
heatmap.2(pcors_binned,
          trace = "none",
          ColSideColors = col_side_colors[labels],
          RowSideColors = row_side_colors[labels],
          col = hmcol,
          cexRow = 0.1 + 1 / log10(ncol(pcors)),
          cexCol = 0.1 + 1 / log10(ncol(pcors)),
          main = "Common CpGs (1kb bins)")
dev.off()

# ------------------------------------------------------------------------------
# Intercelltype variation
#

# UP TO HERE: Still need to make final plots

x <- filter(pcors_df,
            Donor1 == Donor2,
            Tissue1 == Tissue2) %>%
  mutate(Donor = Donor1,
         Tissue = Tissue1) %>%
  group_by(i = row_number()) %>%
  mutate(pair = paste(sort(c(NeuN1, NeuN2)), collapse = ":")) %>%
  ungroup() %>%
  select(-i)
g <- ggplot(x) +
  geom_boxplot(aes(x = pair, y = cor)) +
  facet_grid(~ Tissue) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../figures/common_CpGs.intercelltype.boxplot.pdf", g)

g <- ggplot(x, aes(x = Donor, y = cor, col = pair)) +
  geom_line(aes(group = pair)) +
  geom_point() +
  facet_grid(~ Tissue) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("../figures/common_CpGs.levels_of_Donor:Tissue.pdf", g)
ggsave("../figures/common_CpGs.levels_of_Donor:Tissue.pdf", g)

TukeyHSD(aov(cor ~ pair, filter(x, Tissue == "BA24")))
TukeyHSD(aov(cor ~ pair, filter(x, Tissue == "BA9")))
TukeyHSD(aov(cor ~ pair, filter(x, Tissue == "HC")))
TukeyHSD(aov(cor ~ pair, filter(x, Tissue == "NAcc")))














x <- filter(cors, design == "Tissue:NeuN") %>%
  mutate(pair = case_when(
    grepl("5085", .$S1) & grepl("5086", .$S2) ~ "5085:5086",
    grepl("5085", .$S1) & grepl("5250", .$S2) ~ "5085:5250",
    grepl("5085", .$S1) & grepl("5284", .$S2) ~ "5085:5284",
    grepl("5085", .$S1) & grepl("5358", .$S2) ~ "5085:5358",
    grepl("5085", .$S1) & grepl("5404", .$S2) ~ "5085:5405",
    grepl("5085", .$S1) & grepl("5628", .$S2) ~ "5085:5628",
    grepl("5086", .$S1) & grepl("5248", .$S2) ~ "5086:5248",
    grepl("5086", .$S1) & grepl("5250", .$S2) ~ "5086:5250",
    grepl("5086", .$S1) & grepl("5284", .$S2) ~ "5086:5284",
    grepl("5086", .$S1) & grepl("5358", .$S2) ~ "5086:5358",
    grepl("5086", .$S1) & grepl("5404", .$S2) ~ "5086:5404",
    grepl("5086", .$S1) & grepl("5540", .$S2) ~ "5086:5540",
    grepl("5086", .$S1) & grepl("5552", .$S2) ~ "5086:5552",
    grepl("5086", .$S1) & grepl("5569", .$S2) ~ "5086:5569",
    grepl("5086", .$S1) & grepl("5570", .$S2) ~ "5086:5570",
    grepl("5086", .$S1) & grepl("5628", .$S2) ~ "5086:5628",
    grepl("5248", .$S1) & grepl("5250", .$S2) ~ "5248:5250",
    grepl("5248", .$S1) & grepl("5284", .$S2) ~ "5248:5284",
    grepl("5248", .$S1) & grepl("5358", .$S2) ~ "5248:5358",
    grepl("5248", .$S1) & grepl("5540", .$S2) ~ "5248:5540",
    grepl("5248", .$S1) & grepl("5552", .$S2) ~ "5248:5552",
    grepl("5248", .$S1) & grepl("5569", .$S2) ~ "5248:5569",
    grepl("5248", .$S1) & grepl("5570", .$S2) ~ "5248:5570",
    grepl("5248", .$S1) & grepl("5628", .$S2) ~ "5248:5628",
    grepl("5250", .$S1) & grepl("5284", .$S2) ~ "5250:5284",
    grepl("5250", .$S1) & grepl("5358", .$S2) ~ "5250:5358",
    grepl("5250", .$S1) & grepl("5404", .$S2) ~ "5250:5404",
    grepl("5250", .$S1) & grepl("5540", .$S2) ~ "5250:5540",
    grepl("5250", .$S1) & grepl("5552", .$S2) ~ "5250:5552",
    grepl("5250", .$S1) & grepl("5569", .$S2) ~ "5250:5569",
    grepl("5250", .$S1) & grepl("5570", .$S2) ~ "5250:5570",
    grepl("5250", .$S1) & grepl("5628", .$S2) ~ "5250:5628",
    grepl("5284", .$S1) & grepl("5358", .$S2) ~ "5284:5358",
    grepl("5284", .$S1) & grepl("5404", .$S2) ~ "5284:5404",
    grepl("5284", .$S1) & grepl("5552", .$S2) ~ "5284:5552",
    grepl("5284", .$S1) & grepl("5569", .$S2) ~ "5284:5569",
    grepl("5284", .$S1) & grepl("5570", .$S2) ~ "5284:5570",
    grepl("5284", .$S1) & grepl("5628", .$S2) ~ "5284:5628",
    grepl("5358", .$S1) & grepl("5404", .$S2) ~ "5358:5404",
    grepl("5358", .$S1) & grepl("5540", .$S2) ~ "5358:5540",
    grepl("5358", .$S1) & grepl("5552", .$S2) ~ "5358:5552",
    grepl("5358", .$S1) & grepl("5569", .$S2) ~ "5358:5569",
    grepl("5358", .$S1) & grepl("5570", .$S2) ~ "5358:5570",
    grepl("5358", .$S1) & grepl("5628", .$S2) ~ "5358:5628",
    grepl("5404", .$S1) & grepl("5552", .$S2) ~ "5404:5552",
    grepl("5404", .$S1) & grepl("5570", .$S2) ~ "5404:5570",
    grepl("5404", .$S1) & grepl("5628", .$S2) ~ "5404:5628",
    grepl("5540", .$S1) & grepl("5552", .$S2) ~ "5540:5552",
    grepl("5540", .$S1) & grepl("5569", .$S2) ~ "5540:5569",
    grepl("5540", .$S1) & grepl("5570", .$S2) ~ "5540:5570",
    grepl("5540", .$S1) & grepl("5628", .$S2) ~ "5540:5628",
    grepl("5552", .$S1) & grepl("5569", .$S2) ~ "5552:5569",
    grepl("5552", .$S1) & grepl("5570", .$S2) ~ "5552:5570",
    grepl("5552", .$S1) & grepl("5628", .$S2) ~ "5552:5628",
    grepl("5569", .$S1) & grepl("5570", .$S2) ~ "5569:5570",
    grepl("5569", .$S1) & grepl("5628", .$S2) ~ "5569:5628",
    grepl("5570", .$S1) & grepl("5628", .$S2) ~ "5570:5628"))

g <- x %>%
  group_by(pair) %>%
  add_tally() %>%
  # NOTE: Only retain pairs of donors with at least 3 Tissue:NeuN pairs
  filter(n >= 3) %>%
  ggplot() +
  geom_boxplot(aes(x = pair, y = cor)) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../figures/common_CpGs.Tissue:NeuN.pdf", g)

g <- x %>%
  group_by(pair) %>%
  add_tally() %>%
  # NOTE: Only retain pairs of donors with at least 3 Tissue:NeuN pairs
  filter(n >= 3) %>%
  ggplot() +
  geom_point(aes(x = level, y = cor)) +
  geom_line(aes(x = level, y = cor, col = pair, group = pair)) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../figures/common_CpGs.levels_of_Tissue:NeuN.pdf", g)

g <- x %>%
  extract(.,
          col = level,
          into = c("Tissue", "NeuN"),
                   regex = "([[:alnum:]]+)\\:([[:alnum:]]+)",
          remove = FALSE) %>%
  ggplot(.) +
  geom_boxplot(aes(x = NeuN, y = cor, fill = NeuN)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("bulk" = "white",
                                 "neg" = "purple",
                                 "pos" = "darkgreen")) +
  ggtitle("Common CpGs") +
  facet_grid(~ Tissue)
ggsave("../figures/common_CpGs.levels_of_Tissue:NeuN.boxplots.pdf", g)

# Interindividual correlation vs.
g <- filter(cors,
       design == "Donor:Tissue",
       grepl("neg|pos", S1) & grepl("neg|pos", S2)) %>%
  extract(.,
          col = level,
          into = c("Donor", "Tissue"),
          regex = "([[:alnum:]]+)\\:([[:alnum:]]+)",
          remove = FALSE) %>%
  mutate(group = "Donor matched",
         NeuN = "pos:neg") %>%
  select(Tissue, cor, group, NeuN) %>%
  bind_rows(filter(cors,
                   design == "Tissue:NeuN",
                   grepl("neg|pos", S1) & grepl("neg|pos", S2)) %>%
              extract(.,
                      col = level,
                      into = c("Tissue", "NeuN"),
                      regex = "([[:alnum:]]+)\\:([[:alnum:]]+)",
                      remove = FALSE) %>%
              mutate(group = "NeuN matched") %>%
              select(Tissue, cor, group, NeuN)) %>%
  ggplot(aes(x = NeuN, y = cor, colour = group)) +
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(~ Tissue) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../figures/common_CpGs.interindividual_vs_intercelltype.correlation.pdf", g)

# NOTE: Can run an ANOVA testing for differences between NeuN within Tissue
# y <- x %>%
#   extract(.,
#           col = level,
#           into = c("Tissue", "NeuN"),
#           regex = "([[:alnum:]]+)\\:([[:alnum:]]+)",
#           remove = FALSE)
# TukeyHSD(aov(cor ~ Tissue*NeuN, y))
# TukeyHSD(aov(cor ~ Tissue*NeuN, filter(y, Tissue == "BA24"))

x <- filter(cors, design == "Donor:NeuN") %>%
  mutate(pair = case_when(
    (grepl("BA24", .$S1) & grepl("BA9", .$S2)) |
      (grepl("BA24", .$S2) & grepl("BA9", .$S1)) ~ "BA24:BA9",
    (grepl("BA24", .$S1) & grepl("HC", .$S2)) |
      (grepl("BA24", .$S2) & grepl("HC", .$S1)) ~ "BA24:HC",
    (grepl("BA24", .$S1) & grepl("NA", .$S2)) |
      (grepl("BA24", .$S2) & grepl("NA", .$S1)) ~ "BA24:NAcc",
    (grepl("BA9", .$S1) & grepl("HC", .$S2)) |
      (grepl("BA9", .$S2) & grepl("HC", .$S1)) ~ "BA9:HC",
    (grepl("BA9", .$S1) & grepl("NA", .$S2)) |
      (grepl("BA9", .$S2) & grepl("NA", .$S1)) ~ "BA9:NAcc",
    (grepl("HC", .$S1) & grepl("NA", .$S2)) |
      (grepl("HC", .$S2) & grepl("NA", .$S1))~ "HC:NAcc"))

g <- ggplot(x) +
  geom_boxplot(aes(x = pair, y = cor)) +
  ggtitle("Common CpGs")
ggsave("../figures/common_CpGs.Donor:NeuN.pdf", g)

g <- x %>%
  ggplot(.) +
  geom_point(aes(x = level, y = cor, col = pair)) +
  geom_line(aes(x = level, y = cor, col = pair, group = pair)) +
  ggtitle("Common CpGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../figures/common_CpGs.levels_of_Donor:NeuN.pdf", g)


### ============================================================================
### Using CpGs in 'promoters' (8 kb upstream, 2kb downstream of first exon) of
### protein coding genes
###

# TODO?

# ------------------------------------------------------------------------------
# Get co-ordinates of 'promoters'
#
