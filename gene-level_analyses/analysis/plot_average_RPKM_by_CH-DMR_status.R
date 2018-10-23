# Plot average expression level of genes with and without CA-DMRs
# Peter Hickey
# 2018-03-08

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
  file.path("flow-sorted-brain-gene-level_analyses", "objects",
            "gene-level_analyses_data.rda"))

### ----------------------------------------------------------------------------
### Data wrangling
###

CA_DMRs <- c(list_of_CH_DMRs[["mCA (+)"]], list_of_CH_DMRs[["mCA (-)"]])

x <- genes_df %>%
  inner_join(data_frame(gene = genes_df$gene,
                        `CH-DMR` = overlapsAny(genes[genes_df$gene], CA_DMRs)))

count(x, DE, `CH-DMR`)
# # A tibble: 4 x 3
#   DE    `CH-DMR`     n
#   <lgl> <lgl>    <int>
# 1 FALSE FALSE    29197
# 2 FALSE TRUE      1202
# 3 TRUE  FALSE     2263
# 4 TRUE  TRUE       689

### ----------------------------------------------------------------------------
### Plots
###

# NOTE: Omitting 10% of genes with largest RPKM from plot
g <- ggplot(x, aes(x = `CH-DMR`, y = ave_rpkm)) +
  facet_grid(DE ~ .,
             margins = "DE",
             labeller = labeller(DE = label_both)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, quantile(genes_df$ave_rpkm, 0.90))) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "orange"))
save_plot("../figures/average_RPKM_by_CH-DMR_status.pdf",
          g,
          base_height = 10)

### ----------------------------------------------------------------------------
### Hypothesis test
###

t.test(ave_rpkm ~ `CH-DMR`, x, alternative = "greater")
# Welch Two Sample t-test
#
# data:  ave_rpkm by CH-DMR
# t = 3.3948, df = 31587, p-value = 0.0003438
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   1.105706      Inf
# sample estimates:
#   mean in group FALSE  mean in group TRUE
# 2.6102990           0.4652318
t.test(ave_rpkm ~ `CH-DMR`, filter(x, gene_type == "PC", !is.na(logFC)),
       alternative = "greater")
# Welch Two Sample t-test
#
# data:  ave_rpkm by CH-DMR
# t = 18.493, df = 6627.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   1.140411      Inf
# sample estimates:
#   mean in group FALSE  mean in group TRUE
# 1.8299151           0.5781508
t.test(ave_rpkm ~ `CH-DMR`, filter(x, gene_type == "PC", DE, !is.na(logFC)),
       alternative = "greater")

wilcox.test(ave_rpkm ~ `CH-DMR`, x, alternative = "greater")
#
# Wilcoxon rank sum test with continuity correction
#
# data:  ave_rpkm by CH-DMR
# W = 35164000, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0

wilcox.test(ave_rpkm ~ `CH-DMR`, filter(x, gene_type == "PC", !is.na(logFC)),
            alternative = "greater")
# Wilcoxon rank sum test with continuity correction
#
# data:  ave_rpkm by CH-DMR
# W = 12396000, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
wilcox.test(ave_rpkm ~ `CH-DMR`,
            filter(x, gene_type == "PC", DE, !is.na(logFC)),
            alternative = "greater")
#
# Wilcoxon rank sum test with continuity correction
#
# data:  ave_rpkm by CH-DMR
# W = 650420, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
