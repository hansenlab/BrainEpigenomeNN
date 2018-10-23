# Plot CH-DMR meanDiffs, stratified by various features
# Peter Hickey
# 2018-02-19


### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(cowplot)

### ----------------------------------------------------------------------------
### Load data
###

list_of_candidate_CH_DMRs <- readRDS(
  "../objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x <- x[x$fwer <= 50]
  x$meanDiff <- x$NA_pos - x$BA9_pos
  x
})

### ----------------------------------------------------------------------------
### Boxplots
###

# - boxplots of meanDiff(mCA+) over pos_CA-DMRs stratified by whether pos_CA-DMR overlaps pos_CT-DMR for each strand

x <- data_frame(meanDiff_mCA_pos = list_of_CH_DMRs[["mCA (+)"]]$meanDiff,
                ol = overlapsAny(list_of_CH_DMRs[["mCA (+)"]],
                                 list_of_CH_DMRs[["mCA (-)"]]))
ggplot(x, aes(x = ol, y = abs(meanDiff_mCA_pos))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  ylim(0, 0.25)


#   - Additionally, stratify by whether pos_CA-DMR overlaps gene+, gene-, both, or none
# - boxplots of meanDiff(mCA+) and nCA+ and nCA- over pos_CA-DMRs stratified by whether pos_CA-DMR overlaps a neg_CA-DMR
#   - Additionally, stratify by whether pos_CA-DMR overlaps gene+, gene-, both, or none
