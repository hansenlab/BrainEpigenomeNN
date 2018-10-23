# Plot some example CH-DMRs
# Peter Hickey
# 2017-12-19

library(bsseq)
library(matrixStats)
library(limma)
library(scales)
library(dplyr)

### ============================================================================
### Load data
###

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "neg_CA", "pos_CT", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
                                 ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq[, BSseq$NeuN == "pos"]
})
# NOTE: Fix up colours
list_of_CH_BSseq <- lapply(list_of_CH_BSseq, function(BSseq) {
  BSseq$col <- BSseq$Tissue_color
  BSseq
})

# ------------------------------------------------------------------------------
# mCG
#

CG_small_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.small.sorted.somatic.all"))
CG_small_BSseq <- CG_small_BSseq[, CG_small_BSseq$NeuN == "pos"]

CG_large_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                  "BS.fit.large.sorted.somatic.all"))
CG_large_BSseq <- CG_large_BSseq[, CG_large_BSseq$NeuN == "pos"]

list_of_CG_BSseq <- list("mCG (S)" = CG_small_BSseq,
                         "mCG (L)" = CG_large_BSseq)
# NOTE: Fix up colours
list_of_CG_BSseq <- lapply(list_of_CG_BSseq, function(BSseq) {
  BSseq$col <- BSseq$Tissue_color
  BSseq
})

# ------------------------------------------------------------------------------
# mC
#

list_of_BSseq <- c(list_of_CG_BSseq, list_of_CH_BSseq)

# ------------------------------------------------------------------------------
# DMRs and blocks
#

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.rds")
# NOTE: CG-DMRs and CG-blocks have used FWER <= 0.05 rather than FWER < 0.05
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x[x$fwer <= 50]
})
DMRs_and_blocks <- c(list("mCG (S)" = dmrs_pos,
                          "mCG (L)" = makeGRangesFromDataFrame(sig_block_dmrs)),
                     list_of_CH_DMRs)

# ------------------------------------------------------------------------------
# Gene models for plots
#

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds")

### ============================================================================
### Plots
###

source("../analysis/EDA_functions.R")

# ------------------------------------------------------------------------------
# DEG-focused (i.e. DEGs with CH-DMRs)
#

# Find DEGs where at least 50% of gene is covered by a CA_pos-DMR
tmp <- filter(rna_atac_meth, DE, db == "PC") %>%
  arrange(desc(abs(expLogFC))) %>%
  pull(gene)
regions <- subsetByOverlaps(unstrand(gencode_features$genes[tmp]),
                            DMRs_and_blocks[["mCA (+)"]])
# Get percentage of overlap (pol), i.e. how much of the gene is covered by a
# CA_pos-DMR
pol <- sapply(regions, function(r) {
  sum(width(
    GenomicRanges::intersect(r, unstrand(DMRs_and_blocks[["mCA (+)"]])))) /
    width(r)
})

pdf("../figures/DEGs_with_majority_covered_by_CA_pos-DMRs.pdf")
plotManyRegions2(
  list_of_CG_BSseq = list_of_CG_BSseq,
  list_of_CH_BSseq = list_of_CH_BSseq,
  regions = regions[pol > 0.5],
  list_of_CG_addRegions = DMRs_and_blocks[c("mCG (S)", "mCG (L)")],
  list_of_CH_addRegions = DMRs_and_blocks[names(list_of_CH_BSseq)],
  ATAC_SE = NULL,
  geneTrack = RefSeq_exons3,
  extend = 10000)
dev.off()
