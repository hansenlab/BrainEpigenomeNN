# Plot WGBS and ATAC-seq data for DEGs without various epigenetic marks
# Peter Hickey
# 2017-12-07

library(bsseq)

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

load("../objects/assays-and-features.rda")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds")

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
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

load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
list_of_candidate_CH_DMRs <- readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.rds")
# NOTE: CG-DMRs and CG-blocks have used FWER <= 0.05 rather than FWER < 0.05
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x[x$fwer <= 50]
})
DMRs_and_blocks <- c(list("mCG (S)" = dmrs_pos,
                          "mCG (L)" = makeGRangesFromDataFrame(sig_block_dmrs)),
                     list_of_CH_DMRs)

# ------------------------------------------------------------------------------
# plotManyRegions2()
# NOTE: The DEGs to plot were created in mCG_mCG_ATAC_RNA.R
#

library(scales)
source("../../nonCG/analysis/EDA_functions.R")

DEGs_with_CH_DMRs_but_no_CG_DMR_or_DAR <-
  readRDS("../tmp/DEGs_with_CH-DMRs_but_no_CG-DMR_or_DAR.rds")
regions <- unstrand(
  gencode_features$genes[DEGs_with_CH_DMRs_but_no_CG_DMR_or_DAR])

pdf("../figures/DEGs_with_CH-DMRs_but_no_CG-DMR_or_DAR.pdf")
plotManyRegions2(
  list_of_CG_BSseq = list_of_CG_BSseq,
  list_of_CH_BSseq = list_of_CH_BSseq,
  regions = regions,
  list_of_CG_addRegions = DMRs_and_blocks[c("mCG (S)", "mCG(L)")],
  list_of_CH_addRegions = DMRs_and_blocks[names(list_of_CH_BSseq)],
  ATAC_SE = NULL,
  geneTrack = RefSeq_exons3,
  extend = 10000,
  main = DEGs_with_CH_DMRs_but_no_CG_DMR_or_DAR)
dev.off()

DEGs_without_CG_DMR_or_DAR_or_CH_DMR <-
  readRDS("../tmp/DEGs_without_CG-DMR_or_DAR_or_CH-DMR.rds")
regions <- unstrand(
  gencode_features$genes[DEGs_without_CG_DMR_or_DAR_or_CH_DMR])

pdf("../figures/DEGs_without_CG-DMR_or_DAR_or_CH-DMR.pdf")
plotManyRegions2(
  list_of_CG_BSseq = list_of_CG_BSseq,
  list_of_CH_BSseq = list_of_CH_BSseq,
  regions = regions,
  list_of_CG_addRegions = DMRs_and_blocks[c("mCG (S)", "mCG(L)")],
  list_of_CH_addRegions = DMRs_and_blocks[names(list_of_CH_BSseq)],
  ATAC_SE = NULL,
  geneTrack = RefSeq_exons3,
  extend = 10000,
  main = DEGs_without_CG_DMR_or_DAR_or_CH_DMR)
dev.off()
