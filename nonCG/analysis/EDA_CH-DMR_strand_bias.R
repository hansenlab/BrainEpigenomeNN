# EDA of relationship between CH-DMRs over genes and the strand of the
# underlying gene
# 2018-02-07
# Peter Hickey

library(bsseq)
library(dplyr)
library(DelayedMatrixStats)
library(ggplot2)
library(limma)
library(scales)

### ============================================================================
### Load data
###

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$successful_permutations <= 0.05, ]
})
strand(list_of_CH_DMRs[["mCA (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCA (-)"]]) <- "-"
strand(list_of_CH_DMRs[["mCT (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCT (-)"]]) <- "-"

g <- gencode_features$genes[rna_atac_meth$gene]


table(overlapsAny(list_of_CH_DMRs[["mCA (+)"]], g))
# FALSE  TRUE
#  4482  2449
table(overlapsAny(invertStrand(list_of_CH_DMRs[["mCA (+)"]]), g))
# FALSE  TRUE
#  4183  2748
fisher.test(
  table(overlapsAny(list_of_CH_DMRs[["mCA (+)"]], g),
        overlapsAny(invertStrand(list_of_CH_DMRs[["mCA (+)"]]), g)))

table(overlapsAny(list_of_CH_DMRs[["mCA (-)"]], g))
# FALSE  TRUE
#  3107  1752
table(overlapsAny(invertStrand(list_of_CH_DMRs[["mCA (-)"]]), g))
# FALSE  TRUE
#  2742  2117
fisher.test(
  table(overlapsAny(list_of_CH_DMRs[["mCA (-)"]], g),
        overlapsAny(invertStrand(list_of_CH_DMRs[["mCA (-)"]]), g)))

# Look for CH-DMRs on the opposite strand to a single gene that don't
# overlap (or are near?) a CH-DMR on the other strand

i <- overlapsAny(invertStrand(list_of_CH_DMRs[["mCA (+)"]]), g)
z <- subsetByOverlaps(granges(list_of_CH_DMRs[["mCA (+)"]][i]),
                      list_of_CH_DMRs[["mCA (-)"]],
                      ignore.strand = TRUE,
                      invert = TRUE,
                      maxgap = 100000)

z <- resize(z, width(z) + 10 * width(z), fix = "center")
as.character(z[4])
