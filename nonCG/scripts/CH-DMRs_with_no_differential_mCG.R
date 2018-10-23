# Find CH-DMRs that do not overlap a CG-DMR and/or those with low CG density
# Peter Hickey
# 2017-11-29

library(SummarizedExperiment)
library(matrixStats)

# ------------------------------------------------------------------------------
# Load data
#

list_of_SEs <- readRDS("../objects/list_of_SE.candidate_CH-DMRs.rds")
load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../Objects/All_Annotated_DMRs_GRanges.rda")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
blocks_pos <- makeGRangesFromDataFrame(sig_block_dmrs)

# ------------------------------------------------------------------------------
# Only retain candidate CH-DMRs with FWER < 0.05
#

# NOTE: CW-DMRs are union of other CH-DMRs with FWER <= 0.05. Since these don't
#       have an FWER attached to each element, we don't apply this filter
k <- match("mCW", names(list_of_SEs))
list_of_SEs[-k] <- lapply(list_of_SEs[-k], function(se) {
  se[rowRanges(se)$fwer <= 50, ]
})

# Add number of CpGs in each CH-DMR
list_of_SEs <- lapply(list_of_SEs, function(se) {
  rowRanges(se)$n_CpGs <- countOverlaps(se, cpgs)
  se
})

# ------------------------------------------------------------------------------
# CH-DMRs with no CG-DMR and/or CG-block
#

sapply(list_of_SEs, function(se) {
  sum(!overlapsAny(se, dmrs_pos)) / length(se)
})
#   mCA (+)   mCT (+)   mCA (-)   mCT (-)       mCW
# 0.7754120 0.7192386 0.7541649 0.7362319 0.7759456
sapply(list_of_SEs, function(se) {
  sum(!overlapsAny(se, blocks_pos)) / length(se)
})
#   mCA (+)   mCT (+)   mCA (-)   mCT (-)       mCW
# 0.7236401 0.5764786 0.6849637 0.5721739 0.7425730
sapply(list_of_SEs, function(se) {
  sum(!overlapsAny(se, union(dmrs_pos, blocks_pos))) / length(se)
})
#   mCA (+)   mCT (+)   mCA (-)   mCT (-)       mCW
# 0.6046376 0.4649898 0.5657839 0.4771014 0.6203377

sapply(list_of_SEs, function(se) {
  sum(width(intersect(rowRanges(se), dmrs_pos)))
})
# mCA (+) mCT (+) mCA (-) mCT (-)     mCW
# 1940966  557111 1571949  610169 2229726
sapply(list_of_SEs, function(se) {
  sum(width(intersect(rowRanges(se), blocks_pos)))
})
#  mCA (+)  mCT (+)  mCA (-)  mCT (-)      mCW
# 11123852  3243285 10109215  3555514 13292174
sapply(list_of_SEs, function(se) {
  sum(width(intersect(rowRanges(se), union(dmrs_pos, blocks_pos))))
})
#  mCA (+)  mCT (+)  mCA (-)  mCT (-)      mCW
# 11880397  3394109 10643504  3714495 14183789

# ------------------------------------------------------------------------------
# Distance from CH-DMR to nearest CG-DMR and/or CG-block for those that don't
# overlap
#

sapply(list_of_SEs, function(se) {
  d <- distanceToNearest(subsetByOverlaps(se, dmrs_pos, invert = TRUE))
  quantile(mcols(d)$distance, 0:10 / 10)
})
#        mCA (+)    mCT (+)  mCA (-)    mCT (-)       mCW
# 0%        57.0       94.0       67      102.0       1.0
# 10%     1600.0     1996.3     1621     1997.7    1525.8
# 20%     2678.4     3735.0     2585     3263.2    2644.0
# 30%     4432.0     6957.7     4091     6013.6    4308.0
# 40%     7969.0    12519.0     7196     9051.0    7690.0
# 50%    13508.0    23431.5    13018    18495.0   13718.0
# 60%    28093.2    48997.6    26734    52608.4   28878.0
# 70%    73616.0   185744.0    76734   179486.4   71399.0
# 80%   202655.0   867606.6   250929   825256.2  194446.0
# 90%   565202.8  3232748.7   749067  2337479.3  497461.2
# 100% 5910784.0 21282489.0 12367187 14959461.0 5125163.0

sapply(list_of_SEs, function(se) {
  d <- distanceToNearest(subsetByOverlaps(se, blocks_pos, invert = TRUE))
  quantile(mcols(d)$distance, 0:10 / 10)
})
        # mCA (+)    mCT (+)   mCA (-)    mCT (-)       mCW
# 0%         57.0      284.0      67.0      102.0       1.0
# 10%      1844.7     2726.0    1704.6     2570.2    1732.0
# 20%      3197.0     5211.2    3035.0     4493.0    2970.0
# 30%      5755.0    11066.0    5412.8     8158.0    5474.8
# 40%     10412.0    21288.8    9903.6    18660.4    9872.0
# 50%     22117.5    78829.0   21779.0    58663.0   20471.5
# 60%     51069.0   242168.2   63110.0   204838.6   45837.4
# 70%    122409.0   716487.8  159962.8   624908.8  107284.0
# 80%    270212.8  1942463.8  384447.0  1594528.4  242888.4
# 90%    589950.4  4033855.0  868705.8  2929321.8  510871.0
# 100% 11350172.0 17145090.0 8354936.0 17198603.0 4773570.0
sapply(list_of_SEs, function(se) {
  d <- distanceToNearest(subsetByOverlaps(se, union(dmrs_pos, blocks_pos),
                                          invert = TRUE))
  quantile(mcols(d)$distance, 0:10 / 10)
})
#         mCA (+)    mCT (+)    mCA (-)    mCT (-)       mCW
# 0%         57.0      284.0       67.0      102.0       1.0
# 10%      1869.5     2729.0     1710.2     2799.0    1806.0
# 20%      3174.0     5334.0     3039.8     4686.0    3008.0
# 30%      5849.0    10624.8     5460.0     8453.2    5513.0
# 40%     10578.0    21108.0     9714.0    18693.0    9989.6
# 50%     21643.0    59728.0    20359.0    58708.0   20198.5
# 60%     50093.0   213388.8    59104.4   209858.4   44267.2
# 70%    136606.5   830619.4   170651.8   802211.4  122283.0
# 80%    322731.0  2585612.4   444599.4  1818191.8  289165.8
# 90%    765532.0  4728041.0  1028662.8  3400554.2  651607.8
# 100% 11350172.0 25066543.0 12367187.0 18016289.0 7401574.0

# ------------------------------------------------------------------------------
# Find 'lonely CH-DMRs'
#

# NOTE: Minimum CpG density in CG-DMRs and CG-blocks
min_cpg_density <- min(c(countOverlaps(dmrs_pos, cpgs) / width(dmrs_pos),
                         countOverlaps(blocks_pos, cpgs) / width(blocks_pos)))

# CpG density for those CH-DMRs that overlap a CG-DMR or CG-block
sapply(list_of_SEs, function(se) {
  se <- subsetByOverlaps(se, union(dmrs_pos, blocks_pos))
  table(countOverlaps(se, cpgs) / width(se) < min_cpg_density)
})
#       mCA (+) mCT (+) mCA (-) mCT (-)  mCW
# FALSE    2641     769    1992     873 2834
# TRUE       70      18      41      29   67

# CpG density for those CH-DMRs that don't overlap a CG-DMR or CG-block
sapply(list_of_SEs, function(se) {
  se <- subsetByOverlaps(se, union(dmrs_pos, blocks_pos), invert = TRUE)
  table(countOverlaps(se, cpgs) / width(se) < min_cpg_density)
})
#       mCA (+) mCT (+) mCA (-) mCT (-)  mCW
# FALSE    3705     625    2416     745 4262
# TRUE      441      59     233      78  478


# CpG density for those CH-DMRs that don't overlap a CG-DMR or CG-block
sapply(list_of_SEs, function(se) {
  se <- subsetByOverlaps(se, union(dmrs_pos, blocks_pos), invert = TRUE)
  quantile(100 * countOverlaps(se, cpgs) / width(se), 0:10 / 10)
})

# (1) CH-DMRs that don't overlap a CG-DMR or CG-block and
# with max(meanDiff) < 0.1 in mCG (S) and max(meanDiff) < 0.1 in mCG (L) but a
# high CpG density
list_of_lonely1_ch_dmrs <- lapply(list_of_SEs, function(se) {
  se <- subsetByOverlaps(se, union(dmrs_pos, blocks_pos), invert = TRUE)
  rr_s <- rowRanges(assay(se, "mCG (S)"))
  rr_l <- rowRanges(assay(se, "mCG (L)"))
  i <- ((rowDiffs(rr_s) < 0.1) | rowAlls(rr_s, value = NA)) &
    ((rowDiffs(rr_s) < 0.1) | rowAlls(rr_l, value = NA)) &
    ((countOverlaps(se, cpgs) / width(se)) >= min_cpg_density)
  se[as.vector(i), ]
})

# (2): CH-DMRs that don't overlap a CG-DMR or CG-block and with low CpG density
list_of_lonely2_ch_dmrs <- lapply(list_of_SEs, function(se) {
  se <- subsetByOverlaps(se, union(dmrs_pos, blocks_pos), invert = TRUE)
  i <- (countOverlaps(se, cpgs) / width(se)) < min_cpg_density
  se[as.vector(i), ]
})

sapply(list_of_lonely1_ch_dmrs, function(x) sum(width(x)) / 10 ^ 6)
#  mCA (+)   mCT (+)   mCA (-)   mCT (-)       mCW
# 9.526109  1.665468  7.902278  1.871662 12.546308

sapply(list_of_lonely2_ch_dmrs, function(x) sum(width(x)) / 10 ^ 6)
# mCA (+)  mCT (+)  mCA (-)  mCT (-)      mCW
# 1.426148 0.209057 0.905199 0.266083 1.634169

saveRDS(list_of_lonely1_ch_dmrs, "../objects/list_of_lonely1_ch_dmrs.rds")
saveRDS(list_of_lonely2_ch_dmrs, "../objects/list_of_lonely2_ch_dmrs.rds")

# ------------------------------------------------------------------------------
# CH-DMRs and gene expression
#

# TODO: Move to its own script
# TODO: list_of_CH_DMRs and list_of_SEs have different nubmer of rows due to
#       different FWER cutoff (<50 vs. <=50)

logFC <- tapply(degs, overlapsAny(degs, list_of_CH_DMRs[[5]]),
                 function(x) x$logFC)
t.test(abs(logFC[["TRUE"]]), abs(logFC[["FALSE"]]))
logFC <- tapply(degs, overlapsAny(degs, dmrs_pos),
                 function(x) x$logFC)
t.test(abs(logFC[["TRUE"]]), abs(logFC[["FALSE"]]))

width <- tapply(degs, overlapsAny(degs, list_of_CH_DMRs[[5]]),
                function(x) log10(width(x)))
t.test(width[["TRUE"]], width[["FALSE"]])
width <- tapply(degs, overlapsAny(degs, dmrs_pos),
                function(x) log10(width(x)))
t.test(width[["TRUE"]], width[["FALSE"]])

AveExpr <- tapply(degs, overlapsAny(degs, list_of_CH_DMRs[[5]]),
                 function(x) x$AveExpr)
t.test(AveExpr[["TRUE"]], AveExpr[["FALSE"]])
AveExpr <- tapply(degs, overlapsAny(degs, dmrs_pos),
                  function(x) x$AveExpr)
t.test(AveExpr[["TRUE"]], AveExpr[["FALSE"]])
