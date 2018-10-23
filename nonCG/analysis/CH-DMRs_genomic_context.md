CH-DMRs genomic context
================
Peter Hickey
29 September 2017

-   [Size of CH-DMRs](#size-of-ch-dmrs)
-   [EDA of overlap between sets of CH-DMRs, CG-blocks, CG-DMRs, OCRs, and DARs](#eda-of-overlap-between-sets-of-ch-dmrs-cg-blocks-cg-dmrs-ocrs-and-dars)
-   [TODOs](#todos)

``` r
library(GenomicRanges)
library(dplyr)
library(purrr)


load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../../FlowSortingProject/Objects/All_BLOCK_POS_DMRs_fwer50.rda")
CG_DMRs <- dmrs_pos
CG_blocks <- makeGRangesFromDataFrame(sig_block_dmrs)
list_of_candidate_CH_DMRs <-
  readRDS("../objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$successful_permutations <= 0.05, ]
})
```

Size of CH-DMRs
===============

``` r
lengths(list_of_candidate_CH_DMRs)
#> mCA (+) mCA (-) mCT (+) mCT (-) 
#>   24424   24771   11040   11387
lengths(list_of_CH_DMRs)
#> mCA (+) mCA (-) mCT (+) mCT (-) 
#>    6931    4859    1499    1740

# Proportion of candidate CH-DMRs that pass FWER <= 0.05 cutoff
100 * lengths(list_of_CH_DMRs) / lengths(list_of_candidate_CH_DMRs)
#>  mCA (+)  mCA (-)  mCT (+)  mCT (-) 
#> 28.37783 19.61568 13.57790 15.28058

# Size of CH-DMRs in Mb
sapply(list_of_CH_DMRs, function(xx) sum(width(xx))) / 10 ^ 6
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 31.573596 26.534304  6.890655  7.532227
```

EDA of overlap between sets of CH-DMRs, CG-blocks, CG-DMRs, OCRs, and DARs
==========================================================================

``` r


# Proportion of CH-DMRs that overlap CH-DMRs in other context/strand
# Read as '% colname overlapping rowname'
sapply(list_of_CH_DMRs, function(x) {
  sapply(list_of_CH_DMRs, function(y) {
    sum(overlapsAny(x, y)) / length(x)
  })
})
#>           mCA (+)   mCA (-)   mCT (+)   mCT (-)
#> mCA (+) 1.0000000 0.7721753 0.9679787 0.8977011
#> mCA (-) 0.5525898 1.0000000 0.8725817 0.9005747
#> mCT (+) 0.1868417 0.2389381 1.0000000 0.4591954
#> mCT (-) 0.2017025 0.2792756 0.5256838 1.0000000

# Proportion of CH-DMRs that overlap CH-DMRs (in bp) in other context/strand
sapply(list_of_CH_DMRs, function(x) {
  sapply(list_of_CH_DMRs, function(y) {
    sum(width(GenomicRanges::intersect(x, y))) / sum(width(x))
  })
})
#>           mCA (+)   mCA (-)   mCT (+)   mCT (-)
#> mCA (+) 1.0000000 0.7190613 0.9459538 0.8740157
#> mCA (-) 0.6042958 1.0000000 0.8635879 0.9059403
#> mCT (+) 0.2064459 0.2242639 1.0000000 0.4540900
#> mCT (-) 0.2085060 0.2571670 0.4963692 1.0000000

# Proportion of CH-DMRs that overlap CG-DMRs
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, CG_DMRs)) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.2239215 0.2412019 0.2808539 0.2632184

# Proportion of CH-DMRs that overlap CG-DMRs (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(GenomicRanges::intersect(x, CG_DMRs))) / sum(width(x))
})
#>    mCA (+)    mCA (-)    mCT (+)    mCT (-) 
#> 0.06178736 0.05993860 0.08176015 0.08137354

# Proportion of CH-DMRs that overlap CG-blocks
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, CG_blocks)) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.2761506 0.3107635 0.4209473 0.4252874

# Proportion of CH-DMRs that overlap CG-blocks (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(GenomicRanges::intersect(x, CG_blocks))) / sum(width(x))
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.3537497 0.3846144 0.4738515 0.4727696

# Proportion of CH-DMRs that overlap DARs
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, dars_pos)) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.4870870 0.5143034 0.5283522 0.5132184

# Proportion of CH-DMRs that overlap DARs (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(GenomicRanges::intersect(x, dars_pos))) / sum(width(x))
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.1988993 0.1964957 0.2381520 0.2356922

# Proportion of CH-DMRs that overlap OCRs
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, ocrs_overall)) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.9300245 0.9438156 0.9499666 0.9505747

# Proportion of CH-DMRs that overlap OCRs (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(GenomicRanges::intersect(x, ocrs_overall))) / sum(width(x))
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.4329867 0.4302852 0.4822573 0.4740727

# Proportion of CH-DMRs that overlap CG-DMRs or CG-blocks
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, GenomicRanges::union(CG_DMRs, CG_blocks))) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.3947482 0.4282774 0.5336891 0.5212644

# Proportion of CH-DMRs that overlap CG-DMRs or CG-blocks (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(
    GenomicRanges::intersect(x, GenomicRanges::union(CG_DMRs, CG_blocks)))) /
    sum(width(x))
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.3778934 0.4051878 0.4961666 0.4942421

# Proportion of CH-DMRs that overlap CG-DMRs or CG-blocks or DARs
sapply(list_of_CH_DMRs, function(x) {
  sum(overlapsAny(x, GenomicRanges::union(
    GenomicRanges::union(CG_DMRs, CG_blocks), dars_pos))) / length(x)
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.6303564 0.6620704 0.7158105 0.6982759

# Proportion of CH-DMRs that overlap CG-DMRs or CG-blocks or DARs (in bp)
sapply(list_of_CH_DMRs, function(x) {
  sum(width(
    GenomicRanges::intersect(x,
                             GenomicRanges::union(
                               GenomicRanges::union(CG_DMRs, CG_blocks),
                               dars_pos)))) /
    sum(width(x))
})
#>   mCA (+)   mCA (-)   mCT (+)   mCT (-) 
#> 0.4756289 0.4944935 0.5860427 0.5785179

# Distance distribution for CH-DRMs from nearest CG-DMR
sapply(list_of_CH_DMRs, function(x) {
  quantile(mcols(distanceToNearest(x, CG_DMRs))$distance, 0:10 / 10)
})
#>      mCA (+)   mCA (-)   mCT (+)   mCT (-)
#> 0%         0       0.0       0.0       0.0
#> 10%        0       0.0       0.0       0.0
#> 20%        0       0.0       0.0       0.0
#> 30%     5309    4284.4    1201.2    1394.6
#> 40%    14583   12313.2    6437.4    5617.2
#> 50%    28932   25636.0   14202.0   13686.5
#> 60%    52354   46951.2   25656.0   26263.4
#> 70%    89986   80787.4   49100.0   50165.7
#> 80%   160328  143381.6   96949.2   93337.2
#> 90%   318471  286155.4  204905.0  200552.4
#> 100% 4361354 2024940.0 1765493.0 1789665.0

# Distance distribution for CH-DRMs from nearest CG-DMR or CG-block
sapply(list_of_CH_DMRs, function(x) {
  quantile(
    mcols(distanceToNearest(
      x, GenomicRanges::union(CG_DMRs, CG_blocks)))$distance,
    0:10 / 10)
})
#>      mCA (+)   mCA (-)   mCT (+)   mCT (-)
#> 0%         0       0.0       0.0       0.0
#> 10%        0       0.0       0.0       0.0
#> 20%        0       0.0       0.0       0.0
#> 30%        0       0.0       0.0       0.0
#> 40%      858       0.0       0.0       0.0
#> 50%    13371    8584.0       0.0       0.0
#> 60%    34519   26664.8    7442.6    6777.0
#> 70%    70794   59979.4   24445.4   26139.3
#> 80%   134080  117110.4   74949.8   70593.2
#> 90%   291788  264966.2  186534.0  177487.5
#> 100% 4361354 2024940.0 1765493.0 1789665.0
```

TODOs
=====

-   \[ \] mCH is reportedly enriched in regions of low CG density (<https://www.ncbi.nlm.nih.gov/pubmed/24362762>). Are CH-DMRs enriched in regions of low CG density?
