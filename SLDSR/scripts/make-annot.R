# Make .annot file for use with ldsc containing brain-derived genomic features
# Peter Hickey
# 2018-01-30

### ============================================================================
### NOTEs
###

# - Adapted from tutorial at https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#building-on-top-of-the-finucane-et-al-baseline-model

### ============================================================================
### Setup
###

library(readr)
library(GenomicRanges)
library(R.utils)
library(rtracklayer)
library(AnnotationHub)

options("mc.cores" = 40)

seqlevels <- 1:22

### ============================================================================
### Load and prepare data
###

# ------------------------------------------------------------------------------
# Our data
#

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
load("../../Objects/All_Annotated_DMRs_GRanges.rda")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
load("../../Objects/Annotated_POSvNEG_BLOCKs_GRanges.rda")

# Non-differential features
OCRs <- list(
  "NAcc_pos_OCRs" = ocrs_NAcc_pos,
  "NAcc_neg_OCRs" = ocrs_NAcc_neg,
  "BA9_pos_OCRs" = ocrs_BA9_pos,
  "BA9_neg_OCRs" = ocrs_BA9_neg,
  "overall_OCRs" = ocrs_overall)
non_DARs <- subsetByOverlaps(ocrs_overall,
                             c(dars_pos, dars_neg, dars_pos_vs_neg),
                             invert = TRUE)
set.seed(3000)
OCR_subsets <- mclapply(
  list("POS_non-DARs" = dars_pos, "POSvsNEG_non-DARs" = dars_pos_vs_neg),
  function(dars) {
    # Keep sampling until get non-DARs with about the same total width as DARs
    tmp <- sample(non_DARs, length(dars))
    while (sum(width(tmp)) < sum(width(dars))) {
      message(sum(width(dars)) - sum(width(tmp)))
      tmp <- c(tmp, sample(
        subsetByOverlaps(non_DARs, tmp, invert = TRUE, type = "equal"),
        1000))
    }
    list("sameN" = sample(non_DARs, length(dars)),
         "sameWidth" = tmp)
  }, mc.cores = getOption("mc.cores"))
our_non_differential_features <- c(OCRs, unlist(OCR_subsets, recursive = FALSE))
# NOTE: Rename to avoid '.' in names
names(our_non_differential_features) <-
  gsub("\\.", "_", names(our_non_differential_features))

# Differential features
POS_CG_DMRs <- Annotated_POS_DMRs_gr
POSvsNEG_CG_DMRs <- Annotated_POSvNEG_DMRs_gr
POS_CG_blocks <- makeGRangesFromDataFrame(sig_block_dmrs,
                                          keep.extra.columns = TRUE)
POSvsNEG_CG_blocks <- Annotated_POSvNEG_BLOCKs_gr
POS_DARs <- dars_pos
POS_bigDARs <- dars_pos[abs(dars_pos$logFC) >= 1]
POSvsNEG_DARs <- dars_pos_vs_neg
POSvsNEG_bigDARs <- POSvsNEG_DARs[abs(POSvsNEG_DARs$logFC) >= 1]
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- endoapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$successful_permutations <= 0.05]
})
CH_DMRs <- reduce(unlist(GRangesList(list_of_CH_DMRs)), ignore.strand = TRUE)

our_differential_features <- list(
  `POS_CG-DMRs` = POS_CG_DMRs,
  `POSvsNEG_CG-DMRs` = POSvsNEG_CG_DMRs,
  `POS_CG-blocks` = POS_CG_blocks,
  `POSvsNEG_CG-blocks` = POSvsNEG_CG_blocks,
  POS_DARs = POS_DARs,
  POS_bigDARs = POS_bigDARs,
  POSvsNEG_DARs = POSvsNEG_DARs,
  POSvsNEG_bigDARs = POSvsNEG_bigDARs,
  `CH-DMRs` = CH_DMRs,
  `POS_CG-DMRs_union_POS_bigDARs` = union(POS_CG_DMRs, POS_bigDARs)
)

# ------------------------------------------------------------------------------
# Other people's data
#

Luo_small_CG_DMRs <- readRDS("../../Objects/Luo_small_CG_DMRs.rds")
load("../../Objects/Lister_hs_CGDMRs.rda")

# Non-differential features
brain_enhancers <- unlinked_enhancers["H3K27ac"]

# Subset of chromHMM track with regulatory region-like states
ah <- AnnotationHub()
# Ten brain samples with chromHMM from RoadMap; identified using
# query(ah, c("chromHMM", "brain"))
list_of_chromHMM <- ah[c("AH46920", "AH46921", "AH46922", "AH46923", "AH46924",
                         "AH46925", "AH46926", "AH46927", "AH46934", "AH46935")]
chromHMM_regulatory_regions <- lapply(list_of_chromHMM, function(x) {
  x <- x[[1]]
  x[x$name %in%
      c("Bivalent Enhancer", "Bivalent/Poised TSS", "Genic enhancers",
        "Flanking Active TSS", "Active TSS", "Strong transcription",
        "Enhancers")]
})

chromHMM_regulatory_regions <- list(
  "chromHMM_union" = reduce(unlist(GRangesList(chromHMM_regulatory_regions))))
others_non_differential_features <- c(brain_enhancers,
                                      chromHMM_regulatory_regions)

# Differential features
Lister_CG_DMRs <- makeGRangesFromDataFrame(Lister_DMRs,
                                           keep.extra.columns = TRUE)
Lister_CG_DMRs <-
  sort(sortSeqlevels(keepSeqlevels(Lister_CG_DMRs, paste0("chr", seqlevels),
                                   pruning.mode = "coarse")))
# NOTE: Lister's CG-DMRs form a motley crew of 4 sets (from the Supplementary
#       Material):
#       - NeuN+ hyper-mCG: This set includes CG-DMRs where mCG/CG is larger in
#         NeuN+ compared with NeuN– and Fetal; larger than Fetal only; or
#         larger than NeuN– only
#       - NeuN– hyper-mCG: This set includes CG-DMRs where mCG/CG is larger in
#         NeuN– compared with NeuN+, and those in which mCG/CG is larger in
#         NeuN– than in fetal.
#       - NeuN+ hypo-mCG: This set includes CG-DMRs where mCG/CG is smaller in
#         NeuN+ compared with NeuN– and Fetal; or smaller compared with Fetal
#         only.
#       - NeuN– hypo-mCG: This set includes CG-DMRs where mCG/CG is smaller in
#         NeuN– compared with NeuN+ and Fetal, or Fetal only.
#       Notably, CG-DMRs may overlap one another, so we merge to form disjoint
#       regions
Lister_CG_DMRs <- disjoin(Lister_CG_DMRs)
set.seed(7000)
CG_DMR_subsets <- mcmapply(
  function(others_dmrs, our_dmrs) {
    # Keep sampling until get others DMRs with about the same total width as
    # our DMRs
    tmp <- sample(others_dmrs, length(our_dmrs))
    while (sum(width(tmp)) < sum(width(our_dmrs))) {
      message(sum(width(our_dmrs)) - sum(width(tmp)))
      tmp <- c(tmp, sample(
        subsetByOverlaps(others_dmrs, tmp, invert = TRUE, type = "equal"),
        1000))
    }
    list("sameN" = sample(others_dmrs, length(our_dmrs)),
         "sameWidth" = tmp)
  }, others_dmrs = list("Luo_small_CG-DMRs" = Luo_small_CG_DMRs,
                        "Lister_CG-DMRs" = Lister_CG_DMRs),
  our_dmrs = list("POS_CG-DMRs" = POS_CG_DMRs,
                  "POSvsNEG_CG-DMRs" = POSvsNEG_CG_DMRs),
  SIMPLIFY = FALSE,
  mc.cores = getOption("mc.cores"))
others_differential_features <-
  c(list("Luo_small_CG-DMRs" = Luo_small_CG_DMRs,
         "Lister_CG-DMRs" = Lister_CG_DMRs))
# NOTE: Rename to avoid '.' in names
names(others_differential_features) <-
  gsub("\\.", "_", names(others_differential_features))

# ------------------------------------------------------------------------------
# CNS (central nervous system) cell type group
#

# NOTE: Annotation files for CNS were made available by the LDSC authors.
#       However, the BED files were not, so I construct an approximation to the
#       BED files based on the annotation files
CNS_annot <- mclapply(seqlevels, function(sl) {
  read_tsv(file.path(extdir, "Phase1", "cell_type_groups",
                     paste0("CNS.", sl, ".annot.gz")))
}, mc.cores = getOption("mc.cores"))
CNS_BED <- unlist(GRangesList(mclapply(CNS_annot, function(x) {
  y <- Rle(x$CNS)
  s <- cumsum(runLength(y))[runValue(y) == 1] -
    runLength(y)[runValue(y) == 1] + 1
  e <- cumsum(runLength(y))[runValue(y) == 1]
  GRanges(seqnames = paste0("chr", x[["CHR"]][1]),
          ranges = IRanges(x[["BP"]][s],
                           x[["BP"]][e]))
})))

# ------------------------------------------------------------------------------
# List of all brain-derived features
#

categories <- c(our_non_differential_features,
                our_differential_features,
                others_non_differential_features,
                others_differential_features,
                list("CNS" = CNS_BED))
stopifnot(all(sapply(categories, isDisjoint)))
stopifnot(all(!grepl("\\.", names(categories))))

saveRDS(categories, "../objects/categories.rds")

### ============================================================================
### Make annotations
###

# NOTE: Don't re-make the CNS annotation
k <- grep("CNS", names(categories), invert = TRUE)
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("../extdata/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(categories)[k]] <- mclapply(names(categories), function(cn) {
    stopifnot(isDisjoint(categories[[cn]]))
    as.integer(overlapsAny(cds_gr, categories[[cn]]))
  }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(categories[k]), function(n) {
    fl <- paste0("../output/ldsc/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)

# ------------------------------------------------------------------------------
# Annot file for 'base' category
# NOTE: Not added to 'categories' since it is not a brain-derived category
#

mclapply(seqlevels, function(sl) {
  x <- read_tsv(file.path(extdir, "Phase1", "baseline",
                          paste0("baseline.", sl, ".annot.gz")))
  write_tsv(x[, 1:5], paste0("../output/ldsc/base.Phase1.", sl, ".annot.gz"))
}, mc.cores = getOption("mc.cores"))

### ============================================================================
### Complex set operation features
###

POS_CG_DMRs <- categories[["POS_CG-DMRs"]]
chromHMM <- categories[["chromHMM_union"]]
CNS <- categories[["CNS"]]
H3K27ac <- categories[["H3K27ac"]]

intersect <- GenomicRanges::intersect
union <- GenomicRanges::union
setdiff <- GenomicRanges::setdiff

# Brain only
chromHMM_only <- setdiff(chromHMM, CNS_H3K27ac)
CNS_only <- setdiff(CNS, chromHMM_H3K27ac)
H3K27ac_only <- setdiff(H3K27ac, CNS_chromHMM)

# 2+
two_plus <- setdiff(
  setdiff(
    setdiff(union(CNS, union(H3K27ac, chromHMM)), chromHMM_only),
    CNS_only),
  H3K27ac_only)

# POS_CG-DMRs
CG_only <- setdiff(POS_CG_DMRs, union(chromHMM, union(CNS, H3K27ac)))
CG_shared <- intersect(POS_CG_DMRs, union(chromHMM, union(CNS, H3K27ac)))

# 2+ without CG_shared
two_plus_no_CG <- setdiff(two_plus, CG_shared)

# Brain only without CG_shared
chromHMM_only_no_CG <- setdiff(chromHMM_only, CG_shared)
CNS_only_no_CG <- setdiff(CNS_only, CG_shared)
H3K27ac_only_no_CG <- setdiff(H3K27ac_only, CG_shared)

complex_set_op_features <- list(CG_only = CG_only,
                                CG_shared = CG_shared,
                                chromHMM_only = chromHMM_only,
                                H3K27ac_only = H3K27ac_only,
                                CNS_only = CNS_only,
                                two_plus = two_plus,
                                two_plus_no_CG = two_plus_no_CG,
                                chromHMM_only_no_CG = chromHMM_only_no_CG,
                                CNS_only_no_CG = CNS_only_no_CG,
                                H3K27ac_only_no_CG = H3K27ac_only_no_CG)
saveRDS(complex_set_op_features, "../objects/complex_set_op_features.rds")

# Make annotations
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("../extdata/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(complex_set_op_features)] <-
    mclapply(names(complex_set_op_features), function(cn) {
      stopifnot(isDisjoint(complex_set_op_features[[cn]]))
      as.integer(overlapsAny(cds_gr, complex_set_op_features[[cn]]))
    }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(complex_set_op_features), function(n) {
    fl <- paste0("../output/ldsc/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)

### ============================================================================
### Non-differential features without each differential feature
###

differential_features <- categories[c("POS_CG-DMRs", "POSvsNEG_CG-DMRs",
                                      "CH-DMRs",
                                      "POS_DARs", "POSvsNEG_DARs",
                                      "POS_bigDARs", "POSvsNEG_bigDARs")]
non_differential_features <- categories[c("CNS", "H3K27ac", "chromHMM_union")]

n_df <- names(differential_features)
names(n_df) <- n_df
n_ndf <- names(non_differential_features)
names(n_ndf) <- n_ndf
ndf_excluding_df <- unlist(lapply(n_ndf, function(n1) {
  ndf <- non_differential_features[[n1]]
  lapply(n_df, function(n2) {
    df <- differential_features[[n2]]
    setdiff(ndf, df)
  })
}))
names(ndf_excluding_df) <- gsub("\\.", "_excluding_", names(ndf_excluding_df))

saveRDS(
  ndf_excluding_df,
  "../objects/non-differential_features_excluding_differential_features.rds")

# Make annotations
mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("../extdata/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(ndf_excluding_df)] <-
    mclapply(names(ndf_excluding_df), function(cn) {
      stopifnot(isDisjoint(ndf_excluding_df[[cn]]))
      as.integer(overlapsAny(cds_gr, ndf_excluding_df[[cn]]))
    }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(ndf_excluding_df), function(n) {
    fl <- paste0("../output/ldsc/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)
