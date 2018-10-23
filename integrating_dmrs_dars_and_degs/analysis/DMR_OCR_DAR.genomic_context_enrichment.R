# Plot genomic context of CG-DMRs, CG-blocks, CH-DMRs, and DARs
# Peter Hickey
# 2017-10-30

# ==============================================================================
# Setup
#

library(GenomicRanges)
library(readr)
library(dplyr)
library(tidyr)
library(viridis)
library(gplots)
library(RColorBrewer)
library(rtracklayer)

# ------------------------------------------------------------------------------
# Load data
#

load("../objects/assays-and-features.rda")
source("../scripts/functions.R")
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
load("../../Objects/Annotated_POSvNEG_BLOCKs_GRanges.rda")
load("../../Objects/All_NeuNpos_nonNA_DMRs_fwer50.rda")
load("../objects/dmrs_pos_forEnrich.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_candidate_CH_DMRs <- c(
  list_of_candidate_CH_DMRs)

pos_CA <- readRDS("../../nonCG/extdata/flow-sorted-brain-wgbs/objects/pos_CA.loci.GRanges.rds")
neg_CA <- readRDS("../../nonCG/extdata/flow-sorted-brain-wgbs/objects/neg_CA.loci.GRanges.rds")
pos_CT <- readRDS("../../nonCG/extdata/flow-sorted-brain-wgbs/objects/pos_CT.loci.GRanges.rds")
neg_CT <- readRDS("../../nonCG/extdata/flow-sorted-brain-wgbs/objects/neg_CT.loci.GRanges.rds")

list_of_CH_loci <- list(
  "mCA (+)" = pos_CA,
  "mCA (-)" = neg_CA,
  "mCT (+)" = pos_CT,
  "mCT (-)" = neg_CT,
  "union" = reduce(unlist(GRangesList(pos_CA, neg_CA, pos_CT, neg_CT),
                          use.names = FALSE),
                   ignore.strand = TRUE))

# ------------------------------------------------------------------------------
# Construct list of CG-DMRs
#

list_of_CG_DMRs <- c(list(
  "CG-DMRs (POS)" = dmrs_pos,
  "CG-DMRs (POSvsNEG)" = Annotated_POSvNEG_DMRs_gr),
  setNames(as.list(split(dmrs_pos, dmrs_pos$cluster3)),
           paste0("CG-DMRs (POS, cluster ",
                  names(split(dmrs_pos, dmrs_pos$cluster3)), ")")))

# ------------------------------------------------------------------------------
# Construct list of CH-DMRs
#

list_of_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(dmrs) {
  dmrs[dmrs$fwer / dmrs$successful_permutations <= 0.05, ]
})
list_of_CH_DMRs <- c(
  list_of_CH_DMRs,
  list("CH-DMRs (union)" =
         reduce(unlist(GRangesList(list_of_CH_DMRs), use.names = FALSE),
                ignore.strand = TRUE)))

# ------------------------------------------------------------------------------
# Construct list of CG-blocks
#

list_of_CG_blocks <- list(
  "CG-blocks (POS)" = GRanges(sig_block_dmrs)[sig_block_dmrs$maxDiff > 0.1],
  "CG-blocks (POSvsNEG)" =
    Annotated_POSvNEG_BLOCKs_gr[Annotated_POSvNEG_BLOCKs_gr$maxDiff > 0.1])

# ------------------------------------------------------------------------------
# Load/construct bases in DARs and OCRs
#

list_of_DARs <- list(
  "DARs (POS)" = dars_pos,
  "DARs (POSvsNEG)" = dars_pos_vs_neg)

list_of_OCRs <- list(
  "union" = ocrs_overall,
  "NAcc (POS)" = ocrs_NAcc_pos,
  "NAcc (NEG)" = ocrs_NAcc_neg,
  "BA9 (POS)" = ocrs_BA9_pos,
  "BA9 (NEG)" = ocrs_BA9_neg)

# ------------------------------------------------------------------------------
# Construct intersection of DMRs and DARs
#

list_of_CG_DMR_and_DAR_intersection <- list(
  "CG-DMR DAR intersection (POS)" =
    GenomicRanges::intersect(list_of_CG_DMRs[["CG-DMRs (POS)"]],
              list_of_DARs[["DARs (POS)"]]),
  "CG-DMR DAR intersection (POSvsNEG)" =
    GenomicRanges::intersect(list_of_CG_DMRs[["CG-DMRs (POSvsNEG)"]],
              list_of_DARs[["DARs (POSvsNEG)"]]))

# ------------------------------------------------------------------------------
# Construct bodies of non-DEGs (POS)
#

non_degs <- unflattened_features$genes[non_deg_names]

# ------------------------------------------------------------------------------
# Construct promoters of DEGs and non-DEGs (POS)
#

deg_promoters <- promoters[match(deg_names, unlist(promoters$GENEID))]
non_deg_promoters <- promoters[match(non_deg_names, unlist(promoters$GENEID))]

# ------------------------------------------------------------------------------
# Load/construct features
#

# NOTE: All features are genic features are with respect to protein coding
#       genes in GENCODE v19
mySession <- browserSession("UCSC")
genome(mySession) <- "hg19"
repeat_masker <- getTable(ucscTableQuery(mySession,
                                         track = "rmsk",
                                         table = "rmsk"))
repeat_masker <- makeGRangesFromDataFrame(df = repeat_masker,
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "genoName",
                                          start.field = "genoStart",
                                          end.field = "genoEnd",
                                          strand = "strand")
repeat_masker <- sort(keepSeqlevels(repeat_masker, paste0("chr", 1:22),
                                    pruning.mode = "coarse"))
repeat_masker_by_repclass <- as.list(split(repeat_masker,
                                           repeat_masker$repClass))

features <- c(
  cgi_features,
  gencode_features[["pc_transcripts"]][c("promoter", "five_utr", "exonic",
                                         "intronic", "three_utr",
                                         "intergenic")],
  unlinked_enhancers[c("H3K27ac", "FANTOM5")],
  list("OCR (union)" = list_of_OCRs[["union"]]),
  list_of_CG_DMRs[c("CG-DMRs (POS)", "CG-DMRs (POSvsNEG)")],
  list("CH-DMRs (POS)" = list_of_CH_DMRs[["union"]]),
  repeat_masker_by_repclass[c("DNA", "LINE", "Low_complexity", "LTR",
                              "Satellite", "Simple_repeat", "SINE")])

# ==============================================================================
# Compute enrichment log odds ratios
#

CG_DMRs_OR_df <- bind_rows(
  lapply(names(list_of_CG_DMRs), function(name) {
    CG_DMRs <- list_of_CG_DMRs[[name]]
    bind_rows(FT(x = subsetByOverlaps(cpgs, CG_DMRs),
                 notx = subsetByOverlaps(cpgs, CG_DMRs, invert = TRUE),
                 features = features,
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                      c(degs, non_degs, ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                          invert = TRUE),
                                         c(degs, non_degs, ignore.mcols = TRUE)),
                 features = list("DEGs" = degs),
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                      c(deg_promoters, non_deg_promoters,
                                        ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                          invert = TRUE),
                                         c(deg_promoters, non_deg_promoters,
                                           ignore.mcols = TRUE)),
                 features = list("DEG promoters" = deg_promoters),
                 db = name))
  }))
saveRDS(CG_DMRs_OR_df, "../objects/CG_DMRs_OR_df.rds")

CG_blocks_OR_df <- bind_rows(
  lapply(names(list_of_CG_blocks), function(name) {
    CG_DMRs <- list_of_CG_blocks[[name]]
    bind_rows(FT(x = subsetByOverlaps(cpgs, CG_DMRs),
                 notx = subsetByOverlaps(cpgs, CG_DMRs, invert = TRUE),
                 features = features,
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                      c(degs, non_degs, ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                          invert = TRUE),
                                         c(degs, non_degs, ignore.mcols = TRUE)),
                 features = list("DEGs" = degs),
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs),
                                      c(deg_promoters, non_deg_promoters,
                                        ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMRs,
                                                          invert = TRUE),
                                         c(deg_promoters, non_deg_promoters,
                                           ignore.mcols = TRUE)),
                 features = list("DEG promoters" = deg_promoters),
                 db = name))
  }))
saveRDS(CG_blocks_OR_df, "../objects/CG_blocks_OR_df.rds")

# TODO: This breaks if CH-DMRs (union) is compared against itself
CH_DMRs_OR_df <- bind_rows(
  lapply(names(list_of_CH_DMRs), function(name) {
    CH_DMRs <- list_of_CH_DMRs[[name]]
    CH <- list_of_CH_loci[[name]]
    bind_rows(FT(x = subsetByOverlaps(CH, CH_DMRs),
                 notx = subsetByOverlaps(CH, CH_DMRs, invert = TRUE),
                 features = features,
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(CH, CH_DMRs),
                                      c(degs, non_degs, ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(CH, CH_DMRs,
                                                          invert = TRUE),
                                         c(degs, non_degs, ignore.mcols = TRUE)),
                 features = list("DEGs" = degs),
                 db = name),
              FT(x = subsetByOverlaps(subsetByOverlaps(CH, CH_DMRs),
                                      c(deg_promoters, non_deg_promoters,
                                        ignore.mcols = TRUE)),
                 notx = subsetByOverlaps(subsetByOverlaps(CH, CH_DMRs,
                                                          invert = TRUE),
                                         c(deg_promoters, non_deg_promoters,
                                           ignore.mcols = TRUE)),
                 features = list("DEG promoters" = deg_promoters),
                 db = name))
  }))
saveRDS(CH_DMRs_OR_df, "../objects/CH_DMRs_OR_df.rds")

# NOTE: This computates enrichment with respect to CpGs inside and outside of
#       CG-DMR:DAR intersection
CG_DMR_DAR_intersection_OR_df1 <- bind_rows(
  lapply(names(list_of_CG_DMR_and_DAR_intersection), function(name) {
    CG_DMR_and_DAR_intersection <- list_of_CG_DMR_and_DAR_intersection[[name]]
    bind_rows(FT(
      x = subsetByOverlaps(cpgs, CG_DMR_and_DAR_intersection),
      notx = subsetByOverlaps(cpgs, CG_DMR_and_DAR_intersection, invert = TRUE),
      features = features,
      db = name),
      FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMR_and_DAR_intersection),
                              c(degs, non_degs, ignore.mcols = TRUE)),
         notx = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMR_and_DAR_intersection,
                                                  invert = TRUE),
                                 c(degs, non_degs, ignore.mcols = TRUE)),
         features = list("DEGs" = degs),
         db = name),
      FT(x = subsetByOverlaps(subsetByOverlaps(cpgs, CG_DMR_and_DAR_intersection),
                              c(deg_promoters, non_deg_promoters,
                                ignore.mcols = TRUE)),
         notx = subsetByOverlaps(subsetByOverlaps(cpgs,
                                                  CG_DMR_and_DAR_intersection,
                                                  invert = TRUE),
                                 c(deg_promoters, non_deg_promoters,
                                   ignore.mcols = TRUE)),
         features = list("DEG promoters" = deg_promoters),
         db = name))

  }))
saveRDS(CG_DMR_DAR_intersection_OR_df1,
        "../objects/CG_DMR_DAR_intersection_OR_df1.rds")

# NOTE: This computates enrichment with respect to bases inside CG-DMR:DAR
#       intersection and inside null peaks
CG_DMR_DAR_intersection_OR_df2 <- bind_rows(
  lapply(names(list_of_CG_DMR_and_DAR_intersection), function(name) {
    CG_DMR_and_DAR_intersection <- list_of_CG_DMR_and_DAR_intersection[[name]]
    null_peaks <- subsetByOverlaps(list_of_OCRs[["union"]],
                                   CG_DMR_and_DAR_intersection,
                                   invert = TRUE)
    FT3(x = CG_DMR_and_DAR_intersection,
        y = null_peaks,
        features = features,
        db = name)
  })
)
saveRDS(CG_DMR_DAR_intersection_OR_df2,
        "../objects/CG_DMR_DAR_intersection_OR_df2.rds")

OCRs_vs_rest_of_genome_OR_df <- bind_rows(
  lapply(names(list_of_OCRs), function(name) {
    OCRs <- list_of_OCRs[[name]]
    FT2(x = OCRs,
        features = features,
        db = name,
        sl = sl)
  })
)
saveRDS(OCRs_vs_rest_of_genome_OR_df,
        "../objects/OCRs_vs_rest_of_genome_OR_df.rds")

DARs_vs_null_peaks_OR_df <- bind_rows(
  lapply(names(list_of_DARs), function(name) {
    DARs <- list_of_DARs[[name]]
    null_peaks <- subsetByOverlaps(list_of_OCRs[["union"]], DARs,
                                   invert = TRUE)
    FT3(x = DARs,
        y = null_peaks,
        features = features,
        db = name)
  })
)
saveRDS(DARs_vs_null_peaks_OR_df,
        "../objects/DARs_vs_null_peaks_OR_df.rds")

DARs_with_logFCgeq1_vs_null_peaks_OR_df <- bind_rows(
  lapply(names(list_of_DARs), function(name) {
    DARs <- list_of_DARs[[name]]
    DARs <- DARs[abs(DARs$logFC) > 1]
    null_peaks <- subsetByOverlaps(list_of_OCRs[["union"]], DARs,
                                   invert = TRUE)
    FT3(x = DARs,
        y = null_peaks,
        features = features,
        db = paste(name, "(|logFC| > 1)"))
  })
)
saveRDS(DARs_with_logFCgeq1_vs_null_peaks_OR_df,
        "../objects/DARs_with_logFCgeq1_vs_null_peaks_OR_df.rds")

# NOTE: Replace non-significant with logOR = NA
mCG_OR_df <- bind_rows(CG_DMRs_OR_df, CG_blocks_OR_df) %>%
  mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
  dplyr::select(db, feature, `log2(OR)`) %>%
  spread(feature, `log2(OR)`)
saveRDS(mCG_OR_df, "../objects/mCG_OR_df.rds")

ATAC_OR_df <- bind_rows(OCRs_vs_rest_of_genome_OR_df,
                        DARs_vs_null_peaks_OR_df,
                        DARs_with_logFCgeq1_vs_null_peaks_OR_df) %>%
  filter(!grepl("logFC", db)) %>%
  mutate(`log2(OR)` = ifelse(lower < 1 & upper > 1, NA, log2(estimate)),
         DAR = grepl("DAR", db)) %>%
  dplyr::select(db, feature, `log2(OR)`, DAR) %>%
  spread(feature, `log2(OR)`) %>%
  arrange(DAR)
saveRDS(ATAC_OR_df, "../objects/ATAC_OR_df.rds")

mCG_ATAC_intersection_OR_df1 <- bind_rows(CG_DMR_DAR_intersection_OR_df1,
                                          CG_DMRs_OR_df) %>%
  filter(!grepl("cluster", db)) %>%
  mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
  dplyr::select(db, feature, `log2(OR)`) %>%
  spread(feature, `log2(OR)`)
saveRDS(mCG_ATAC_intersection_OR_df1,
        "../objects/mCG_ATAC_intersection_OR_df1.rds")

mCG_ATAC_intersection_OR_df2 <- bind_rows(CG_DMR_DAR_intersection_OR_df2,
                                          DARs_vs_null_peaks_OR_df) %>%
  filter(!grepl("logFC", db)) %>%
  mutate(`log2(OR)` = ifelse(lower < 1 & upper > 1, NA, log2(estimate)),
         DAR = grepl("DAR", db)) %>%
  dplyr::select(db, feature, `log2(OR)`, DAR) %>%
  spread(feature, `log2(OR)`) %>%
  arrange(DAR)
saveRDS(mCG_ATAC_intersection_OR_df2,
        "../objects/mCG_ATAC_intersection_OR_df2.rds")

ATAC_OR_df2 <- rbind(ATAC_OR_df,mCG_ATAC_intersection_OR_df2[1:2,])

# ==============================================================================
# Heatmap of enrichment log odds ratios
#

hmcol <- colorRampPalette(c("blue", "white", "red"))(100)

mCG_OR_matrix_all <- mCG_OR_df %>%
  dplyr::select(-db, -`CG-DMRs (POS)`, -`CG-DMRs (POSvsNEG)`) %>%
  as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]

mCG_OR_matrix <- mCG_OR_matrix_all[!grepl("cluster",
                                          rownames(mCG_OR_matrix_all)), ]
mCG_clusters_OR_matrix <-
  mCG_OR_matrix_all[grepl("cluster", rownames(mCG_OR_matrix_all)), ]

ATAC_OR_matrix <- ATAC_OR_df %>%
  dplyr::select(`CG-DMRs (POS)`, `CG-DMRs (POSvsNEG)`, `CGI`,
                `exonic`, `FANTOM5`, `five_utr`, `H3K27ac`, `intergenic`,
                `intronic`, `OpenSea`, `promoter`, `Shelves`,
                `Shores`, `three_utr`) %>%
  as.matrix()
rownames(ATAC_OR_matrix) <- ATAC_OR_df[, 1]


ATAC_OR_matrix2 <- ATAC_OR_df2 %>%
        dplyr::select(`CG-DMRs (POS)`, `CG-DMRs (POSvsNEG)`, `CGI`,
                      `exonic`, `FANTOM5`, `five_utr`, `H3K27ac`, `intergenic`,
                      `intronic`, `OpenSea`, `promoter`, `Shelves`,
                      `Shores`, `three_utr`) %>%
        as.matrix()
rownames(ATAC_OR_matrix2) <- ATAC_OR_df2[, 1]



mCG_ATAC_intersection_OR_matrix1 <- mCG_ATAC_intersection_OR_df1 %>%
  dplyr::select(-db, -`CG-DMRs (POS)`, -`CG-DMRs (POSvsNEG)`,
                -`OCR (union)`) %>%
  as.matrix()
rownames(mCG_ATAC_intersection_OR_matrix1) <- mCG_ATAC_intersection_OR_df1[, 1]

mCG_ATAC_intersection_OR_matrix2 <- mCG_ATAC_intersection_OR_df2 %>%
  dplyr::select(`CG-DMRs (POS)`, `CG-DMRs (POSvsNEG)`, `CGI`,
                `exonic`, `FANTOM5`, `five_utr`, `H3K27ac`, `intergenic`,
                `intronic`, `OpenSea`, `promoter`, `Shelves`,
                `Shores`, `three_utr`) %>%
  as.matrix()
rownames(mCG_ATAC_intersection_OR_matrix2) <- mCG_ATAC_intersection_OR_df2[, 1]

pdf("../figures/DMR_OCR_DAR.genomic_context_enrichment.pdf")
heatmap.2(t(mCG_OR_matrix),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_OR_matrix),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_OR_matrix),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_clusters_OR_matrix),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_clusters_OR_matrix),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_clusters_OR_matrix),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(ATAC_OR_matrix),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(ATAC_OR_matrix),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(ATAC_OR_matrix),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix1),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix1),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix1),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix2),
          Rowv = NULL,
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix2),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(mCG_ATAC_intersection_OR_matrix2),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

heatmap.2(t(ATAC_OR_matrix2),
          Colv = NULL,
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")
dev.off()

