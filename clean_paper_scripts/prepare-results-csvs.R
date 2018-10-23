# Prepare CSV files of RNA-seq, ATAC-seq, and WGBS results
# Peter Hickey
# 2018-03-22

library(readr)
library(dplyr)
library(GenomicRanges)
outdir <- "~/Dropbox/Brain/SupplementaryTablesDrafts"

# NOTE: Don't do any rounding of P-values or logFC because this can change the
#       set of genes that satisfy the P-value or logFC cutoff

# RNA-seq (differential expression analysis)

NAvsBA9_pos_rna <-
  read_csv("../RNA-seq/extdata/topTable.NA_posvsBA9_pos.RNA-seq.csv.gz")
colnames(NAvsBA9_pos_rna)[c(1, 6)] <- c("gene_id", "BH_P")
NAvsBA9_pos_rna <- NAvsBA9_pos_rna[, c("gene_id", "logFC", "BH_P")]
write_csv(NAvsBA9_pos_rna, file.path(outdir, "NAvsBA9_pos.RNA-seq.csv"))

NAvsBA9_neg_rna <-
  read_csv("../RNA-seq/extdata/topTable.NA_negvsBA9_neg.RNA-seq.csv.gz")
colnames(NAvsBA9_neg_rna)[c(1, 6)] <- c("gene_id", "BH_P")
NAvsBA9_neg_rna <- NAvsBA9_neg_rna[, c("gene_id", "logFC", "BH_P")]
write_csv(NAvsBA9_neg_rna, file.path(outdir, "NAvsBA9_neg.RNA-seq.csv"))

# ATAC-seq (differential accessibility analysis)

NAvsBA9_pos_atac <-
  read_csv("../ATAC-seq/extdata/topTable.NA_posvsBA9_pos.ATAC-seq.csv.gz")
NAvsBA9_pos_atac <- NAvsBA9_pos_atac[, c(2, 3, 4, 7, 11)]
colnames(NAvsBA9_pos_atac)[c(1, 5)] <- c("chromosome", "BH_P")
NAvsBA9_pos_atac <- mutate(NAvsBA9_pos_atac,
                           logFC = signif(logFC, 3))
write_csv(NAvsBA9_pos_atac, file.path(outdir, "NAvsBA9_pos.ATAC-seq.csv"))

NAvsBA9_neg_atac <-
  read_csv("../ATAC-seq/extdata/topTable.NA_negvsBA9_neg.ATAC-seq.csv.gz")
NAvsBA9_neg_atac <- NAvsBA9_neg_atac[, c(2, 3, 4, 7, 11)]
colnames(NAvsBA9_neg_atac)[c(1, 5)] <- c("chromosome", "BH_P")
write_csv(NAvsBA9_neg_atac, file.path(outdir, "NAvsBA9_neg.ATAC-seq.csv"))

POSvsNEG_atac <-
  read_csv("../ATAC-seq/extdata/topTable.ave_pos_vs_ave_neg.ATAC-seq.csv.gz")
POSvsNEG_atac <- POSvsNEG_atac[, c(2, 3, 4, 7, 11)]
colnames(POSvsNEG_atac)[c(1, 5)] <- c("chromosome", "BH_P")
write_csv(POSvsNEG_atac, file.path(outdir, "POSvsNEG.ATAC-seq.csv"))

# WGBS (CG-DMRs)

load("../Objects/All_Annotated_DMRs_GRanges.rda")
dmrs_pos <- as.data.frame(Annotated_POS_DMRs_gr)[c(1:3, 9, 11, 13:17, 25:30)]
dmrs_pos$perm_P <- dmrs_pos$fwer / 1000
dmrs_pos$fwer <- NULL
j <- match(c("seqnames", "BA24", "BA9", "HC", "NA.", "NAvsBA9pos",
             "NAvsBA24pos", "NAvsHCpos", "BA24vsBA9pos", "HCvsBA9pos",
             "HCvsBA24pos"),
           colnames(dmrs_pos))
colnames(dmrs_pos)[j] <- c("chromosome", "BA24_pos", "BA9_pos", "HC_pos",
                           "NAcc_pos", "NAcc_pos_vs_BA9_pos",
                           "NAcc_pos_vs_BA24_pos", "NAcc_pos_vs_HC_pos",
                           "BA24_pos_vs_BA9_pos", "HC_pos_vs_BA9_pos",
                           "HC_pos_vs_BA24_pos")
write_csv(dmrs_pos, file.path(outdir, "CG-DMRs_pos.csv"))

dmrs_neg <- as.data.frame(Annotated_NEG_DMRs_gr)[c(1:3, 9, 11, 13:17, 25:30)]
dmrs_neg$perm_P <- dmrs_neg$fwer / 1000
dmrs_neg$fwer <- NULL
j <- match(c("seqnames", "BA24", "BA9", "HC", "NA.", "NAvsBA9neg",
             "NAvsBA24neg", "NAvsHCneg", "BA24vsBA9neg", "HCvsBA9neg",
             "HCvsBA24neg"),
           colnames(dmrs_neg))
colnames(dmrs_neg)[j] <- c("chromosome", "BA24_neg", "BA9_neg", "HC_neg",
                           "NAcc_neg", "NAcc_neg_vs_BA9_neg",
                           "NAcc_neg_vs_BA24_neg", "NAcc_neg_vs_HC_neg",
                           "BA24_neg_vs_BA9_neg", "HC_neg_vs_BA9_neg",
                           "HC_neg_vs_BA24_neg")
write_csv(dmrs_neg, file.path(outdir, "CG-DMRs_neg.csv"))

dmrs <- as.data.frame(Annotated_POSvNEG_DMRs_gr)[c(1:3, 9, 11, 13:21, 24:36)]
dmrs$perm_P <- dmrs$fwer / 1000
dmrs$fwer <- NULL
j <- match(c("seqnames", "NA_neg", "NA_pos", "NAvsBA9pos", "NAvsBA24pos",
             "NAvsHCpos", "BA24vsBA9pos", "HCvsBA9pos", "HCvsBA24pos",
             "NAvsBA9neg", "NAvsBA24neg", "NAvsHCneg", "BA24vsBA9neg",
             "HCvsBA9neg", "HCvsBA24neg", "POSvsNEG"),
           colnames(dmrs))
colnames(dmrs)[j] <- c("chromosome", "NAcc_neg", "NAcc_pos",
                       "NAcc_pos_vs_BA9_pos", "NAcc_pos_vs_BA24_pos",
                       "NAcc_pos_vs_HC_pos", "BA24_pos_vs_BA9_pos",
                       "HC_pos_vs_BA9_pos", "HC_pos_vs_BA24_pos",
                       "NAcc_neg_vs_BA9_neg", "NAcc_neg_vs_BA24_neg",
                       "NAcc_neg_vs_HC_neg", "BA24_neg_vs_BA9_neg",
                       "HC_neg_vs_BA9_neg", "HC_neg_vs_BA24_neg", "pos_vs_neg")
write_csv(dmrs, file.path(outdir, "CG-DMRs.csv"))

load("../Objects/All_NeuNpos_nonNA_DMRs_fwer50.rda")
dmrs <- as.data.frame(sig_pos_nonNA_dmrs_gr)
dmrs$perm_P <- dmrs$fwer / 1000
dmrs$fwer <- NULL
j <- match(c("seqnames"), colnames(dmrs))
colnames(dmrs)[j] <- "chromosome"
write_csv(dmrs, file.path(outdir, "non-NAcc_CG-DMRs.csv"))

# WGBS (CG-blocks)

load("../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
load("../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
blocks_pos <- sig_block_dmrs[, c(1:3, 7, 10, 12:16)]
blocks_pos$perm_P <- blocks_pos$fwer / 1000
blocks_pos$fwer <- NULL
j <- match(c("chr", "BA24", "BA9", "HC", "NA"), colnames(blocks_pos))
colnames(blocks_pos)[j] <- c("chromosome", "BA24_pos", "BA9_pos", "HC_pos",
                             "NAcc_pos")
ol <- findOverlaps(query = unflattened_features_pc_transcripts$genes,
                   subject = GRanges(blocks_pos),
                   type = "within")
blocks_pos$coversPCGene <- countSubjectHits(ol) > 0
write_csv(blocks_pos, file.path(outdir, "CG-blocks_pos.csv"))

load("../Objects/Annotated_POSvNEG_BLOCKs_GRanges.rda")
blocks <- as.data.frame(Annotated_POSvNEG_BLOCKs_gr)[, c(1:3, 9, 11, 13:21,
                                                         24:36)]
blocks$perm_P <- blocks$fwer / 1000
blocks$fwer <- NULL
j <- match(c("seqnames", "NA_neg", "NA_pos", "NAvsBA9pos", "NAvsBA24pos",
             "NAvsHCpos", "BA24vsBA9pos", "HCvsBA9pos", "HCvsBA24pos",
             "NAvsBA9neg", "NAvsBA24neg", "NAvsHCneg", "BA24vsBA9neg",
             "HCvsBA9neg", "HCvsBA24neg", "POSvsNEG"),
           colnames(blocks))
colnames(blocks)[j] <- c("chromosome", "NAcc_neg", "NAcc_pos",
                       "NAcc_pos_vs_BA9_pos", "NAcc_pos_vs_BA24_pos",
                       "NAcc_pos_vs_HC_pos", "BA24_pos_vs_BA9_pos",
                       "HC_pos_vs_BA9_pos", "HC_pos_vs_BA24_pos",
                       "NAcc_neg_vs_BA9_neg", "NAcc_neg_vs_BA24_neg",
                       "NAcc_neg_vs_HC_neg", "BA24_neg_vs_BA9_neg",
                       "HC_neg_vs_BA9_neg", "HC_neg_vs_BA24_neg", "pos_vs_neg")
write_csv(blocks, file.path(outdir, "CG-blocks.csv"))

# WGBS (CA-DMRs and CT-DMRs)
list_of_candidate_CH_DMRs <- readRDS(
  "../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds"
)
list_of_CH_DMRs <- lapply(
  list_of_candidate_CH_DMRs,
  function(x) {
    x[x$fwer <= 50]
  }
)
strand(list_of_CH_DMRs[["mCA (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCA (-)"]]) <- "-"
strand(list_of_CH_DMRs[["mCT (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCT (-)"]]) <- "-"
# CA-DMRs
CA_DMRs <- unname(
  c(list_of_CH_DMRs[["mCA (+)"]], list_of_CH_DMRs[["mCA (-)"]]))
CA_DMRs_df <- as.data.frame(
  CA_DMRs)[, c("seqnames", "strand", "start", "end", "n", "areaStat",
               "BA24_pos", "BA9_pos", "HC_pos", "NA_pos", "fwer")]
CA_DMRs_df$perm_P <- CA_DMRs_df$fwer / 1000
CA_DMRs_df$fwer <- NULL
colnames(CA_DMRs_df)[match(c("seqnames", "NA_pos"), colnames(CA_DMRs_df))] <-
  c("chromosome", "NAcc_pos")
write_csv(CA_DMRs_df, file.path(outdir, "CA-DMRs_pos.csv"))

CT_DMRs <- unname(
  c(list_of_CH_DMRs[["mCT (+)"]], list_of_CH_DMRs[["mCT (-)"]]))
CT_DMRs_df <- as.data.frame(
  CT_DMRs)[, c("seqnames", "start", "end", "n", "areaStat",
               "BA24_pos", "BA9_pos", "HC_pos", "NA_pos", "fwer")]
CT_DMRs_df$perm_P <- CT_DMRs_df$fwer / 1000
CT_DMRs_df$fwer <- NULL
colnames(CT_DMRs_df)[match(c("seqnames", "NA_pos"), colnames(CT_DMRs_df))] <-
  c("chromosome", "NAcc_pos")
write_csv(CT_DMRs_df, file.path(outdir, "CT-DMRs_pos.csv"))
