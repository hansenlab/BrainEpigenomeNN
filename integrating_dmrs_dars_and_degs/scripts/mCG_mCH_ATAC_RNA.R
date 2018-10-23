# Integrating mCG, mCH, ATAC-seq, and RNA-seq data
# Peter Hickey
# 2017-12-06

library(SummarizedExperiment)
library(dplyr)

# ------------------------------------------------------------------------------
# Load data
#

load("../objects/assays-and-features.rda")
list_of_SE.DMRs_and_blocks <-
  readRDS("../../nonCG/objects/list_of_SEs.DMRs_and_blocks.rds")
gene_body_mC_gene_level <-
  readRDS("../../nonCG/objects/gene_body_mC.gene_level.rds")

# TODO: A promoter_mC_promoter_level object
# TODO: An enhancer_mC_enhancer_level object
# TODO: Some sort of
#       {gene, promoter, enhancer}_ATAC_{gene, promoter, enhancer}_level objects

# ------------------------------------------------------------------------------
# Construct genomic features
#

# TODO: CH_DMRs may not be strictly differential between NAcc and BA9
CH_DMRs <- Reduce(GenomicRanges::union,
                  lapply(list_of_SE.DMRs_and_blocks[3:6],
                         function(x) rowRanges(x[["sample_level"]])))
diff_epi <- list("CG-DMRs" = dmrs_pos,
                 "DARs" = dars_pos,
                 "CH-DMRs" = CH_DMRs)
diff_epi_union <- Reduce(GenomicRanges::union, diff_epi)

# ------------------------------------------------------------------------------
# Construct SE
#

# NOTE: Only includes genes tested for DE
se <- SummarizedExperiment(
  assays = lapply(gene_body_mC_gene_level, function(x) {
    x[rna_atac_meth$gene, ]
  }),
  rowRanges = gencode_features$genes[rna_atac_meth$gene])
rowData(se) <- inner_join(rna_atac_meth,
                          data_frame(gene = rownames(se))) %>%
  dplyr::select(DE, expLogFC, gene_symbol)

# ------------------------------------------------------------------------------
# Make table of DEGs with/without CG-DMRs, CH-DMRs, DARs
#

gene_body_ol_tbl <- sapply(diff_epi, function(x) {
  overlapsAny(degs[deg_names], x, ignore.strand = TRUE)
})
colnames(gene_body_ol_tbl) <- paste0("gb", colnames(gene_body_ol_tbl))
promoter_ol_tbl <- sapply(diff_epi, function(x) {
  overlapsAny(promoters_by_gene[deg_names], x, ignore.strand = TRUE)
})
colnames(promoter_ol_tbl) <- paste0("p", colnames(promoter_ol_tbl))
enhancer_ol_tbl <- sapply(diff_epi, function(x) {
  overlapsAny(fantom5_enhancers_by_gene_all_genes[deg_names], x,
              ignore.strand = TRUE)
})
colnames(enhancer_ol_tbl) <- paste0("e", colnames(enhancer_ol_tbl))
ol_tbl <- as_data_frame(cbind(
  data_frame(gene_id = deg_names),
  gene_body_ol_tbl, promoter_ol_tbl, enhancer_ol_tbl,
  data_frame(hasEnhancer = deg_names %in% names(fantom5_enhancers_by_gene))))

# DEGs with a CG-DMR, DAR, or CH-DMR in gene body, promoter, or enhancer
filter(ol_tbl,
       `gbCH-DMRs` | `pCH-DMRs` | `eCH-DMRs`,
       !`gbCG-DMRs`,
       !`pCG-DMRs`,
       !`eCG-DMRs`,
       !`gbDARs`,
       !`pDARs`,
       !`eDARs`) %>%
  pull(gene_id) %>%
  saveRDS(file = "../tmp/DEGs_with_CH-DMRs_but_no_CG-DMR_or_DAR.rds")

# DEGs without a CG-DMR, DAR, or CH-DMR in gene body, promoter, or enhancer
filter(ol_tbl,
       !`gbCH-DMRs`,
       !`pCH-DMRs`,
       !`eCH-DMRs`,
       !`gbCG-DMRs`,
       !`pCG-DMRs`,
       !`eCG-DMRs`,
       !`gbDARs`,
       !`pDARs`,
       !`eDARs`) %>%
  pull(gene_id) %>%
  saveRDS(file = "../tmp/DEGs_without_CG-DMR_or_DAR_or_CH-DMR.rds")

# DEGs without a CG-DMR, DAR, or CH-DMR in gene body, promoter, or enhancer
# **that have an enhancer**
filter(ol_tbl,
       !`gbCH-DMRs`,
       !`pCH-DMRs`,
       !`eCH-DMRs`,
       !`gbCG-DMRs`,
       !`pCG-DMRs`,
       !`eCG-DMRs`,
       !`gbDARs`,
       !`pDARs`,
       !`eDARs`,
       hasEnhancer)

# ------------------------------------------------------------------------------
# Plot RNA-seq logFC vs. gene body mC
#

par(mfrow = c(2, 3))
lapply(assayNames(se), function(an) {
  a <- assay(se, an)
  cor <- cor(a[, "NAcc"] - a[, "BA9"], rowData(se)$expLogFC,
             use = "complete.obs")
  plot(a[, "NAcc"] - a[, "BA9"], rowData(se)$expLogFC,
       main = paste0(an, ": cor = ", round(cor, 2)))
})

# ------------------------------------------------------------------------------
# Are genes with a DiffEpi (CG-DMR, CH-DMRs, or DAR) in the gene body more
# likely to be differentially expressed than not?
#

data_frame(diff_mark = c("CG-DMRs", "DARs", "CH-DMRs", "Any"),
           gene_feature = c("body"),
           OR = c(fisher.test(
             do.call(rbind, tapply(X = rowData(se)$DE,
                                   INDEX = overlapsAny(se, dmrs_pos),
                                   FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se)$DE,
                                     INDEX = overlapsAny(se, dars_pos),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se)$DE,
                                     INDEX = overlapsAny(se, CH_DMRs),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se)$DE,
                                     INDEX = overlapsAny(se, diff_epi),
                                     FUN = table)))$estimate))

# ------------------------------------------------------------------------------
# Are genes with a DiffEpi (CG-DMR, CH-DMRs, or DAR) in a promoter more
# likely to be differentially expressed than not?
#

# TODO: Use the SummarizedExperiment object with promoter level mC (not yet
#       made)
se2 <- se
rowRanges(se2) <- promoters(gencode_features$transcripts_by_gene[rownames(se)],
                            2000, 2000)
rowData(se2) <- rowData(se)

data_frame(diff_mark = c("CG-DMRs", "DARs", "CH-DMRs", "Any"),
           gene_feature = c("promoter"),
           OR = c(fisher.test(
             do.call(rbind, tapply(X = rowData(se2)$DE,
                                   INDEX = overlapsAny(se2, dmrs_pos),
                                   FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se2)$DE,
                                     INDEX = overlapsAny(se2, dars_pos),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se2)$DE,
                                     INDEX = overlapsAny(se2, CH_DMRs),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se2)$DE,
                                     INDEX = overlapsAny(se2, diff_epi),
                                     FUN = table)))$estimate))

# ------------------------------------------------------------------------------
# Are genes with a DiffEpi (CG-DMR, CH-DMRs, or DAR) in an enhancer more
# likely to be differentially expressed than not?
#

# TODO: Use the SummarizedExperiment object with enhancer level mC (not yet
#       made)
se3 <- se[rownames(se) %in% names(fantom5_enhancers_by_gene), ]
rowRanges(se3) <- fantom5_enhancers_by_gene[rownames(se3)]
rowData(se3) <- rowData(
  se[rownames(se) %in% names(fantom5_enhancers_by_gene), ])

data_frame(diff_mark = c("CG-DMRs", "DARs", "CH-DMRs", "Any"),
           gene_feature = c("enhancer"),
           OR = c(fisher.test(
             do.call(rbind, tapply(X = rowData(se3)$DE,
                                   INDEX = overlapsAny(se3, dmrs_pos),
                                   FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se3)$DE,
                                     INDEX = overlapsAny(se3, dars_pos),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se3)$DE,
                                     INDEX = overlapsAny(se3, CH_DMRs),
                                     FUN = table)))$estimate,
             fisher.test(
               do.call(rbind, tapply(X = rowData(se3)$DE,
                                     INDEX = overlapsAny(se3, diff_epi),
                                     FUN = table)))$estimate))

# ------------------------------------------------------------------------------
# Are genes with a DiffEpi (CG-DMR, CH-DMRs, or DAR) in the gene body, a
# promoter, or an enhancer more likely to be differentially expressed than not?
#

# TODO: Compute ORs
# TODO: Find DEGs without any diff_epi marks
