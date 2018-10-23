# Construct data needed for gene-level analyses
# Peter Hickey
# 2018-03-08

library(GenomicRanges)
library(limma)
library(dplyr)
library(edgeR)
library(tximport)
library(sva)

extdir <- "../extdata"

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
big_dars_pos <- dars_pos[abs(dars_pos$logFC) > 1]
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_CH_DMRs <- lapply(
  X = list_of_candidate_CH_DMRs,
  FUN = function(x) x[x$fwer <= 50])
list_of_CH_DMRs <- lapply(list_of_CH_DMRs, function(x) {
  x$meanDiff <- x$NA_pos - x$BA9_pos
  x
})
strand(list_of_CH_DMRs[["mCA (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCA (-)"]]) <- "-"
strand(list_of_CH_DMRs[["mCT (+)"]]) <- "+"
strand(list_of_CH_DMRs[["mCT (-)"]]) <- "-"
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")
blocks_pos <- GRanges(sig_block_dmrs)
blocks_pos$meanDiff <- blocks_pos$`NA` - blocks_pos$BA9

# Exons and introns by gene
tx_by_gene <- unflattened_features$transcripts_by_gene
# Exons
unlisted_tx_by_gene <- unlist(tx_by_gene)
exons_by_tx <- unflattened_features$exons_by_transcript
unlisted_exons_by_tx <- unlist(exons_by_tx)
i <- match(names(unlisted_exons_by_tx), unlisted_tx_by_gene$tx_name)
exons_by_gene <- split(unlisted_exons_by_tx, names(unlisted_tx_by_gene)[i])
# Introns
introns_by_tx <- unflattened_features$introns_by_transcript
unlisted_introns_by_tx <- unlist(introns_by_tx)
i <- match(names(unlisted_introns_by_tx), unlisted_tx_by_gene$tx_name)
introns_by_gene <- split(unlisted_introns_by_tx, names(unlisted_tx_by_gene)[i])
# NOTE: Add empty GRanges() for genes without introns
intronless_genes <- setdiff(names(exons_by_gene), names(introns_by_gene))
genes_without_introns <- setNames(GRangesList(
  lapply(seq_along(intronless_genes), function(i) GRanges())),
  intronless_genes)
introns_by_gene <- c(introns_by_gene, genes_without_introns)

# All GENCODE genes with logFC for NAcc_pos vs. BA9_pos comparison
# NOTE: Will be NA if gene not tested for differential expression
marraylm <- readRDS("../../RNA-seq/objects/fit_with_sv.rds")
genes <- gencode_features$genes
genes$logFC <- NA_real_
tt <- topTable(marraylm, "NA_posvsBA9_pos", number = Inf)
idx <- match(rownames(tt), names(genes))
genes$logFC[idx] <- tt$logFC

# Compute measures of gene expression in NAcc_pos and BA9_pos
txi_gene <- readRDS(
  "../../RNA-seq/objects/txi-gene.flow-sorted-brain-rna-seq.rds")
txi_gene[1:3] <- lapply(txi_gene[1:3], function(df) {
  df[rownames(df) %in% names(gencode_features$transcripts_by_gene), ]
})
cd <- readRDS("../../RNA-seq/objects/colData-flow-sorted-brain-rna-seq.rds")
# Add combination of TISSUE and NEUN as its own variable
cd$TISSUE_NEUN <- factor(paste0(cd$TISSUE, "_", cd$NEUN))
dgelist <- DGEList(txi_gene$counts, samples = cd)
dgelist <- dgelist[names(genes), ]
rpkm_by_group <- rpkmByGroup(
  dgelist,
  group = dgelist$samples$TISSUE_NEUN,
  gene.length = width(gencode_features$genes[rownames(dgelist)]))
ave_rpkm <- rpkmByGroup(
  dgelist[, dgelist$samples$TISSUE_NEUN %in% c("NA_pos", "BA9_pos")],
  group = rep("NA_pos_and_BA9_pos",
              sum(dgelist$samples$TISSUE_NEUN %in% c("NA_pos", "BA9_pos"))),
  gene.length = width(gencode_features$genes[rownames(dgelist)]))

# Whether gene overlaps another gene (ignoring strand)
overlaps_another_gene <- overlapsAny(
  genes, drop.self = TRUE, ignore.strand = TRUE)

genes_df <- data_frame(
  gene = names(genes),
  logFC = genes$logFC,
  rpkm_NAcc_pos = rpkm_by_group[, "NA_pos"],
  rpkm_BA9_pos = rpkm_by_group[, "BA9_pos"],
  ave_rpkm = ave_rpkm[, "NA_pos_and_BA9_pos"],
  DE = names(genes) %in% deg_names,
  CGI_promoter = names(genes) %in% names(cgi_promoters_by_gene),
  width = width(genes),
  gene_type = ifelse(names(genes) %in% pc_genes,
                     "PC",
                     ifelse(names(genes) %in% lnc_genes, "lncRNA", NA)),
  overlaps_another_gene)

save(genes, promoters_by_gene, exons_by_gene, introns_by_gene,
     dmrs_NAvsBA9pos, blocks_pos, dars_pos, big_dars_pos, list_of_CH_DMRs,
     genes_df,
     file =
       file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
                 "gene-level_analyses_data.rda"))
