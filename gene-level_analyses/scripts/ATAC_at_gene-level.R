# ATAC-seq read counts at gene-level
# Peter Hickey
# 2018-03-07

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicAlignments)
library(dplyr)
library(edgeR)

extdir <- "../extdata"

options("mc.cores" = 11)

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
gr_files <- grep(
  "[0-9]+-(NA|BA9)-pos\\.rep1\\_gr\\.rds",
  list.files(file.path(extdir, "flow-sorted-brain-atac", "objects"),
             full.names = TRUE),
  value = TRUE)

### ----------------------------------------------------------------------------
### Compute and save results
###

se <- do.call(cbind, mclapply(gr_files, function(f) {
  gr <- readRDS(f)
  gr <- gr[!gr$isDuplicate]
  se <- summarizeOverlaps(genes, gr, ignore.strand = TRUE)
  colnames(se) <- basename(f)
  se
}))
se$TISSUE_NEUN <- ifelse(grepl("NA", colnames(se)), "NAcc_pos", "BA9_pos")

dgelist <- DGEList(assay(se), samples = colData(se))
dgelist <- dgelist[names(genes), ]
rpkm_by_group <- rpkmByGroup(
  dgelist,
  group = dgelist$samples$TISSUE_NEUN,
  gene.length = width(genes))
se <- SummarizedExperiment(
  rpkm_by_group,
  rowRanges = genes)

saveRDS(
  se,
  file.path(extdir,
            "flow-sorted-brain-gene-level_analyses", "objects",
            "ATAC_at_gene_level.rds"))
