# ATAC-seq read counts in bins along genome
# 2018-03-13

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicAlignments)
library(dplyr)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

options("mc.cores" = 11)

binsize <- 10 ^ 4
si <- keepSeqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg19),
                    paste0("chr", 1:22))
bins <- tileGenome(seqlengths(si),
                   tilewidth = binsize,
                   cut.last.tile.in.chrom = TRUE)

### ----------------------------------------------------------------------------
### Load data
###

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
  se <- summarizeOverlaps(bins, gr, ignore.strand = TRUE)
  colnames(se) <- basename(f)
  se
}))
se$TISSUE_NEUN <- ifelse(grepl("NA", colnames(se)), "NAcc_pos", "BA9_pos")

dgelist <- DGEList(assay(se), samples = colData(se))
rpkm_by_group <- rpkmByGroup(
  dgelist,
  group = dgelist$samples$TISSUE_NEUN,
  gene.length = width(bins))
se <- SummarizedExperiment(
  rpkm_by_group,
  rowRanges = bins)

saveRDS(se, "../objects/ATAC_in_bins.rds")
