# Make a GRanges object with unique
# (FANTOM5 TSS-enhancer associated enhancer, GENCODE v19 gene ID)-pairs
# Peter Hickey
# 2017-01-09

library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

###-----------------------------------------------------------------------------
### Import BED file
###

# Import the BED file as a UCSCData object
enhancer_tss_assocations <-
  rtracklayer::import(file.path(extdir, "enhancer_tss_associations.bed"))

# NOTE: X and Y appear in middle of autosomes for some reason, so need to
# re-order seqlevels
seqlevels(enhancer_tss_assocations, force = TRUE) <-
  paste0("chr", c(1:22, "X", "Y"))
# Sort by co-ordinates
enhancer_tss_assocations <- sort(enhancer_tss_assocations)

###-----------------------------------------------------------------------------
### Extract co-ordinates of putative enhancers and their linked TSS +/- 500 bp
### blocks
###

putative_enhancers <- as(
  vapply(X = strsplit(enhancer_tss_assocations$name, split = ";"),
         FUN = "[[",
         FUN.VALUE = character(1L), 1L),
  "GRanges")
start(putative_enhancers) <- start(putative_enhancers) + 1L
blocks <- mcols(enhancer_tss_assocations)[["blocks"]]

leftmost <- enhancer_tss_assocations
start(leftmost) <- vapply(X = start(leftmost) + start(blocks),
                          FUN = "[[",
                          FUN.VALUE = integer(1L),
                          1L) - 1L
width(leftmost) <- vapply(width(blocks), "[[", integer(1L), 1L) - 1L
leftmost <- granges(leftmost)

rightmost <- enhancer_tss_assocations
start(rightmost) <- vapply(X = start(rightmost) + start(blocks),
                           FUN = "[[",
                           FUN.VALUE = integer(1L),
                           2L)
rightmost <- granges(rightmost)

tss_is_leftmost <- width(leftmost) == 1000L
tss_block <- rightmost
tss_block[tss_is_leftmost] <- leftmost[tss_is_leftmost]
putative_enhancers$TSSBlock <- tss_block

###-----------------------------------------------------------------------------
### Add SeqInfo and drop non-autosomal seqlevels
###

seqinfo(putative_enhancers) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
putative_enhancers <- keepSeqlevels(putative_enhancers, paste0("chr", 1:22))

###-----------------------------------------------------------------------------
### Find GENCODE v19 protein-coding genes with a TSS overlapping the
### TSS +/- 500 bp block of the FANTOM5 enhancers list. Then, create a GRanges
### object with unique (putative enhancer, GENCODE gene ID) pairs
### NOTE: Requiring that the TSS matches the midpoint of the TSS +/- 500 bp
###       block is too stringent a criteria
###

load("../objects/unflattened-GENCODE-v19-features.rda")
tx_by_gene <- unflattened_features_pc_transcripts$transcripts_by_gene
tss_by_gene <- promoters(tx_by_gene, upstream = 0, downstream = 1)
ol <- findOverlaps(putative_enhancers$TSSBlock, tss_by_gene)
unique_pairs_idx <- lapply(as.list(ol), function(i) {
  !duplicated(names(tss_by_gene)[i])
})
putative_enhancers_with_gencode_links <- putative_enhancers[queryHits(ol)]
putative_enhancers_with_gencode_links$GENCODEGeneID <-
  names(tss_by_gene)[subjectHits(ol)]
putative_enhancers_with_gencode_links <-
  putative_enhancers_with_gencode_links[unlist(unique_pairs_idx)]
putative_enhancers_with_gencode_links <-
  sort(putative_enhancers_with_gencode_links)

saveRDS(putative_enhancers_with_gencode_links,
        "../objects/FANTOM5_enhancers_with_GENCODE_v19_links.rds")
