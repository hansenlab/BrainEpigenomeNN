# Smooth non-CG data
# Peter Hickey
# 2017-06-11

### ============================================================================
### Setup
###

library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

options("mc.cores" = 11)

# ------------------------------------------------------------------------------
# Load BSseq objects
#

pos_CH_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/pos_CH-flow-sorted-brain-wgbs")
neg_CH_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/neg_CH-flow-sorted-brain-wgbs")
unstranded_CH_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/unstranded_CH-flow-sorted-brain-wgbs")

# NOTE: The BSseq objects contain CH loci from which we now want to select CA
#       or CT loci. Furthermore, the BSseq objects are HDF5Array-backed. Rather
#       than subset the BSseq object then load the data from the relevant loci,
#       it is **much** faster to load the data and then subset (the trade off
#       is increased memory usage). Therefore, I use this (slightly convoluted)
#       strategy in what follows.

# ------------------------------------------------------------------------------
# Functions
#

smooth <- function(strand, context, BSseq, loci, ns, h, size) {
  message(strand)
  list_of_BSseq <- lapply(seqlevels(BSseq), function(sl) {
      message(sl)
    BSseq <- keepSeqlevels(BSseq, sl, pruning.mode = "coarse")
    message("Loading assays into memory")
    assays(BSseq) <- endoapply(assays(BSseq), as.matrix)
    message("Subsetting to just ", context, " loci")
    BSseq <- subsetByOverlaps(BSseq, loci, type = "equal")
    message("Smoothing")
    BSseq <- BSmooth(BSseq = BSseq,
                     ns = ns,
                     h = h,
                     maxGap = 10 ^ 8,
                     parallelBy = "sample",
                     mc.preschedule = FALSE,
                     mc.cores = getOption("mc.cores"),
                     keep.se = FALSE,
                     verbose = 2L)
    message("Saving BSseq object")
    saveHDF5SummarizedExperiment(
      x = BSseq,
      dir = file.path("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects",
                      paste0(strand, "_", context, "_", size,
                             "-flow-sorted-brain-wgbs.", sl)),
      verbose = TRUE)
  })
  BSseq <- do.call(rbind, list_of_BSseq)
  saveHDF5SummarizedExperiment(
    x = BSseq,
    dir = file.path("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects",
                    paste0(strand, "_", context, "_", size,
                           "-flow-sorted-brain-wgbs")),
    verbose = TRUE)
}

### ============================================================================
### mCA
###

pos_CA <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("CA", s)))
pos_CA_gr <- resize(GRanges(seqnames = Rle(names(pos_CA), lengths(pos_CA)),
                            ranges = unlist(as(pos_CA, "IRangesList"),
                                            use.names = FALSE),
                            strand = "+",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
rm(pos_CA)

neg_CA <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("TG", s)))
neg_CA_gr <- resize(GRanges(seqnames = Rle(names(neg_CA), lengths(neg_CA)),
                            ranges = unlist(as(neg_CA, "IRangesList"),
                                            use.names = FALSE),
                            strand = "-",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
rm(neg_CA)

# ------------------------------------------------------------------------------
# Pos strand (small smooth)
#

# NOTE: This errored after 43 hours on chr1 using ns = 150, h = 1000
# NOTE: This errored using 15 cores with "Error in names(object) <- nm : 'names'
#       attribute [45] must be the same length as the vector [44] (ns = 200,
#       h = 3000)
# NOTE: Succeeded on chr1 after 14.5 hours using ns = 200, h = 3000 with 11
#       cores
pos_CA_BSseq_small_smooth <- smooth(strand = "pos",
                                    context = "CA",
                                    BSseq = pos_CH_BSseq,
                                    loci = pos_CA_gr,
                                    ns = 200,
                                    h = 3000,
                                    size = "small")

# ------------------------------------------------------------------------------
# Neg strand (small smooth)
#

# NOTE: This worked on chr1 after 14 hours with 11 cores
neg_CA_BSseq_small_smooth <- smooth(strand = "neg",
                                    context = "CA",
                                    BSseq = neg_CH_BSseq,
                                    loci = neg_CA_gr,
                                    ns = 200,
                                    h = 3000,
                                    size = "small")

# ------------------------------------------------------------------------------
# Unstranded (small smooth)
#

# TODO

### ============================================================================
### mCT
###

pos_CT <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("CT", s)))
pos_CT_gr <- resize(GRanges(seqnames = Rle(names(pos_CT), lengths(pos_CT)),
                            ranges = unlist(as(pos_CT, "IRangesList"),
                                            use.names = FALSE),
                            strand = "+",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
rm(pos_CT)

neg_CT <- bsapply(new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19,
                      FUN = function(s) matchPattern("AG", s)))
neg_CT_gr <- resize(GRanges(seqnames = Rle(names(neg_CT), lengths(neg_CT)),
                            ranges = unlist(as(neg_CT, "IRangesList"),
                                            use.names = FALSE),
                            strand = "-",
                            seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19)),
                    width = 1,
                    fix = "start")
rm(neg_CT)

# ------------------------------------------------------------------------------
# Pos strand (small smooth)
#

# NOTE: Succeeded on chr1 after 14 hours using ns = 200, h = 3000 with 11 cores
pos_CT_BSseq_small_smooth <- smooth(strand = "pos",
                                    context = "CT",
                                    BSseq = pos_CH_BSseq,
                                    loci = pos_CT_gr,
                                    ns = 200,
                                    h = 3000)

# ------------------------------------------------------------------------------
# Neg strand (small smooth)
#

# NOTE: Succeeded on chr1 after 14 hours using ns = 200, h = 3000 with 11 cores
neg_CT_BSseq_small_smooth <- smooth(strand = "neg",
                                    context = "CT",
                                    BSseq = neg_CH_BSseq,
                                    loci = neg_CT_gr,
                                    ns = 200,
                                    h = 3000)

# ------------------------------------------------------------------------------
# Unstranded (small smooth)
#

# TODO
