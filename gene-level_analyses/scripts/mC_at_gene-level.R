# mC levels at gene-level
# Peter Hickey
# 2018-03-22

### ----------------------------------------------------------------------------
### Setup
###

library(bsseq)
library(HDF5Array)
library(dplyr)
library(DelayedMatrixStats)

options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))

CG_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))
CG_BSseq <- CG_BSseq[, CG_BSseq$NeuN == "pos"]
CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "pos_CA_small-flow-sorted-brain-wgbs"))
CA_pos_BSseq <- CA_pos_BSseq[, CA_pos_BSseq$NeuN == "pos"]
CA_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "neg_CA_small-flow-sorted-brain-wgbs"))
CA_neg_BSseq <- CA_neg_BSseq[, CA_neg_BSseq$NeuN == "pos"]

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
mCAtGeneLevel <- function(genes, BSseq, n, min_cov = 1) {
  stopifnot(is(BSseq, "BSseq"))

  message("Computing meth")
  meth <- as.matrix(getMeth(BSseq, type = "smooth"))
  message("Applying coverage filter")
  keep <- rowMins(getCoverage(BSseq)) >= min_cov
  # NOTE: Constructing a SE is not strictly necessary, but does help
  #       keep GRanges and meth parallel
  message("Constructing initial SummarizedExperiment")
  se <- SummarizedExperiment(
    assay = meth[keep, , drop = FALSE],
    rowRanges = rowRanges(BSseq)[keep],
    colData = colData(BSseq))
  # Compute average methylation for each gene
  message("Finding overlaps with genes")
  ol <- findOverlaps(genes, se)
  ol_list <- as(ol, "List")
  m <- assay(se, withDimnames = FALSE)
  message("Computing methylation per gene")
  ave_meth <- do.call(rbind, lapply(ol_list, function(i) {
    colMeans2(m, rows = i, na.rm = TRUE)
  }))
  message("Constructing final SummarizedExperiment")
  val <- SummarizedExperiment(
    assay = ave_meth,
    rowRanges = genes,
    colData = colData(BSseq))
  assayNames(val) <- n
  val
}

### ----------------------------------------------------------------------------
### Compute and save results
###

# TODO: Report bug to Herve when running on CA_BSseq, specifically the
#       `meth <- as.matrix(getMeth(BSseq, type = "smooth"))` step
#       Error in validObject(.Object) :
#       invalid class “IRanges” object: 'start(x)', 'end(x)', and 'width(x)' cannot contain NAs
# In addition: Warning messages:
#   1: In successiveIRanges(rep.int(width, k), from = offset + 1L) :
#   integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'
# 2: In width(x) - 1L + start(x) : NAs produced by integer overflow
# Timing stopped at: 2916 228.8 3156
list_of_se <- list("mCG" = mCAtGeneLevel(genes, CG_BSseq, "mCG", 1),
                   "mCA_pos" = mCAtGeneLevel(genes, CA_pos_BSseq, "mCA", 1),
                   "mCA_neg" = mCAtGeneLevel(genes, CA_neg_BSseq, "mCA", 1))
se_mCA <- list_of_se[["mCA_pos"]]
neg_strand <- which(strand(se_mCA) == "-")
assay(se_mCA, withDimnames = FALSE)[neg_strand, ] <-
  assay(list_of_se[["mCA_neg"]], withDimnames = FALSE)[neg_strand, ]
se <- list_of_se[["mCG"]]
assay(se, "mCA", withDimnames = FALSE) <- assay(se_mCA, withDimnames = FALSE)

saveRDS(
  se,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "mC_at_gene_level.rds"))

# NOTE: Unstrand the genes so that methylation levels on the opposite strand
#       are estimate
list_of_se <- list(
  "mCG" = mCAtGeneLevel(unstrand(genes), CG_BSseq, "mCG", 1),
  "mCA_pos" = mCAtGeneLevel(unstrand(genes), CA_pos_BSseq, "mCA", 1),
  "mCA_neg" = mCAtGeneLevel(unstrand(genes), CA_neg_BSseq, "mCA", 1))
# TODO: Make a SE with 3 assays: mCG, mCA (same strand), mCA (opposite strand)
se <- list_of_se[["mCG"]]
strand(se) <- strand(genes)
neg_strand <- which(strand(se) == "-")
a <- assay(list_of_se[["mCA_pos"]], withDimnames = FALSE)
a[neg_strand, ] <- assay(
  list_of_se[["mCA_neg"]], withDimnames = FALSE)[neg_strand, ]
assay(se, "mCA (same strand)", withDimnames = FALSE) <- a
a <- assay(list_of_se[["mCA_neg"]], withDimnames = FALSE)
a[neg_strand, ] <- assay(
  list_of_se[["mCA_pos"]], withDimnames = FALSE)[neg_strand, ]
assay(se, "mCA (opposite strand)", withDimnames = FALSE) <- a

saveRDS(
  se,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "mC_at_gene_level.same_strand_and_opposite_strand.rds"))
