# mC levels at first exon, first intron, rest of exons, rest of introns
# Peter Hickey
# 2018-03-06

### ----------------------------------------------------------------------------
### Setup
###

library(bsseq)
library(HDF5Array)
library(dplyr)
library(tidyr)
library(DelayedMatrixStats)

options("mc.cores" = 2)
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
CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "pos_CA_small-flow-sorted-brain-wgbs"))
CA_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "neg_CA_small-flow-sorted-brain-wgbs"))

conditions_df <- data_frame(condition = c("NAcc_pos", "BA9_pos"),
                            tissue = c("NA", "BA9"),
                            NeuN = c("pos", "pos"))

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
mCAtGeneFeatures <- function(genes, promoters_by_gene, exons_by_gene,
                             introns_by_gene, BSseq,
                             conditions_df, min_cov = 1) {
  stopifnot(is(BSseq, "BSseq"))

  # Ensure genes, exons, and introns are all ordered the same way and that
  # exons and introns within each gene are sorted
  message("Sorting genes, promoters, exons, and introns")
  o <- names(genes)
  promoters_by_gene <- sort(promoters_by_gene[o])
  exons_by_gene <- sort(exons_by_gene[o])
  introns_by_gene <- sort(introns_by_gene[o])
  strand_is_minus <- strand(genes) == "-"

  # NOTE: Internal helper function
  .mCAtGeneFeature <- function(k, features) {
    message("Subsetting BSseq: ", conditions_df$condition[k])
    BSseq <- BSseq[, BSseq$Tissue %in% conditions_df$tissue[k] &
                     BSseq$NeuN %in% conditions_df$NeuN[k]]
    message("Computing meth: ", conditions_df$condition[k])
    meth <- as.matrix(getMeth(BSseq, type = "smooth"))
    keep <- rowMins(getCoverage(BSseq)) >= min_cov
    # NOTE: Constructing a SE is not strictly necessary, but does help
    #       keep GRanges and meth parallel
    se <- SummarizedExperiment(
      assay = as.matrix(rowMeans2(meth[keep, , drop = FALSE])),
      rowRanges = rowRanges(BSseq)[keep],
      colData = DataFrame(row.names = conditions_df$condition[k]))
    # Compute average methylation in each feature
    bind_rows(lapply(names(features), function(n) {
      message("\nComputing mC at ", n, "L ")
      feature <- features[[n]]
      stopifnot(is(feature, "GRangesList"))
      ol <- findOverlaps(feature, se)
      ol_list <- as(ol, "List")
      m <- assay(se, withDimnames = FALSE)
      val <- vapply(ol_list, function(i) {
        mean2(m, idxs = i)
      }, numeric(1))
      data_frame(condition = conditions_df$condition[k],
                 gene = names(feature),
                 nC = countOverlaps(feature, se),
                 mC = val,
                 feature = n)
    }))
  }

  message("Extracting promoters")
  # NOTE: Nothing to do

  message("Extracting first exon")
  i <- ifelse(strand_is_minus, lengths(exons_by_gene), 1L)
  first_exon <- exons_by_gene[as(i, "List")]

  message("Extracting remaining exons")
  ii <- IntegerList(mapply(function(l, ii) {
    setdiff(seq_len(l), ii)
  }, l = lengths(exons_by_gene), ii = i))
  remaining_exons <- exons_by_gene[ii]

  message("Extracting first intron")
  # NOTE: Some genes have 0 introns
  i <- ifelse(lengths(introns_by_gene),
              ifelse(strand_is_minus, lengths(introns_by_gene), 1L),
              NA_integer_)
  first_intron <- introns_by_gene
  first_intron[!is.na(i)] <- introns_by_gene[as(na.omit(i), "List")]
  first_intron[is.na(i)] <- GRangesList(
    lapply(seq_len(sum(is.na(i))), function(i) GRanges()))

  message("Extracting remaining introns")
  ii <- IntegerList(mapply(function(l, ii) {
    setdiff(seq_len(l), ii)
  }, l = lengths(introns_by_gene), ii = i))
  remaining_introns <- introns_by_gene[ii]

  features <- list("promoters" = promoters_by_gene,
                   "first_exon" = first_exon,
                   "remaining_exons" = remaining_exons,
                   "first_intron" = first_intron,
                   "remaining_introns" = remaining_introns)

  bind_rows(mclapply(seq_len(nrow(conditions_df)), function(k) {
    .mCAtGeneFeature(k, features)
  }))
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  mCAtGeneFeatures(
    genes, promoters_by_gene, exons_by_gene, introns_by_gene, CG_BSseq,
    conditions_df) %>%
    mutate(context = "CpG"),
  mCAtGeneFeatures(
    genes, promoters_by_gene, exons_by_gene, introns_by_gene, CA_pos_BSseq,
    conditions_df) %>%
    mutate(context = "CpA (+)"),
  mCAtGeneFeatures(
    genes, promoters_by_gene, exons_by_gene, introns_by_gene, CA_neg_BSseq,
    conditions_df) %>%
    mutate(context = "CpA (-)"),
  mCAtGeneFeatures(
    genes, promoters_by_gene, exons_by_gene, introns_by_gene,
    invertStrand(CA_pos_BSseq),
    conditions_df) %>%
    mutate(context = "CpA (+) (opposite strand)"),
  mCAtGeneFeatures(
    genes, promoters_by_gene, exons_by_gene, introns_by_gene,
    invertStrand(CA_neg_BSseq),
    conditions_df) %>%
    mutate(context = "CpA (-) (opposite strand)"))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "mC_at_gene_features.rds"))
