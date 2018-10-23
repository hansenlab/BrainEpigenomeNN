# DiffEpi measure at promoter, first exon, first intron, rest of exons,
# rest of introns
# Peter Hickey
# 2018-03-08

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)

options("mc.cores" = 8)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))

### ----------------------------------------------------------------------------
### Functions
###

# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene feature
DiffEpiMeasureAtGeneFeatures <- function(genes, promoters_by_gene,
                                         exons_by_gene, introns_by_gene,
                                         feature, yname) {
  # Ensure genes, exons, and introns are all ordered the same way and that
  # exons and introns within each gene are sorted
  message("Sorting genes, promoters, exons, and introns")
  o <- names(genes)
  promoters_by_gene <- sort(promoters_by_gene[o])
  exons_by_gene <- sort(exons_by_gene[o])
  introns_by_gene <- sort(introns_by_gene[o])
  strand_is_minus <- strand(genes) == "-"

  message("Extracting promoters")
  message("Finding overlaps")
  ol <- findOverlaps(promoters_by_gene, feature)
  val <- rep(NA_real_, length(promoters_by_gene))
  message("Computing median DiffEpi_measure")
  # TODO: Should I take all DiffEpi over promoters?
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  promoters_df <- data_frame(
    feature = "promoters",
    gene = names(promoters_by_gene),
    DiffEpi_measure = val)

  message("Extracting first exon")
  i <- ifelse(strand_is_minus, lengths(exons_by_gene), 1L)
  first_exon <- exons_by_gene[as(i, "List")]
  message("Finding overlaps")
  ol <- findOverlaps(first_exon, feature)
  val <- rep(NA_real_, length(first_exon))
  message("Computing median DiffEpi_measure")
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  first_exon_df <- data_frame(
    feature = "first_exon",
    gene = names(first_exon),
    DiffEpi_measure = val)

  message("Extracting remaining exons")
  ii <- IntegerList(mapply(function(l, ii) {
    setdiff(seq_len(l), ii)
  }, l = lengths(exons_by_gene), ii = i))
  remaining_exons <- exons_by_gene[ii]
  message("Finding overlaps")
  ol <- findOverlaps(remaining_exons, feature)
  val <- rep(NA_real_, length(remaining_exons))
  message("Computing median DiffEpi_measure")
  # TODO: Should I take all DiffEpi over remaining exons?
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  remaining_exons_df <- data_frame(
    feature = "remaining_exons",
    gene = names(remaining_exons),
    DiffEpi_measure = val)

  message("Extracting first intron")
  # NOTE: Some genes have 0 introns
  i <- ifelse(lengths(introns_by_gene),
              ifelse(strand_is_minus, lengths(introns_by_gene), 1L),
              NA_integer_)
  first_intron <- introns_by_gene
  first_intron[!is.na(i)] <- introns_by_gene[as(na.omit(i), "List")]
  first_intron[is.na(i)] <- GRangesList(
    lapply(seq_len(sum(is.na(i))), function(i) GRanges()))
  message("Finding overlaps")
  ol <- findOverlaps(first_intron, feature)
  val <- rep(NA_real_, length(first_intron))
  message("Computing median DiffEpi_measure")
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  first_intron_df <- data_frame(
    feature = "first_intron",
    gene = names(first_intron),
    DiffEpi_measure = val)

  message("Extracting remaining introns")
  ii <- IntegerList(mapply(function(l, ii) {
    setdiff(seq_len(l), ii)
  }, l = lengths(introns_by_gene), ii = i))
  remaining_introns <- introns_by_gene[ii]
  message("Finding overlaps")
  ol <- findOverlaps(remaining_introns, feature)
  val <- rep(NA_real_, length(remaining_introns))
  message("Computing median DiffEpi_measure")
  # TODO: Should I take all DiffEpi over remaining introns?
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("Constructing result")
  remaining_introns_df <- data_frame(
    feature = "remaining_introns",
    gene = names(remaining_introns),
    DiffEpi_measure = val)

  bind_rows(promoters_df, first_exon_df, first_intron_df,
            remaining_exons_df, remaining_introns_df)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

features <- list(dmrs_NAvsBA9pos, blocks_pos, dars_pos, big_dars_pos,
                 list_of_CH_DMRs[["mCA (+)"]],
                 list_of_CH_DMRs[["mCA (-)"]],
                 invertStrand(list_of_CH_DMRs[["mCA (+)"]]),
                 invertStrand(list_of_CH_DMRs[["mCA (-)"]]))
ynames <- list("meanDiff", "meanDiff", "logFC", "logFC", "meanDiff",
               "meanDiff", "meanDiff", "meanDiff")
diff_epis <- list("CG-DMRs", "CG-blocks", "DARs", "bigDARs", "CA_pos-DMRs",
                  "CA_neg-DMRs", "CA_pos-DMRs (opposite strand)",
                  "CA_neg-DMRs (opposite strand)")

val <- bind_rows(mcmapply(function(feature, yname, diff_epi) {
  DiffEpiMeasureAtGeneFeatures(
    genes, promoters_by_gene,exons_by_gene, introns_by_gene, feature,
    yname) %>%
    mutate(DiffEpi = diff_epi)
}, feature = features, yname = ynames, diff_epi = diff_epis, SIMPLIFY = FALSE))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "DiffEpi_measure_at_gene_features.rds"))
