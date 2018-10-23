# DiffEpi measure upstream and downstream of gene body
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
# Each row is the (median) DiffEpi measure for that gene-bin combination
DiffEpiMeasureUpstreamAndDownstreamOfGene <- function(genes, feature, yname,
                                                      binsize = 10 ^ 3,
                                                      max_distance = 10 ^ 6) {
  message("Upstream of genes")
  # NOTE: This created out-of-bounds ranges, but this doesn't matter
  upstream <- suppressWarnings(flank(genes, max_distance, start = TRUE))
  tiled_upstream <- tile(upstream, width = binsize)
  names(tiled_upstream) <- names(genes)
  # NOTE: Reverse order of bins for negative strand genes so that bins go in
  #       ascending order from upstream to downstream
  tiled_upstream[strand(upstream) == "-"] <-
    endoapply(tiled_upstream[strand(upstream) == "-"], rev)
  unlisted_tiled_upstream <- unlist(tiled_upstream)
  message("\tFinding overlaps")
  ol <- findOverlaps(unlisted_tiled_upstream, feature)
  val <- rep(NA_real_, length(unlisted_tiled_upstream))
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("\tConstructing result")
  upstream_df <- data_frame(
    bin = rep(-rev(seq_len(max_distance / binsize)),
              times = length(tiled_upstream)),
    gene = rep(names(tiled_upstream), each = max_distance / binsize),
    DiffEpi_measure = val)

  message("Downstream of genes")
  # NOTE: This created out-of-bounds ranges, but this doesn't matter
  downstream <- suppressWarnings(flank(genes, max_distance, start = FALSE))
  tiled_downstream <- tile(downstream, width = binsize)
  names(tiled_downstream) <- names(genes)
  # NOTE: Reverse order of bins for negative strand genes so that bins go in
  #       ascending order from upstream to downstream
  tiled_downstream[strand(downstream) == "-"] <-
    endoapply(tiled_downstream[strand(downstream) == "-"], rev)
  unlisted_tiled_downstream <- unlist(tiled_downstream)
  message("\tFinding overlaps")
  ol <- findOverlaps(unlisted_tiled_downstream, feature)
  val <- rep(NA_real_, length(unlisted_tiled_downstream))
  # NOTE: Use median measurement if there are multiple hits per bin
  tmp <- tapply(
    mcols(feature)[[yname]][subjectHits(ol)],
    queryHits(ol),
    median)
  val[as.numeric(names(tmp))] <- tmp
  message("\tConstructing result")
  downstream_df <- data_frame(
    bin = rep(seq_len(max_distance / binsize),
              times = length(tiled_downstream)),
    gene = rep(names(tiled_downstream), each = max_distance / binsize),
    DiffEpi_measure = val)

  bind_rows(upstream_df, downstream_df)
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
  DiffEpiMeasureUpstreamAndDownstreamOfGene(
    genes, feature, yname) %>%
    mutate(DiffEpi = diff_epi)
}, feature = features, yname = ynames, diff_epi = diff_epis, SIMPLIFY = FALSE))

saveRDS(
  val,
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "DiffEpi_measure_upstream_and_downstream_of_gene.rds"))
