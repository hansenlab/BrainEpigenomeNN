# DiffEpi measure around TSS and TES of genes
# Peter Hickey
# 2018-03-08

### ----------------------------------------------------------------------------
### Setup
###

library(GenomicRanges)
library(dplyr)

options("mc.cores" = 20)

extdir <- "../extdata"

### ----------------------------------------------------------------------------
### Load data
###

load(
  file.path(
    extdir, "flow-sorted-brain-gene-level_analyses", "objects",
    "gene-level_analyses_data.rda"))
tss <- resize(genes, width = 1, fix = "start")
tes <- resize(genes, width = 1, fix = "end")

### ----------------------------------------------------------------------------
### Functions
###

# TODO: This uses 1 bp bins; should they be diff(s) sized?
# Returns a tidy data frame
# Each row is the (median) DiffEpi measure for that gene-bin combination
DiffEpiMeasureAroundFixedPoint <- function(fixed_points, feature, yname,
                                           s = c(seq(0, 1000, 1),
                                                 seq(1010, 10000, 10),
                                                 seq(10100, 100000, 100),
                                                 seq(101000, 1000000, 1000))) {
  stopifnot(all(s >= 0))
  # NOTE: Don't compute s == 0 for downstream (computed in upstream to avoid
  #       zero-width ranges)
  s_downstream <- s[s != 0]
  s_upstream <- s
  downstream_vals <- mclapply(s_downstream, function(ss) {
    message("Downstream step = ", ss)
    # NOTE: This can create 'out-of-bounds' ranges. I suppress these warnings
    x <- suppressWarnings(
      resize(
        promoters(fixed_points, downstream = ss, upstream = 0),
        width = 1,
        fix = "end"))
    ol <- findOverlaps(x, feature)
    val <- rep(NA_real_, length(x))
    # NOTE: Use median measurement if there are multiple hits per gene
    tmp <- tapply(
      mcols(feature)[[yname]][subjectHits(ol)],
      queryHits(ol),
      median)
    val[as.numeric(names(tmp))] <- tmp
    val
  })
  downstream_df <- data_frame(
    gene = rep(names(genes), times = length(s_downstream)),
    distance = rep(s_downstream, each = length(fixed_points)),
    DiffEpi_measure = do.call(c, downstream_vals)) %>%
    arrange(distance)

  upstream_vals <- mclapply(s_upstream, function(ss) {
    message("Upstream step = ", ss)
    x <- suppressWarnings(
      resize(
        promoters(fixed_points, downstream = 0, upstream = ss),
        width = 1,
        fix = "start"))
    ol <- findOverlaps(x, feature)
    val <- rep(NA_real_, length(x))
    # NOTE: Use median measurement if there are multiple hits per gene
    tmp <- tapply(
      mcols(feature)[[yname]][subjectHits(ol)],
      queryHits(ol),
      median)
    val[as.numeric(names(tmp))] <- tmp
    val
  })
  upstream_df <- data_frame(
    gene = rep(names(fixed_points), times = length(s_upstream)),
    distance = rep(-s_upstream, each = length(fixed_points)),
    DiffEpi_measure = do.call(c, upstream_vals))

  bind_rows(downstream_df, upstream_df)
}

### ----------------------------------------------------------------------------
### Compute and save results
###

val <- bind_rows(
  DiffEpiMeasureAroundFixedPoint(
    tss, dmrs_NAvsBA9pos, "meanDiff") %>%
    mutate(DiffEpi = "CG-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tss, dmrs_NAvsBA9pos, "meanDiff") %>%
    mutate(DiffEpi = "CG-blocks"),
  DiffEpiMeasureAroundFixedPoint(
    tss, dars_pos, "logFC") %>%
    mutate(DiffEpi = "DARs"),
  DiffEpiMeasureAroundFixedPoint(
    tss, big_dars_pos, "logFC") %>%
    mutate(DiffEpi = "bigDARs"),
  DiffEpiMeasureAroundFixedPoint(
    tss, list_of_CH_DMRs[["mCA (+)"]], "meanDiff") %>%
    mutate(DiffEpi = "CA_pos-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tss, list_of_CH_DMRs[["mCA (-)"]], "meanDiff") %>%
    mutate(DiffEpi = "CA_neg-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tss, invertStrand(list_of_CH_DMRs[["mCA (+)"]]), "meanDiff") %>%
    mutate(DiffEpi = "CA_pos-DMRs (opposite strand)"),
  DiffEpiMeasureAroundFixedPoint(
    tss, invertStrand(list_of_CH_DMRs[["mCA (-)"]]), "meanDiff") %>%
    mutate(DiffEpi = "CA_neg-DMRs (opposite strand)"))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "DiffEpi_measure_around_TSS.rds"))

val <- bind_rows(
  DiffEpiMeasureAroundFixedPoint(
    tes, dmrs_NAvsBA9pos, "meanDiff") %>%
    mutate(DiffEpi = "CG-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tes, dmrs_NAvsBA9pos, "meanDiff") %>%
    mutate(DiffEpi = "CG-blocks"),
  DiffEpiMeasureAroundFixedPoint(
    tes, dars_pos, "logFC") %>%
    mutate(DiffEpi = "DARs"),
  DiffEpiMeasureAroundFixedPoint(
    tes, big_dars_pos, "logFC") %>%
    mutate(DiffEpi = "bigDARs"),
  DiffEpiMeasureAroundFixedPoint(
    tes, list_of_CH_DMRs[["mCA (+)"]], "meanDiff") %>%
    mutate(DiffEpi = "CA_pos-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tes, list_of_CH_DMRs[["mCA (-)"]], "meanDiff") %>%
    mutate(DiffEpi = "CA_neg-DMRs"),
  DiffEpiMeasureAroundFixedPoint(
    tes, invertStrand(list_of_CH_DMRs[["mCA (+)"]]), "meanDiff") %>%
    mutate(DiffEpi = "CA_pos-DMRs (opposite strand)"),
  DiffEpiMeasureAroundFixedPoint(
    tes, invertStrand(list_of_CH_DMRs[["mCA (-)"]]), "meanDiff") %>%
    mutate(DiffEpi = "CA_neg-DMRs (opposite strand)"))

saveRDS(
  val,
  file.path(extdir, "flow-sorted-brain-gene-level_analyses", "objects",
            "DiffEpi_measure_around_TES.rds"))
