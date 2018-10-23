# Functions use in integrative analyses
# Peter Hickey
# 2016-11-11


# Overlap `x` with all elements of `features` (a _GRangesList_) to prepare
# data for `UpSetR::upset()`
makeUpSetRList <- function(x, features) {
  feature_names <- names(features)
  setNames(lapply(feature_names, function(feature_name) {
    which(overlapsAny(x, features[[feature_name]]))
  }), feature_names)
}

# Overlap `x` with all elements of `features` (a _GRangesList_) to prepare
# data for `UpSetR::upset()`; overlap computed in bp
makeUpSetRList2 <- function(x, features) {
  feature_names <- names(features)
  setNames(lapply(feature_names, function(feature_name) {
    sum(as.numeric(width(intersect(x, features[[feature_name]]))))
  }), feature_names)
}

# Run `fisher.test()` of `x` against `notx` on all elements of `features` (a
# _GRangesList_) and summarise data as a _data frame_ with odds ratio and 95%
# CI. Adds `offset` to all counts to avoid zero/infinite odds ratio estimates
FT <- function(x, notx, features, db, offset = 1L) {
  x <- lapply(features, function(feature) {
    if (length(feature)) {
      n11 <- sum(overlapsAny(x, feature))
      n12 <- length(x) - n11
      n21 <- sum(overlapsAny(notx, feature))
      n22 <- length(notx) - n21
      m <- matrix(c(n11, n21, n12, n22), ncol = 2) + offset
      ft <- fisher.test(m)
    } else {
      list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
    }
  })

  data.frame(db = db,
             feature = names(features),
             estimate = sapply(x, "[[", "estimate"),
             lower = sapply(lapply(x, "[[", "conf.int"), "[[", 1),
             upper = sapply(lapply(x, "[[", "conf.int"), "[[", 2),
             p.value = sapply(x, "[[", "p.value"),
             stringsAsFactors = FALSE,
             row.names = NULL)
}

# Compute log odds ratio (a la fisher.test()) of bases in `x` vs. bases not in
# `x`over all elements of `features` and sumamrise as a _data frame_.
FT2 <- function(x, features, db, sl) {
  x <- lapply(features, function(feature) {
    if (length(feature)) {
      n11 <- sum(as.numeric(width(GenomicRanges::intersect(x, feature))))
      n.. <- sum(as.numeric(sl))
      n1. <- sum(as.numeric(width(x)))
      n.1 <- sum(as.numeric(width(GenomicRanges::reduce(feature))))
      n12 <- n1. - n11
      n21 <- n.1 - n11
      n22 <- n.. - n11 - n12 - n21
      m <- matrix(c(n11, n21, n12, n22), ncol = 2)
      lor <- log2(m[1, 1]) + log2(m[2, 2]) - log2(m[1, 2]) - log2(m[2, 1])
      lor_se <- sqrt((1 / m[1, 1]) + (1 / m[1, 2]) + (1 / m[2, 1]) +
                       (1 / m[2, 2]))
      lor_lower <- lor - 2 * lor_se
      lor_upper <- lor + 2 * lor_se
      or <- 2 ^ lor
      or_lower <- 2 ^ lor_lower
      or_upper <- 2 ^ lor_upper
      return(list(estimate = or, lower = or_lower, upper = or_upper))
    } else {
      list(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
    }
  })
  data.frame(db = db,
             feature = names(features),
             estimate = sapply(x, "[[", "estimate"),
             lower = sapply(x, "[[", "lower"),
             upper = sapply(x, "[[", "upper"),
             stringsAsFactors = FALSE,
             row.names = NULL)
}

# Compute log odds ratio (a la fisher.test()) of bases in `x` vs. bases in
# `y` over all elements of `features` and sumamrise as a _data frame_.
FT3 <- function(x, y, features, db) {
  x <- lapply(features, function(feature) {
    if (length(feature)) {
      n1. <- sum(as.numeric(width(x)))
      n2. <- sum(as.numeric(width(y)))
      n.. <- n1. + n2.
      n.1 <- sum(as.numeric(width(GenomicRanges::reduce(feature))))
      n.2 <- n.. - n.1
      n11 <- sum(as.numeric(width(GenomicRanges::intersect(x, feature))))
      n12 <- n1. - n11
      n21 <- sum(as.numeric(width(GenomicRanges::intersect(y, feature))))
      n22 <- n2. - n21
      m <- matrix(c(n11, n21, n12, n22), ncol = 2)
      lor <- log2(m[1, 1]) + log2(m[2, 2]) - log2(m[1, 2]) - log2(m[2, 1])
      lor_se <- sqrt((1 / m[1, 1]) + (1 / m[1, 2]) + (1 / m[2, 1]) +
                       (1 / m[2, 2]))
      lor_lower <- lor - 2 * lor_se
      lor_upper <- lor + 2 * lor_se
      or <- 2 ^ lor
      or_lower <- 2 ^ lor_lower
      or_upper <- 2 ^ lor_upper
      return(list(estimate = or, lower = or_lower, upper = or_upper))
    } else {
      list(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
    }
  })
  data.frame(db = db,
             feature = names(features),
             estimate = sapply(x, "[[", "estimate"),
             lower = sapply(x, "[[", "lower"),
             upper = sapply(x, "[[", "upper"),
             stringsAsFactors = FALSE,
             row.names = NULL)
}

# Take a set of regions and classify as being
# "exonic/promoter/UTR-overlapping", "wholly intronic", "wholly intergenic",
# or "other" (weird cases) based on PC genes from GENCODE
genomicContext <- function(query) {
  epu <- overlapsAny(query, c(gencode_features$pc_transcripts$exonic,
                              gencode_features$pc_transcripts$promoter,
                              gencode_features$pc_transcripts$five_utr,
                              gencode_features$pc_transcripts$three_utr))
  intronic <- overlapsAny(query, gencode_features$pc_transcripts$intronic,
                          type = "within") & !epu
  intergenic <- overlapsAny(query, gencode_features$pc_transcripts$intergenic,
                            type = "within") & !epu
  cbind("exon_promoter_utr" = epu,
        "wholly_intronic" = intronic,
        "wholly_intergenic" = intergenic,
        "other" = !epu & !intronic & !intergenic)
}
