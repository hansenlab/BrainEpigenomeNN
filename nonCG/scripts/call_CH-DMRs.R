# Smooth non-CG data
# Peter Hickey
# 2017-07-04

# NOTE: Run with qrsh -pe local 6 -l mem_free=50G,h_vmem=51G

### ============================================================================
### Setup
###

library(bsseq)

options("DelayedArray.block.size" = 45000000)
options("mc.cores" = 6)

# ------------------------------------------------------------------------------
# Load BSseq objects
#

pos_CA_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/pos_CA_small-flow-sorted-brain-wgbs")
neg_CA_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/neg_CA_small-flow-sorted-brain-wgbs")
pos_CT_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/pos_CT_small-flow-sorted-brain-wgbs")
neg_CT_BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/neg_CT_small-flow-sorted-brain-wgbs")

# ------------------------------------------------------------------------------
# Functions
#

# NOTE: These functions are adapted from PeteHaitch/bsseq@HDF5Array. I'll
#       intergrate with the devel branch of bsseq when I get back from holiday

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/lmFit.R#L125
.lm.series <- function(M, design = NULL) {
  narrays <- nrow(M)

  # NOTE: Redundant if called from .lmFit()
  if (is.null(design)) {
    design <- matrix(1, narrays, 1)
  } else {
    design <- as.matrix(design)
  }

  nbeta <- ncol(design)
  coef.names <- colnames(design)
  if (is.null(coef.names)) {
    coef.names <- paste("x", 1:nbeta, sep = "")
  }

  ngenes <- ncol(M)
  # NOTE: Initisalise with NA_real_ rather than NA (which is logical) to
  #       ensure the storage mode of these matrices is numeric
  stdev.unscaled <- beta <- matrix(NA_real_, ngenes, nbeta,
                                   dimnames = list(colnames(M), coef.names))

  # Check whether QR-decomposition is constant for all genes
  # If so, fit all genes in one sweep
  NoProbeWts <- all(is.finite(M))
  if (NoProbeWts) {
    # NOTE: If wanting to be really crazy, could use .lm.fit(). This would
    #       require some post-processing of the .lm.fit() output and it is
    #       not obvious that it is worth the effort involved
    fit <- lm.fit(design, M)
    if (fit$df.residual > 0) {
      if (is.matrix(fit$effects)) {
        fit$sigma <- sqrt(colMeans(
          fit$effects[(fit$rank + 1):narrays, , drop = FALSE]^2))
      } else {
        fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
      }
    } else {
      fit$sigma <- rep(NA, ngenes)
    }

    # NOTE: In bsseq we only need `sigma`, `cov.coefficients`,
    #       `coefficients`, `stdev.unscaled`
    fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    fit$coefficients <- t(fit$coefficients)
    fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    est <- fit$qr$pivot[1:fit$qr$rank]
    dimnames(fit$cov.coefficients) <- list(coef.names[est], coef.names[est])
    stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)),
                                    ngenes, fit$qr$rank, byrow = TRUE)
    fit$stdev.unscaled <- stdev.unscaled
    dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    fit$pivot <- fit$qr$pivot
    return(fit)
  } else {
    stop("'M' contains non-finite values")
  }
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/lmFit.R#L58
.lmFit <- function(y, design) {
  stopifnot(is.matrix(y))
  if (is.null(design)) {
    design <- matrix(1, ncol(y), 1)
  } else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") {
      stop("design must be a numeric matrix")
    }
    if (nrow(design) != nrow(y)) {
      stop("row dimension of design doesn't match column dimension of ",
           "data object")
    }
  }
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    cat("Coefficients not estimable:", paste(ne, collapse = " "), "\n")
  }

  fit <- .lm.series(M = y, design = design)

  if (NCOL(fit$coef) > 1) {
    n <- rowSums(is.na(fit$coef))
    n <- sum(n > 0 & n < NCOL(fit$coef))
    if (n > 0) {
      warning("Partial NA coefficients for ", n, " probe(s)",
              call. = FALSE)
    }
  }

  fit$genes <- NULL
  fit$method <- "ls"
  fit$design <- design
  new("MArrayLM", fit)
}

# NOTE: y should be transposed from its normal orientation, i.e., samples as
#       rows and loci as columns
.bsseq.lm.fit <- function(y, design, contrasts) {
  fit <- .lmFit(y, design)
  limma::contrasts.fit(fit, contrasts)
}

.BSmooth.fstat <- function(tAllPs, parameters, design, contrasts,
                           verbose = TRUE) {
  ptime1 <- proc.time()
  fitC <- .bsseq.lm.fit(tAllPs, design, contrasts)
  ## Need
  ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
  ## actuall just need
  ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
  ##   rawSds <- fitC$sigma
  ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
  rawSds <- as.matrix(fitC$sigma)
  cor.coefficients <- cov2cor(fitC$cov.coefficients)
  rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
  names(dimnames(rawTstats)) <- NULL
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose) {
    cat(sprintf("done in %.1f sec\n", stime))
  }
  list(rawTstats = rawTstats,
       rawSds = rawSds,
       cor.coefficients = cor.coefficients)
}

BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE,
                          hdf5 = FALSE) {
  stopifnot(is(BSseq, "BSseq"))
  stopifnot(hasBeenSmoothed(BSseq))

  if (verbose) {
    cat("[BSmooth.fstat] fitting linear models ... ")
  }
  tAllPs <- as.matrix(t(getMeth(BSseq, type = "smooth", what = "perBase",
                                confint = FALSE)))
  stats <- .BSmooth.fstat(tAllPs = tAllPs,
                          design = design,
                          contrasts = contrasts,
                          verbose = verbose)
  parameters <- c(getBSseq(BSseq, "parameters"),
                  list(design = design, contrasts = contrasts))

  BSseqStat(gr = granges(BSseq), stats = stats, parameters = parameters)
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/utils.R#L64
.smoothSd <- function(Sds, k, qSd) {
  k0 <- floor(k / 2)
  if (all(is.na(Sds))) {
    return(Sds)
  }
  thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
  addSD <- rep(median(Sds, na.rm = TRUE), k0)
  sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
  sSds
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/BSmooth.fstat.R#L87
.smoothSds <- function(clusterIdx, rawSds, k, qSd, mc.cores = 1) {
  smoothSds <- do.call("c",
                       mclapply(clusterIdx, function(idx) {
                         # NOTE: Need to realise rawSds as an array because
                         #       .smoothSd() works with in-memory data
                         rawSds <- as.array(rawSds[idx, ])
                         .smoothSd(rawSds, k = k, qSd = qSd)
                       }, mc.cores = mc.cores))
  smoothSds <- matrix(smoothSds, ncol = 1)
  smoothSds
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/BSmooth.fstat.R#L137
.computeStat <- function(rawTstats, rawSds, smoothSds, coef, cor.coefficients) {
  # TODO: Do I really need to explicitly subset rawSds and smoothSds
  tstats <- rawTstats[, coef, drop = FALSE] * rawSds[, 1L] / smoothSds[, 1L]
  # TODO: Need to realise tstats in case components are DelayedArray objects
  if (length(coef) > 1) {
    # TODO: Need to check this branch
    cor.coefficients <- cor.coefficients[coef, coef]
    stat <- matrix(as.numeric(limma::classifyTestsF(tstats,
                                                    cor.coefficients,
                                                    fstat.only = TRUE)),
                   ncol = 1)
  } else {
    stat <- tstats
  }
  stat
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/dmrFinder.R#L2
.dmrFinder <- function(dmrStat, cutoff, seqnames_as_char, positions, maxGap,
                       verbose = FALSE) {
  if (length(cutoff) == 1) {
    cutoff <- c(-cutoff, cutoff)
  }
  direction_logical <- dmrStat >= cutoff[2]
  direction <- as.integer(direction_logical)
  idx <- dmrStat <= cutoff[1]
  direction[idx] <- -1L
  direction[is.na(direction)] <- 0L
  regions <- bsseq:::regionFinder3(x = direction,
                                   chr = seqnames_as_char,
                                   positions = positions,
                                   maxGap = maxGap,
                                   verbose = verbose)
  if (is.null(regions$down) && is.null(regions$up)) {
    return(NULL)
  }
  if (verbose) {
    cat("[dmrFinder] creating dmr data.frame\n")
  }
  regions <- do.call(rbind, regions)
  rownames(regions) <- NULL
  regions[["width"]] <- regions[["end"]] - regions[["start"]] + 1L
  regions[["invdensity"]] <- regions[["width"]] / regions[["n"]]
  regions
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/getStats.R#L14
.getRegionStats <- function(stat) {
  areaStat <- sum(stat)
  maxStat <- max(stat)
  c(areaStat, maxStat)
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/getStats.R#L21
.getRegionStats_BSseqStat <- function(ov, stat) {
  # NOTE: Rather than split(ov) [as in the original implementation of
  #       getStats_BSseqStat()], we split() the required element of the
  #       the `stat` matrix-like object. This is slightly more efficient if
  #       stat is a matrix and **much** more efficient if stat is a
  #       DelayedArray.
  # NOTE: split,DelayedArray-method returns a *List and thus realises the
  #       data in memory. And, of course, split.array() returns a list, which
  #       is already realised in memory
  stat <- stat[queryHits(ov), ]
  regionStats <- matrix(NA, ncol = 2, nrow = subjectLength(ov),
                        dimnames = list(NULL, c("areaStat", "maxStat")))
  tmp <- lapply(split(stat, subjectHits(ov)), .getRegionStats)
  regionStats[as.integer(names(tmp)), ] <- do.call(rbind, tmp)
  regionStats
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/BSmooth.fstat.R#L240
.fstat.dmr.pipeline <- function(gr, tAllPs, parameters, design, contrasts,
                                clusterIdx, coef, cutoff, maxGap, k = 101,
                                qSd = 0.75, return_bstat = TRUE,
                                verbose = TRUE, hdf5 = FALSE) {
  stats <- .BSmooth.fstat(tAllPs = tAllPs,
                          design = design,
                          contrasts = contrasts)
  smoothSds <- .smoothSds(clusterIdx = clusterIdx,
                          rawSds = stats[["rawSds"]],
                          k = k,
                          qSd = qSd,
                          mc.cores = 1)
  stat <- .computeStat(rawTstats = stats[["rawTstats"]],
                       rawSds = stats[["rawSds"]],
                       smoothSds = smoothSds,
                       coef = coef,
                       cor.coefficients = stats[["cor.coefficients"]])
  dmrs <- .dmrFinder(dmrStat = stat,
                     cutoff = cutoff,
                     seqnames_as_char = as.character(seqnames(gr)),
                     positions = start(gr),
                     maxGap = maxGap,
                     verbose = max(as.integer(verbose) - 1L, 0L))
  if (!is.null(dmrs)) {
    ov <- findOverlaps(gr, data.frame2GRanges(dmrs))
    dmrs_stats <- .getRegionStats_BSseqStat(ov = ov, stat = stat)
    dmrs <- cbind(dmrs, dmrs_stats)
    dmrs <- dmrs[order(abs(dmrs[["areaStat"]]), decreasing = TRUE), ]
  }
  if (return_bstat) {
    parameters <- c(parameters,
                    list(design = design, contrasts = contrasts))
    bstat <- BSseqStat(gr = gr,
                       stats = c(lapply(c(stats, list(smoothSds = smoothSds),
                                        list(stat = stat)), DelayedArray),
                                 list(stat.type = "fstat")),
                       parameters = parameters)
  } else {
    bstat <- NULL
  }
  list(bstat = bstat, dmrs = dmrs)
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/permutations.R#L28
.getNullDistribution_BSmooth.fstat <- function(permutationMatrix, gr, tAllPs,
                                               design, contrasts, clusterIdx,
                                               coef, cutoff, maxGap, mc.cores,
                                               verbose = TRUE,
                                               # NOTE: These args are specific
                                               # to callDMRs()
                                               strand,
                                               context,
                                               size) {

  message(paste0("[getNullDistribution_BSmooth.fstat] performing ",
                 nrow(permutationMatrix), " permutations\n"))

  nullDist <- mclapply(seq_len(nrow(permutationMatrix)), function(ii) {
    ptime1 <- proc.time()
    # NOTE: More efficient to permute design matrix using
    #       permutationMatrix[ii, ] than to permute the raw data with
    #       tAllPs[permutationMatrix[ii, ]]
    permuted_design <- design[permutationMatrix[ii, ], , drop = FALSE]
    bstat_and_dmrs <- .fstat.dmr.pipeline(gr = gr,
                                          tAllPs = tAllPs,
                                          design = permuted_design,
                                          contrasts = contrasts,
                                          clusterIdx = clusterIdx,
                                          coef = coef,
                                          cutoff = cutoff,
                                          maxGap = maxGap,
                                          return_bstat = FALSE,
                                          verbose = verbose,
                                          hdf5 = hdf5)
    dmrs <- bstat_and_dmrs[["dmrs"]]
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
      message(sprintf("[getNullDistribution_BSmooth.fstat] completing permutation %d in %.1f sec\n", ii, stime))
    }
    saveRDS(dmrs,
            paste0(
              "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/null_DMRs/",
              strand, "_", context, "_", size, ".permutation_", ii,
              ".dmrs.rds"))
    dmrs
  })
  nullDist
}

# NOTE: https://github.com/PeteHaitch/bsseq/blob/HDF5Array/R/BSmooth.fstat.R#L284
fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                           type = "dmrs", mc.cores = 1, verbose = TRUE,
                           context, strand, size) {
  type <- match.arg(type, c("dmrs", "blocks"))
  stopifnot(is(BSseq, "BSseq"))
  stopifnot(hasBeenSmoothed(BSseq))

  permutationMatrix <- bsseq:::permuteAll(nperm, design)
  if (nrow(permutationMatrix) < nperm) {
    warning(paste0("Only ", nrow(permutationMatrix), " unique ",
                   "permutations exist (requested ", nperm),
            " permutations)")
    nperm <- nrow(permutationMatrix)
  }
  saveRDS(permutationMatrix,
          paste0(
            "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/null_DMRs/",
            strand, "_", context, "_", size, ".permutationMatrix.rds"))

  # NOTE: Certain objects are used in identifying the candidate DMRs and can
  #       be reused without changed when identifying DMRs in the permuted
  #       data. Constructing these once saves unnecessary computation
  gr <- rowRanges(BSseq)
  # TODO: Should tAllPs be realised at this point; i.e. can the forked
  #       processed share the same tAllPs object?
  tAllPs <- as.matrix(t(getMeth(BSseq, type = "smooth", what = "perBase",
                                confint = FALSE)))
  parameters <- getBSseq(BSseq, "parameters")
  clusterIdx <- bsseq:::makeClusters(gr, maxGap = maxGap.sd)
  if (is.null(coef)) {
    coef <- seq_len(ncol(design) - 1L)
  }

  message("Finding DMRs")
  bstat_and_dmrs <- .fstat.dmr.pipeline(gr = gr,
                                        tAllPs = tAllPs,
                                        parameters = parameters,
                                        design = design,
                                        contrasts = contrasts,
                                        clusterIdx = clusterIdx,
                                        coef = coef,
                                        cutoff = cutoff,
                                        maxGap = maxGap.dmr,
                                        return_bstat = TRUE,
                                        verbose = verbose)
  saveRDS(bstat_and_dmrs,
          paste0(
            "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/null_DMRs/",
            strand, "_", context, "_", size, ".bstat_and_dmrs.rds"))

  bstat <- bstat_and_dmrs[["bstat"]]
  dmrs <- bstat_and_dmrs[["dmrs"]]

  if (is.null(dmrs)) {
    stop("No DMRs identified. Consider reducing the 'cutoff' from (",
         paste0(cutoff, collapse = ", "), ")")
  }

  message("Running permutations")
  nullDist <- .getNullDistribution_BSmooth.fstat(
    permutationMatrix = permutationMatrix,
    gr = gr,
    tAllPs = tAllPs,
    design = design,
    contrasts = contrasts,
    clusterIdx = clusterIdx,
    coef = coef,
    cutoff = cutoff,
    maxGap = maxGap.dmr,
    mc.cores = mc.cores,
    verbose = verbose,
    strand = strand,
    context = context,
    size = size)

  # Compute FWER for candidate DMRs
  message("Computing FWER")
  fwer <- bsseq:::getFWER.fstat(null = c(list(dmrs), nullDist), type = type)
  dmrs$fwer <- fwer

  # Compute average methylation level in each group (`fac`) for each
  # candidate DMR
  message("Computing condition-average methylation in DMRs")
  # TODO: Could use tAllPs
  meth <- as.matrix(getMeth(BSseq, regions = dmrs, what = "perRegion"))
  meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
  dmrs <- cbind(dmrs, meth)
  dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)

  # Return BSseqStat object computed using unpermuted data, DMRs with FWER,
  # the permutation matrix, and the DMRs in the permuted data
  list(bstat = bstat, dmrs = dmrs, permutationMatrix = permutationMatrix,
       nullDist = nullDist)
}

callDMRs <- function(strand,
                     context,
                     BSseq,
                     size,
                     cutoff,
                     nperm = 1000,
                     coef = NULL,
                     maxGap.sd = 10 ^ 8,
                     maxGap.dmr = 300,
                     type = "dmrs",
                     mc.cores = 1,
                     verbose = TRUE) {
  # Subset to NeuN+ samples
  BSseq <- BSseq[, BSseq$NeuN == "pos"]

  # Construct design matrix and contrasts
  design <- model.matrix(~ BSseq$Tissue)
  colnames(design) <- gsub("BSseq\\$", "", colnames(design))
  contrasts <- diag(rep(1, ncol(design)))[, -1]
  rownames(contrasts) <- colnames(design)
  fac <- colData(BSseq)$Tissue

  message(strand)

  fstat.pipeline(BSseq = BSseq,
                 design = design,
                 contrasts = contrasts,
                 cutoff = cutoff,
                 fac = fac,
                 nperm = nperm,
                 coef = coef,
                 maxGap.sd = maxGap.sd,
                 maxGap.dmr = maxGap.dmr,
                 type = type,
                 mc.cores = mc.cores,
                 verbose = verbose,
                 strand = strand,
                 context = context,
                 size = size)
}

### ============================================================================
### mCA
###

# ------------------------------------------------------------------------------
# Pos strand (small smooth)
#

set.seed(666)
pos_mCA_val <- callDMRs(strand = "pos",
                        context = "CA",
                        BSseq = pos_CA_BSseq,
                        size = "small",
                        cutoff = 16,
                        nperm = 1000,
                        coef = NULL,
                        maxGap.sd = 10 ^ 8,
                        maxGap.dmr = 300,
                        type = "dmrs",
                        mc.cores = getOption("mc.cores"))

# ------------------------------------------------------------------------------
# Neg strand (small smooth)
#

set.seed(3087)
neg_mCA_val <- callDMRs(strand = "neg",
                        context = "CA",
                        BSseq = neg_CA_BSseq,
                        size = "small",
                        cutoff = 16,
                        nperm = 1000,
                        coef = NULL,
                        maxGap.sd = 10 ^ 8,
                        maxGap.dmr = 300,
                        type = "dmrs",
                        mc.cores = getOption("mc.cores"))

### ============================================================================
### mCT
###

# ------------------------------------------------------------------------------
# Pos strand (small smooth)
#

set.seed(21202)
pos_mCT_val <- callDMRs(strand = "pos",
                        context = "CT",
                        BSseq = pos_CT_BSseq,
                        size = "small",
                        cutoff = 16,
                        nperm = 1000,
                        coef = NULL,
                        maxGap.sd = 10 ^ 8,
                        maxGap.dmr = 300,
                        type = "dmrs",
                        mc.cores = getOption("mc.cores"))

# ------------------------------------------------------------------------------
# Neg strand (small smooth)
#

set.seed(1987)
neg_mCT_val <- callDMRs(strand = "neg",
                        context = "CT",
                        BSseq = neg_CT_BSseq,
                        size = "small",
                        cutoff = 16,
                        nperm = 1000,
                        coef = NULL,
                        maxGap.sd = 10 ^ 8,
                        maxGap.dmr = 300,
                        type = "dmrs",
                        mc.cores = getOption("mc.cores"))
