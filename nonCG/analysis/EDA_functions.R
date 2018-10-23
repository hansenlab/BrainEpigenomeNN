# Functions used in EDA of mCH
# Peter Hickey
# 2017-04-27

#===============================================================================
# Functions that calculate stuff
#

# NOTE: corrected t-stat not readily applicable for non-CG methylation, so not
#       using it for now
tstatPipeline <- function(bsseq, group1, group2, name, estimate.var = "paired",
                          stat = "tstat", cutoff = c(-4.6, 4.6), maxGap = 300,
                          local.correct = TRUE, permutations = FALSE,
                          max_perm = 1000, mc.cores = 1) {
  stat <- match.arg(stat)
  bstat <- BSmooth.tstat(BSseq = bsseq,
                         group1 = group1,
                         group2 = group2,
                         estimate.var = estimate.var,
                         local.correct = local.correct)
  dmrs <- dmrFinder(bstat = bstat,
                    cutoff = cutoff,
                    maxGap = maxGap,
                    stat = stat)

  # Permutations
  if (permutations) {
    if (!is.null(dmrs)) {
      idxMatrix <- bsseq:::makeIdxMatrix(group1, group2)
      # NOTE: The first element of idxMatrix is just group1 and group2 unpermuted
      idxMatrix <- lapply(idxMatrix, function(x) x[-1, , drop = FALSE])
      nperm <- min(max_perm, nrow(idxMatrix[[1]]))
      idxMatrix <- lapply(idxMatrix,
                          function(x) x[seq_len(nperm), , drop = FALSE])
      nullDist <- bsseq:::getNullDistribution_BSmooth.tstat(
        BSseq = bsseq,
        idxMatrix1 = idxMatrix[[1]],
        idxMatrix2 = idxMatrix[[2]],
        estimate.var = estimate.var,
        local.correct = local.correct,
        stat = stat,
        maxGap = maxGap,
        cutoff = cutoff,
        mc.cores = mc.cores)
      fwer <- bsseq:::getFWER(null = c(list(dmrs), nullDist), type = "dmrs")
      dmrs$fwer <- fwer
      dmrs$nperm <- nperm
    }
  }
  list(bstat = bstat, dmrs = dmrs)
}

#===============================================================================
# Functions that plot stuff
#

# TODO: Modify ylim for beta-values (and perhaps modify ylim for bstat)
# NOTE: A real hack to load the plotting functions from bsseq into the global
#       environment
source("https://raw.githubusercontent.com/kasperdanielhansen/bsseq/47dcc0f2b181145210c8e2603142b2fc79dd9345/R/plotting.R")

#-------------------------------------------------------------------------------
# plotGeneTrack(): gets the `addNames` option for whether gene names
# are printed
#

plotGeneTrack <- function(gr, geneTrack, cex, addNames = FALSE) {
  geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
  ol <- findOverlaps(geneTrack_gr, gr)
  genes <- geneTrack[queryHits(ol), ]
  plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
       ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)),
       xlab = "", ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
  if (nrow(genes) > 0) {
    for (g in 1:nrow(genes)) {
      geneind2 = which(geneTrack$gene_name == genes$gene_name[g])
      geneind2 = geneind2[which(geneTrack$isoforms[geneind2] == 1)]
      direction = unique(geneTrack$strand[geneind2])
      ES = geneTrack$start[geneind2]
      EE = geneTrack$end[geneind2]
      Exons = cbind(ES, EE)
      if (direction == "+") {
        lines(x = c(min(ES), max(EE)),
              y = c(0.65, 0.65))
        apply(Exons, 1, function(x) {
          polygon(c(x[1], x[2], x[2], x[1]),
                  c(0.45, 0.45, 0.85, 0.85), col = "darkgrey")
        })
        if (addNames) {
          text((max(start(gr), min(ES)) +
                  min(end(gr), max(EE))) / 2, 1.2,
               genes$gene_name[g], cex = cex)
        }
      } else {
        lines(x = c(min(ES), max(EE)),
              y = c(-0.65, -0.65))
        apply(Exons, 1, function(x)
          polygon(c(x[1], x[2], x[2], x[1]),
                  c(-0.45, -0.45, -0.85, -0.85),
                  col = "darkgrey"))
        if (addNames) {
          text((max(start(gr), min(ES)
          ) + min(end(gr), max(EE)
          )) / 2, -1.2, genes$gene_name[g], cex = cex)
        }
      }

    }
  }
}

.bsHighlightRegions <- function(regions, gr, ylim, regionCol, highlightMain) {
  if(is.data.frame(regions))
    regions <- data.frame2GRanges(regions)
  if(highlightMain)
    regions <- c(regions, gr, ignore.mcols = TRUE)
  if(is.null(regions)) return(NULL)
  ## regions <- pintersect(region, rep(gr, length(regions)))
  ## regions <- regions[width(regions) == 0]
  regions <- subsetByOverlaps(regions, gr)
  regions <- pintersect(regions, rep(gr, length(regions)))
  if(length(regions) == 0)
    return(NULL)
  rect(xleft = start(regions), xright = end(regions), ybottom = ylim[1],
       ytop = ylim[2], col = regionCol, border = NA)
}

#-------------------------------------------------------------------------------
# Helper functions for plotRegion2()
#

# Modified version of bsseq:::.plotSmoothData() to allow passing of ylab
.plotSmoothData <- function(BSseq, region, extend, addRegions, col, lty, lwd,
                            regionCol, addTicks, addPoints, pointsMinCov,
                            highlightMain, ylab = "Methylation") {
  gr <- .bsGetGr(BSseq, region, extend)
  BSseq <- subsetByOverlaps(BSseq, gr)

  ## Extract basic information
  sampleNames <- sampleNames(BSseq)
  names(sampleNames) <- sampleNames
  positions <- start(BSseq)
  smoothPs <- getMeth(BSseq, type = "smooth")
  rawPs <- getMeth(BSseq, type = "raw")
  coverage <- getCoverage(BSseq)

  # Realise in memory data that are to be plotted
  if (addPoints) {
    rawPs <- as.array(rawPs)
    coverage <- as.array(coverage)
  }
  smoothPs <- as.array(smoothPs)

  ## get col, lwd, lty
  colEtc <- .bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  plot(positions[1],
       0.5,
       type = "n",
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1),
       xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = ylab)
  axis(side = 2, at = c(0.2, 0.5, 0.8))
  if (addTicks) {
    rug(positions)
  }

  .bsHighlightRegions(regions = addRegions,
                      gr = gr,
                      ylim = c(0,1),
                      regionCol = regionCol,
                      highlightMain = highlightMain)

  if (addPoints) {
    sapply(1:ncol(BSseq), function(sampIdx) {
      abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
    })
  } # This adds vertical grey lines so we can see where points are plotted

  sapply(1:ncol(BSseq), function(sampIdx) {
    .bsPlotLines(positions, smoothPs[, sampIdx], col = colEtc$col[sampIdx],
                 lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx],
                 plotRange = c(start(gr), end(gr)))
  })

  if (addPoints) {
    sapply(1:ncol(BSseq), function(sampIdx) {
      .bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx],
                    col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov)
    })
  }
}

.ATACGetCol <- function(object, col, lty, lwd) {
  if (is.null(col)) {
    if ("col" %in% names(colData(object))) {
      col <- colData(object)[["col"]]
    } else {
      col <- rep("black", nrow(colData(object)))
    }
  }
  if (length(col) != ncol(object)) {
    col <- rep(col, length.out = ncol(object))
  }
  if (is.null(names(col))) {
    names(col) <- colnames(object)
  }
  if (is.null(lty)) {
    if ("lty" %in% names(colData(object))) {
      lty <- colData(object)[["lty"]]
    } else {
      lty <- rep(1, ncol(object))
    }
  }
  if (length(lty) != ncol(object)) {
    lty <- rep(lty, length.out = ncol(object))
  }
  if (is.null(names(lty))) {
    names(lty) <- colnames(object)
  }
  if (is.null(lwd)) {
    if ("lwd" %in% names(colData(object))) {
      lwd <- colData(object)[["lwd"]]
    } else {
      lwd <- rep(1, nrow(colData(object)))
    }
  }
  if (length(lty) != ncol(object)) {
    lty <- rep(lty, length.out = ncol(object))
  }
  if (is.null(names(lwd))) {
    names(lwd) <- colnames(object)
  }
  return(list(col = col, lty = lty, lwd = lwd))
}

.covPlotLines <- function(x, y, col, lty, lwd) {
  if (sum(!is.na(y)) <= 1) {
    return(NULL)
  }
  lines(x, y, type = "s", col = col, lty = lty, lwd = lwd)
}

.plotATAC <- function(ATAC_SE, region, extend, addRegions, col, lty, lwd,
                      regionCol, highlightMain) {
  gr <- bsseq:::.bsGetGr(ATAC_SE, region, extend)
  ATAC_SE <- subsetByOverlaps(ATAC_SE, gr)

  ## Extract basic information
  sampleNames <- colnames(ATAC_SE)
  names(sampleNames) <- sampleNames
  positions <- start(ATAC_SE)
  coverage <- assay(ATAC_SE)

  ## get col, lwd, lty
  colEtc <- .ATACGetCol(object = ATAC_SE, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  plot(x = positions[1],
       y = 0.5,
       type = "n",
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1),
       xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "ATAC")
  axis(side = 2, at = c(0.2, 0.5, 0.8))
  .bsHighlightRegions(regions = addRegions,
                      gr = gr,
                      ylim = c(0, 1),
                      regionCol = regionCol,
                      highlightMain = highlightMain)

  sapply(1:ncol(ATAC_SE), function(sampIdx) {
    # NOTE: alpha is hardcoded
    .covPlotLines(positions, coverage[, sampIdx],
                  col = scales::alpha(colEtc$col[sampIdx], 1 / 6),
                  lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx])
  })
  unq_cols <- unique(colEtc$col)
  sapply(seq_along(unq_cols), function(i) {
    j <- colEtc$col %in% unq_cols[i]
    .covPlotLines(positions, rowMeans(coverage[, j]), col = unq_cols[i],
                  lty = colEtc$lty[j[1]], lwd = colEtc$lwd[j[1]])
  })
}

#-------------------------------------------------------------------------------
# plotRegion2(): Includes CG, CH, and ATAC-seq data on the one plot
#

# NOTE: Plotting of BSseqStat/BSseqTstat isn't supported
plotRegion2 <- function(list_of_CG_BSseq,
                        list_of_CH_BSseq,
                        ATAC_SE = NULL,
                        region = NULL,
                        extend = 0,
                        main = "",
                        list_of_CG_addRegions = NULL,
                        list_of_CH_addRegions = NULL,
                        ATAC_addRegions = NULL,
                        annoTrack = NULL,
                        cex.anno = 1,
                        geneTrack = NULL,
                        cex.gene = 1.5,
                        col = NULL,
                        lty = NULL,
                        lwd = NULL,
                        mainWithWidth = TRUE,
                        regionCol = alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        addNames = FALSE) {

  opar <- par(mar = c(0, 4.1, 0, 0),
              oma = c(5, 0, 4, 2),
              mfrow = c(1, 1))
  on.exit(par(opar))
  n_mCG_panels <- length(list_of_CG_BSseq)
  n_mCH_panels <- length(list_of_CH_BSseq)
  if (is.null(geneTrack)) {
    layout(mat =
             matrix(1:(1 + n_mCG_panels + n_mCH_panels - is.null(ATAC_SE)),
                    ncol = 1),
           heights =
             c(rep(2, 1 + n_mCG_panels + n_mCH_panels - is.null(ATAC_SE)), 1))
  } else {
    layout(mat =
             matrix(1:(2 + n_mCG_panels + n_mCH_panels - is.null(ATAC_SE)),
                    ncol = 1),
           heights =
             c(rep(2, 1 + n_mCG_panels + n_mCH_panels - is.null(ATAC_SE)), 1,
               0.3))
  }

  # Plot CG_BSseq objects
  lapply(names(list_of_CG_BSseq), function(n) {
    CG_BSseq <- list_of_CG_BSseq[[n]]
    .plotSmoothData(BSseq = CG_BSseq,
                    region = region,
                    extend = extend,
                    addRegions = list_of_CG_addRegions[[n]],
                    col = col,
                    lty = lty,
                    lwd = lwd,
                    regionCol = regionCol,
                    addTicks = addTicks,
                    addPoints = addPoints,
                    pointsMinCov = pointsMinCov,
                    highlightMain = highlightMain,
                    ylab = n)
  })
  gr <- bsseq:::.bsGetGr(list_of_CG_BSseq[[1]], region, extend)

  # Plot CH_BSseq objects
  lapply(names(list_of_CH_BSseq), function(n) {
    CH_BSseq <- list_of_CH_BSseq[[n]]
    .plotSmoothData(BSseq = CH_BSseq,
                    region = region,
                    extend = extend,
                    addRegions = list_of_CH_addRegions[[n]],
                    col = col,
                    lty = lty,
                    lwd = lwd,
                    regionCol = regionCol,
                    addTicks = addTicks,
                    addPoints = addPoints,
                    pointsMinCov = pointsMinCov,
                    highlightMain = highlightMain,
                    ylab = n)
  })

  # Plot ATAC object
  if (!is.null(ATAC_SE)) {
    .plotATAC(ATAC_SE = ATAC_SE,
              region = region,
              extend = extend,
              addRegions = ATAC_addRegions,
              col = col,
              lty = lty,
              lwd = lwd,
              regionCol = regionCol,
              highlightMain = highlightMain)
  }

  # Plot annotTrack
  if (!is.null(annoTrack)) {
    bsseq:::plotAnnoTrack(gr, annoTrack, cex.anno)
  }

  # Plot geneTrack
  if (!is.null(geneTrack)) {
    plotGeneTrack(gr, geneTrack, cex.gene, addNames = addNames)
  }

  if (main == "") {
    main <- bsseq:::.bsPlotTitle(gr = region,
                                 extend = extend,
                                 main = main,
                                 mainWithWidth = mainWithWidth)
  }
  mtext(side = 3, text = main, outer = TRUE, cex = 1)
}

# ------------------------------------------------------------------------------
# plotManyRegions2(): Includes CG, CH, and ATAC-seq data on the one plot
#

plotManyRegions2 <- function(list_of_CG_BSseq,
                             list_of_CH_BSseq,
                             ATAC_SE = NULL,
                             regions = NULL,
                             extend = 0,
                             main = "",
                             list_of_CG_addRegions = NULL,
                             list_of_CH_addRegions = NULL,
                             ATAC_addRegions = NULL,
                             annoTrack = NULL,
                             cex.anno = 1,
                             geneTrack = NULL,
                             cex.gene = 1.5,
                             col = NULL,
                             lty = NULL,
                             lwd = NULL,
                             mainWithWidth = TRUE,
                             regionCol = alpha("red", 0.1),
                             addTicks = TRUE,
                             addPoints = FALSE,
                             pointsMinCov = 5,
                             highlightMain = FALSE,
                             addNames = FALSE,
                             verbose = TRUE) {
  cat("[plotManyRegions2] preprocessing ...")
  if (!is.null(regions)) {
    if (is(regions, "data.frame")) {
      gr <- data.frame2GRanges(regions, keepColumns = FALSE)
    } else {
      gr <- regions
    }
    if (!is(gr, "GRanges")) {
      stop("'regions' needs to be either a 'data.frame' (with a single row) ",
           "or a 'GRanges' (with a single element)")
    }
  } else {
    gr <- granges(list_of_CG_BSseq[[1]])
  }
  gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
  list_of_CG_BSseq <- lapply(list_of_CG_BSseq, subsetByOverlaps, gr)
  list_of_CH_BSseq <- lapply(list_of_CH_BSseq, subsetByOverlaps, gr)

  if (length(start(list_of_CG_BSseq[[1]])) == 0 &&
      all(lengths(list_of_CH_BSseq) == 0)) {
    stop("No overlap between BSseq data and regions")
  }
  if (length(main) != length(gr)) {
    main <- rep(main, length = length(gr))
  }
  cat("done\n")
  for (ii in seq(along = gr)) {
    if (verbose) {
      cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n",
                  ii,
                  nrow(regions)))
    }
    plotRegion2(list_of_CG_BSseq = list_of_CG_BSseq,
                list_of_CH_BSseq = list_of_CH_BSseq,
                ATAC_SE = ATAC_SE,
                region = regions[ii,],
                extend = extend,
                col = col,
                lty = lty,
                lwd = lwd,
                main = main[ii],
                list_of_CG_addRegions = list_of_CG_addRegions,
                list_of_CH_addRegions = list_of_CH_addRegions,
                ATAC_addRegions = ATAC_addRegions,
                regionCol = regionCol,
                mainWithWidth = mainWithWidth,
                annoTrack = annoTrack,
                cex.anno = cex.anno,
                geneTrack = geneTrack,
                cex.gene = cex.gene,
                addTicks = addTicks,
                addPoints = addPoints,
                pointsMinCov = pointsMinCov,
                highlightMain = highlightMain,
                addNames = TRUE)
  }
}
