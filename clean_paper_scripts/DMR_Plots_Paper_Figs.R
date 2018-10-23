# Plots of DMRs to be used in paper
# Peter Hickey
# 2017-02-04

###-----------------------------------------------------------------------------
### Setup
###

#-------------------------------------------------------------------------------
# Load packages
#

library(GenomicRanges)
library(bsseq)
library(scales)
library(dplyr)

#-------------------------------------------------------------------------------
# Load objects
#

load("../Objects/objects_for_plotting_DMR_and_block_examples.rda")
load("../integrating-dmrs-daps-and-degs/objects/assays-and-features.rda")

#-------------------------------------------------------------------------------
# Adjust the NeuN_color scheme so as not to clash with Tissue_color scheme
#

# TODO: Check whether these 2 lines are necessary (seem to have already saved
#       an object with these colour changes in
#       ../Objects/objects_for_plotting_DMR_and_block_examples.rda)
# BSseq_dmrs$NeuN_color <- ifelse(BSseq_dmrs$NeuN_color == "deepskyblue",
#                                 "darkgreen", "purple")
# BSseq_blocks$NeuN_color <- ifelse(BSseq_blocks$NeuN_color == "deepskyblue",
#                                   "darkgreen", "purple")
# atac_cov_dmrs <- lapply(atac_cov_dmrs, function(se) {
#   se$NEUN_COLOR <- ifelse(se$NEUN_COLOR == "deepskyblue",
#                           "darkgreen", "purple")
#   se
# })
# atac_cov_blocks <- lapply(atac_cov_blocks, function(se) {
#   se$NEUN_COLOR <- ifelse(se$NEUN_COLOR == "deepskyblue",
#                           "darkgreen", "purple")
#   se
# })

#-------------------------------------------------------------------------------
# Make POS-only objects
#

BSseq_pos_dmrs <- BSseq_dmrs[, BSseq_dmrs$NeuN == "pos"]
BSseq_pos_blocks <- BSseq_blocks[, BSseq_blocks$NeuN == "pos"]
atac_cov_pos_dmrs <- lapply(atac_cov_dmrs, function(se) {
  se[, se$NEUN == "pos"]
})
atac_cov_pos_blocks <- lapply(atac_cov_blocks, function(se) {
  se[, se$NEUN == "pos"]
})

#-------------------------------------------------------------------------------
# Make POS-only, no NAcc objects
#

BSseq_pos_no_NA_dmrs <- BSseq_pos_dmrs[, BSseq_pos_dmrs$Tissue != "NA"]
BSseq_pos_no_NA_blocks <- BSseq_pos_blocks[, BSseq_pos_blocks$Tissue != "NA"]
atac_cov_pos_no_NA_dmrs <- lapply(atac_cov_pos_dmrs, function(se) {
  se[, se$TISSUE != "NA"]
})
atac_cov_pos_no_NA_blocks <- lapply(atac_cov_pos_blocks, function(se) {
  se[, se$TISSUE != "NA"]
})

###-----------------------------------------------------------------------------
### Functions
###

plotAnnoTrack <- function(gr, annoTrack, cex) {
  if (!all(sapply(annoTrack, function(xx) is(xx, "GRanges"))))
    stop("all elements in 'annoTrack' needs to be 'GRanges'")
  plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
       ylim = c(0.5, length(annoTrack) + 0.5),
       xlim = c(start(gr),end(gr)), xlab = "", ylab = "")
  lapply(seq(along = annoTrack), function(ii) {
    jj <- length(annoTrack) + 1 - ii
    ir <- subsetByOverlaps(annoTrack[[ii]], gr)
    if (length(ir) > 0)
      rect(start(ir) - 0.5, jj - 0.1, end(ir), jj + 0.1,
           col = "grey60", border = NA)
    mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1,
          line = 1, cex = cex, adj = 0.3)
  })
}

# NOTE: `flank` is to try to ensure that 'nearby' exons are identified
#       (to avoid corner case region being plotted is entirely within an intron
#       in which case no exons are found by the call to findOverlaps() meaning
#       that the gene model is never plotted)
plotGeneTrack <- function(gr, geneTrack, cex, flank = 20000) {
  geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
  ov <- findOverlaps(geneTrack_gr, resize(gr, width(gr) + 2 * flank,
                                          fix = "center"))
  genes <- geneTrack[queryHits(ov), ]
  plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
       ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)), xlab = "",
       ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
  # NOTE: `gene_names` and `this_gene` are used to ensure gene names are
  #        only printed once
  gene_names <- c()
  if (nrow(genes) > 0) {
    for (g in 1:nrow(genes)) {
      this_gene <- genes$gene_name[g]
      geneind2 = which(geneTrack$gene_name == genes$gene_name[g])
      geneind2 = geneind2[which(geneTrack$isoforms[geneind2] ==
                                  1)]
      direction = unique(geneTrack$strand[geneind2])
      ES = geneTrack$start[geneind2]
      EE = geneTrack$end[geneind2]
      Exons = cbind(ES, EE)
      if (direction == "+") {
        lines(x = c(min(ES), max(EE)), y = c(0.65, 0.65))
        apply(Exons, 1, function(x) {
          polygon(c(x[1], x[2], x[2], x[1]), c(0.35, 0.35, 0.95, 0.95),
                  col = "darkgrey")
        })
        if (!this_gene %in% gene_names) {
          text((max(start(gr), min(ES)) + min(end(gr),
                                              max(EE)))/2, 1.25,
               genes$gene_name[g], cex = cex)
          gene_names <- c(gene_names, this_gene)
        }
      }
      else {
        lines(x = c(min(ES), max(EE)), y = c(-0.65, -0.65))
        apply(Exons, 1, function(x) polygon(c(x[1], x[2],
                                              x[2], x[1]),
                                            c(-0.35, -0.35, -0.95, -0.95),
                                            col = "darkgrey"))
        if (!this_gene %in% gene_names) {
          text((max(start(gr), min(ES)) + min(end(gr),
                                              max(EE)))/2, -1.25,
               genes$gene_name[g], cex = cex)
          gene_names <- c(gene_names, this_gene)
        }
      }
    }
  }
}

.covGetCol <- function(object, col, lty, lwd) {
  if (is.null(col)) {
    if ("col" %in% names(colData(object)))
      col <- colData(object)[["col"]]
    else col <- rep("black", nrow(colData(object)))
  }
  if (length(col) != ncol(object))
    col <- rep(col, length.out = ncol(object))
  if (is.null(names(col)))
    names(col) <- colnames(object)
  if (is.null(lty)) {
    if ("lty" %in% names(colData(object)))
      lty <- colData(object)[["lty"]]
    else lty <- rep(1, ncol(object))
  }
  if (length(lty) != ncol(object))
    lty <- rep(lty, length.out = ncol(object))
  if (is.null(names(lty)))
    names(lty) <- colnames(object)
  if (is.null(lwd)) {
    if ("lwd" %in% names(colData(object)))
      lwd <- colData(object)[["lwd"]]
    else lwd <- rep(1, nrow(colData(object)))
  }
  if (length(lty) != ncol(object))
    lty <- rep(lty, length.out = ncol(object))
  if (is.null(names(lwd)))
    names(lwd) <- colnames(object)
  return(list(col = col, lty = lty, lwd = lwd))
}

.covPlotLines <- function(x, y, col, lty, lwd) {
  if (sum(!is.na(y)) <= 1) {
    return(NULL)
  }
  lines(x, y, type = "s", col = col, lty = lty, lwd = lwd)
}

.plotCovData <- function(cov, region, extend, addRegions, col, lty, lwd,
                         regionCol, highlightMain, k = 1) {
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  cov <- subsetByOverlaps(cov, gr)

  ## Extract basic information
  sampleNames <- colnames(cov)
  names(sampleNames) <- sampleNames
  positions <- start(cov)
  coverage <- assay(cov)

  ## Smooth coverage
  coverage <- do.call(cbind, lapply(seq_len(ncol(coverage)), function(j) {
    runmed(coverage[, j], k)
  }))

  ## get col, lwd, lty
  colEtc <- .covGetCol(object = cov, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  # NOTE: y-lim is hard-coded
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 2), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "ATAC signal")
  axis(side = 2, at = c(0.4, 1, 1.6))

  # NOTE: y-lim is hard-coded
  bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0, 2),
                              regionCol = regionCol,
                              highlightMain = highlightMain)

  sapply(1:ncol(cov), function(sampIdx) {
    # NOTE: alpha is hardcoded
    bsseq:::.bsPlotLines(x = positions,
                         y = coverage[, sampIdx],
                         col = scales::alpha(colEtc$col[sampIdx], 1 / 6),
                         lty = colEtc$lty[sampIdx],
                         lwd = colEtc$lwd[sampIdx],
                         plotRange <- range(positions))
    # .covPlotLines(positions, coverage[, sampIdx],
    #               col = scales::alpha(colEtc$col[sampIdx], 1 / 6),
    #              )
  })
  unq_cols <- unique(colEtc$col)
  sapply(seq_along(unq_cols), function(i) {
    j <- colEtc$col %in% unq_cols[i]
    # .covPlotLines(positions, rowMeans(coverage[, j]), col = unq_cols[i],
    #               lty = colEtc$lty[j[1]], lwd = colEtc$lwd[j[1]])
    bsseq:::.bsPlotLines(x = positions,
                         y = rowMeans(coverage[, j], na.rm = TRUE),
                         col = unq_cols[i],
                         lty =  colEtc$lty[j[1]],
                         lwd = colEtc$lwd[j[1]],
                         plotRange <- range(positions))
  })
}

plotRegionWithATAC <- function(cov,
                               cov.region = NULL,
                               cov.addRegions = NULL,
                               cov.col = NULL,
                               cov.lty = NULL,
                               cov.lwd = NULL,
                               cov.regionCol = alpha("red", 0.1),
                               k = 1,
                               BSseq,
                               region = NULL,
                               extend = 0,
                               main = "",
                               addRegions = NULL,
                               annoTrack = NULL,
                               cex.anno = 1,
                               geneTrack = NULL,
                               cex.gene = 1.5,
                               col = NULL,
                               lty = NULL,
                               lwd = NULL,
                               BSseqStat = NULL,
                               stat = "tstat.corrected",
                               stat.col = "black",
                               stat.lwd = 1,
                               stat.lty = 1,
                               stat.ylim = c(-8, 8),
                               mainWithWidth = TRUE,
                               regionCol = alpha("red", 0.1),
                               addTicks = TRUE,
                               addPoints = FALSE,
                               pointsMinCov = 5,
                               highlightMain = FALSE,
                               v1 = NULL,
                               v2 = NULL) {

  # NOTE: This chunk normally would go in plotManyRegionsWithATAC()
  if (!is.null(region)) {
    if (is(region, "data.frame")) {
      gr <- data.frame2GRanges(region, keepColumns = FALSE)
    } else {
      gr <- region
    }
    if (!is(gr, "GRanges")) {
      stop("'region' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
    }
  } else {
    gr <- granges(BSseq)
  }
  gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
  BSseq <- subsetByOverlaps(BSseq, gr)
  # Subset list_of_se by gr to create a single SummarizedExperiment object
  # where the 'cov' assay is now an ordinary matrix
  # NOTE: repete::os(cov) will be large if sum(width(regions)) is large
  cov <- do.call(rbind, lapply(seq_along(gr), function(i) {
    se <- subsetByOverlaps(cov[[as.character(seqnames(gr[i]))]], gr[i])
    assay(se) <- do.call(cbind, lapply(assay(se), as.matrix))
    se
  }))
  if (!is.null(BSseqStat)) {
    BSseqStat <- subsetByOverlaps(BSseqStat, gr)
  }

  opar <- par(mar = c(0, 4.1, 0, 0), oma = c(5, 0, 4, 2), mfrow = c(1, 1))
  on.exit(par(opar))
  # if (is.null(BSseqStat) && is.null(geneTrack)) {
  #   print("A")
  #   layout(matrix(1:3, ncol = 1), heights = c(2, 2, 1))
  # } else if (is.null(geneTrack)) {
  #   print("B")
  #   layout(matrix(1:4, ncol = 1), heights = c(2, 2, 2, 1))
  # } else {
  #   print("C")
  #   layout(matrix(1:5, ncol = 1), heights = c(2, 2, 2, 1, 0.3))
  # }
  layout(matrix(1:4, ncol = 1), heights = c(3, 3, 1, 1))

  # Plot BSseq object
  bsseq:::.plotSmoothData(BSseq = BSseq,
                          region = region,
                          extend = extend,
                          addRegions = addRegions,
                          col = col,
                          lty = lty,
                          lwd = lwd,
                          regionCol = regionCol,
                          addTicks = addTicks,
                          addPoints = addPoints,
                          pointsMinCov = pointsMinCov,
                          highlightMain = highlightMain)
  if (!is.null(v1) & !is.null(v2)) {
    abline(v = v1)
    abline(v = v2)
  }
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)

  # Plot cov object
  .plotCovData(cov = cov,
               region = region,
               extend = extend,
               addRegions = cov.addRegions,
               col = cov.col,
               lty = cov.lty,
               lwd = cov.lwd,
               regionCol = cov.regionCol,
               highlightMain = highlightMain,
               k = k)
  if (!is.null(v1) & !is.null(v2)) {
    abline(v = v1)
    abline(v = v2)
  }

  if (!is.null(BSseqStat)) {
    BSseqStat <- subsetByOverlaps(BSseqStat, gr)
    if (is(BSseqStat, "BSseqTstat")) {
      stat.values <- getStats(BSseqStat)[, stat]
      stat.type <- "tstat"
    }
    if (is(BSseqStat, "BSseqStat")) {
      stat.type <- getStats(BSseqStat, what = "stat.type")
      if (stat.type == "tstat") {
        stat.values <- getStats(BSseqStat, what = "stat")
      }
      if (stat.type == "fstat") {
        stat.values <- sqrt(getStats(BSseqStat, what = "stat"))
      }
    }
    # NOTE: Need to realise in memory the data to be plotted
    if (is(stat.values, "DelayedArray")) {
      stat.values <- as.array(stat.values)
    }
    plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "",
         ylab = stat.type)
    axis(side = 2, at = c(-5,0,5))
    abline(h = 0, col = "grey60")
    .bsPlotLines(start(BSseqStat), stat.values, lty = stat.lty, col = stat.col,
                 lwd = stat.lwd, plotRange = c(start(gr), end(gr)))
  }

  if (!is.null(annoTrack)) {
    plotAnnoTrack(gr, annoTrack, cex.anno)
  }

  if (!is.null(geneTrack)) {
    plotGeneTrack(gr, geneTrack, cex.gene)
  }

  if (!is.null(main)) {
    main <- bsseq:::.bsPlotTitle(gr = region, extend = extend, main = main,
                                 mainWithWidth = mainWithWidth)
    mtext(side = 3, text = main, outer = TRUE, cex = 1)
  }
}

###-----------------------------------------------------------------------------
### Plot design
###

# - Highlight "significant" ATAC peaks as we do DMRs
# - No rug for blocks (and perhaps neither for DMRs)
# - Customise 'extend' based on Lindsay's notes
# - Brain enhancers
# - RefSeq gene models (check with UCSC genome browser that these look correct)

###-----------------------------------------------------------------------------
### Figure 2
###

pdf("../Graphs/Figure-2-for-editing.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)

#-------------------------------------------------------------------------------
# (A) Venn Diagram.
#

# NOTE: Nothing to do.

#-------------------------------------------------------------------------------
# (B) LRRC4 DMR
#

# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["LRRC4"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(start(dmrs["LRRC4"]) - 2000, start(dmrs["LRRC4"]) - 1000),
      y = c(0.85, 0.85))
text("1 kb", x = start(dmrs["LRRC4"]) - 1500, y = 0.9)
legend(x = end(dmrs["LRRC4"]) - 50,
       y = 1.25,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 1)

#-------------------------------------------------------------------------------
# (C) DLGAP2 DMR (originally GRIK4)
#

# NOTE: Lindsay noted in Prepare_objecs_for_DMR_plots.R
#       "COME BACK TO THIS ONE...not in original 208 DMR list instead of GRIK4
#       can use this one ...... DLGAP2
# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_pos_no_NA_dmrs,
                   cov.region = NULL,
                   # NOTE: No DA ATAC-seq analysis for this comparison
                   cov.addRegions = NULL,
                   cov.col = atac_cov_pos_no_NA_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_pos_no_NA_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_no_NA_dmrs,
                   region = dmrs[["DLGAP2"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_non_NA_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_no_NA_dmrs$Tissue_color,
                   lty = BSseq_pos_no_NA_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(start(dmrs["DLGAP2"]) - 2000, start(dmrs["DLGAP2"]) - 1000),
      y = c(0.6, 0.6))
text("1 kb", x = start(dmrs["DLGAP2"]) - 1500, y = 0.65)
legend(x = end(dmrs["DLGAP2"]) + 700,
       y = 0.35,
       legend = unique(BSseq_pos_no_NA_dmrs$Tissue),
       col = unique(BSseq_pos_no_NA_dmrs$Tissue_color),
       lwd = 2,
       bty = "n")

#-------------------------------------------------------------------------------
# (D) GABRB2 BLOCK
#

# NOTE: GABRB2 is on reverse strand
gabrb2_cutout <- GRanges(
  RefSeq_exons3_blocks[RefSeq_exons3_blocks$gene_name == "GABRB2" &
                         RefSeq_exons3_blocks$exon_number == 4 &
                         RefSeq_exons3_blocks$isoforms == 1, ])


# NOTE: Extend 100kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_blocks,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_blocks[
                       abs(atac_NA_posvsBA9_pos_blocks$logFC) > 1],
                   cov.col = atac_cov_blocks[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_blocks[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_blocks,
                   region = blocks[["GABRB2_block"]],
                   extend = 40000,
                   main = "",
                   addRegions = pos_blocks,
                   annoTrack = list("BrainEnh" = Brain_enh_blocks),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_blocks,
                   cex.gene = 1.5,
                   col = BSseq_pos_blocks$Tissue_color,
                   lty = BSseq_pos_blocks$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE,
                   v1 = start(gabrb2_cutout) - 4000,
                   v2 = end(gabrb2_cutout) + 4000)
lines(x = c(start(blocks["GABRB2_block"]) + 10000,
            start(blocks["GABRB2_block"]) + 60000),
      y = c(0.8, 0.8))
text("50 kb", x = start(blocks["GABRB2_block"]) + 35000, y = 0.85)
legend(x = end(blocks["GABRB2_block"]) - 80000,
       y = 0.7,
       legend = sub("NA", "NAcc", unique(BSseq_pos_blocks$Tissue)),
       col = unique(BSseq_pos_blocks$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

# Cutout
plotRegionWithATAC(cov = atac_cov_blocks,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_blocks[
                     abs(atac_NA_posvsBA9_pos_blocks$logFC) > 1],
                   cov.col = atac_cov_blocks[[1]]$TISSUE_COLOR,
                   cov.lty = 1,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_blocks,
                   region = gabrb2_cutout,
                   extend = 4000,
                   main = "",
                   addRegions = pos_blocks,
                   annoTrack = list("BrainEnh" = Brain_enh_blocks),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_blocks,
                   cex.gene = 1.5,
                   col = BSseq_pos_blocks$Tissue_color,
                   lty = 1,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(end(gabrb2_cutout) + 2000,
            end(gabrb2_cutout) + 3000),
      y = c(0.1, 0.1))
text("1 kb", x = end(gabrb2_cutout) + 2500, y = 0.15)
legend(x = start(gabrb2_cutout) - 4000,
       y = 0.35,
       legend = sub("NA", "NAcc", unique(BSseq_pos_blocks$Tissue)),
       col = unique(BSseq_pos_blocks$Tissue_color),
       lwd = 2,
       bty = "n")

#-------------------------------------------------------------------------------
# (E) RBFOX3 DMR
#

# NOTE: Plot POS vs NEG
# NOTE: Extend 5kb

# TODO: Gene model not being plotted

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_ave_pos_vs_ave_neg_dmrs[
                     abs(atac_ave_pos_vs_ave_neg_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$NEUN_COLOR,
                   cov.lty = 1,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_dmrs,
                   region = dmrs[["RBFOX3"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_vs_neg_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_dmrs$NeuN_color,
                   lty = 1,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(start(dmrs["RBFOX3"]) - 3000, start(dmrs["RBFOX3"]) - 2000),
      y = c(0.1, 0.1))
text("1 kb", x = start(dmrs["RBFOX3"]) - 2500, y = 0.15)
legend(x = end(dmrs["RBFOX3"]) + 300,
       y = .3,
       legend = unique(BSseq_dmrs$NeuN),
       col = unique(BSseq_dmrs$NeuN_color),
       lwd = 2,
       bty = "n")

#-------------------------------------------------------------------------------
# (F) QKI BLOCK
#

# NOTE: Plot POS vs NEG
# NOTE: Extend 100kb and include cutout

qki_first_exon <- GRanges(
  RefSeq_exons3_blocks[RefSeq_exons3_blocks$gene_name == "QKI" &
                         RefSeq_exons3_blocks$exon_number == 0 &
                         RefSeq_exons3_blocks$isoforms == 1, ])
qki_cutout <- shift(GRanges(
  RefSeq_exons3_blocks[RefSeq_exons3_blocks$gene_name == "QKI" &
                         RefSeq_exons3_blocks$exon_number == 7 &
                         RefSeq_exons3_blocks$isoforms == 1, ]),
  shift = 3000)

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_blocks,
                   cov.region = NULL,
                   cov.addRegions = atac_ave_pos_vs_ave_neg_blocks[
                     abs(atac_ave_pos_vs_ave_neg_blocks$logFC) > 1],
                   cov.col = atac_cov_blocks[[1]]$NEUN_COLOR,
                   cov.lty = 1,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_blocks,
                   region = blocks[["QKI_block"]],
                   extend = 40000,
                   main = "",
                   addRegions = pos_vs_neg_blocks,
                   annoTrack = list("BrainEnh" = Brain_enh_blocks),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_blocks,
                   cex.gene = 1.5,
                   col = BSseq_blocks$NeuN_color,
                   lty = 1,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE,
                   v1 = start(qki_cutout) - 4000,
                   v2 = end(qki_cutout) + 4000)
lines(x = c(start(blocks["QKI_block"]) + 20000,
            start(blocks["QKI_block"]) + 70000),
      y = c(0.25, 0.25))
text("50 kb", x = start(blocks["QKI_block"]) + 45000, y = 0.3)
legend(x = end(blocks["QKI_block"]) - 110000,
       y = 0.2,
       legend = unique(BSseq_blocks$NeuN),
       col = unique(BSseq_blocks$NeuN_color),
       lwd = 2,
       bty = "n")

# Cutout
plotRegionWithATAC(cov = atac_cov_blocks,
                   cov.region = NULL,
                   cov.addRegions = atac_ave_pos_vs_ave_neg_blocks[
                     abs(atac_ave_pos_vs_ave_neg_blocks$logFC) > 1],
                   cov.col = atac_cov_blocks[[1]]$NEUN_COLOR,
                   cov.lty = 1,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_blocks,
                   region = qki_cutout,
                   extend = 4000,
                   main = "",
                   addRegions = pos_vs_neg_blocks,
                   annoTrack = list("BrainEnh" = Brain_enh_blocks),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_blocks,
                   cex.gene = 1.5,
                   col = BSseq_blocks$NeuN_color,
                   lty = 1,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(end(qki_cutout) + 2000,
            end(qki_cutout) + 3000),
      y = c(0.1, 0.1))
text("1 kb", x = end(qki_cutout) + 2500, y = 0.15)
legend(x = start(qki_cutout) - 4000,
       y = 0.2,
       legend = unique(BSseq_blocks$NeuN),
       col = unique(BSseq_blocks$NeuN_color),
       lwd = 2,
       bty = "n")

dev.off()

###-----------------------------------------------------------------------------
### Figure 3
###

pdf("../Graphs/Figure-3-for-editing.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)

#-------------------------------------------------------------------------------
# (A) MEF2C BLOCK
#

# NOTE: Extend 100kb
# NOTE: Decided to remove this from Fig 3
#       (https://jhu-genomics.slack.com/archives/hansen_gtex/p1476372569000129)
#       but will still create figure here
# NOTE: LINC00461 and MIR9−2 overlap (and appear to share exons) which makes
#       it hard to plot (http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr5%3A87879404-88287364&hgsid=547790305_9iS9OaTyQXiaP7VOggaf2M55aPOA)

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_blocks,
                   cov.region = NULL,
                   cov.addRegions =  atac_NA_posvsBA9_pos_blocks[
                       abs(atac_NA_posvsBA9_pos_blocks$logFC) > 1],
                   cov.col = atac_cov_blocks[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_blocks[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_blocks,
                   region = blocks[["MEF2C_block"]],
                   extend = 40000,
                   main = "",
                   addRegions = pos_blocks,
                   annoTrack = list("BrainEnh" = Brain_enh_blocks),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_blocks,
                   cex.gene = 1.5,
                   col = BSseq_pos_blocks$Tissue_color,
                   lty = BSseq_pos_blocks$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE,
                   v1 = start(dmrs[["MEF2C_small_DMR"]]) - 4000,
                   v2 = end(dmrs[["MEF2C_small_DMR"]]) + 4000)
lines(x = c(start(blocks[["MEF2C_block"]]) + 100000,
            start(blocks[["MEF2C_block"]]) + 150000),
      y = c(0.45, 0.45))
text("50 kb", x = start(blocks[["MEF2C_block"]]) + 125000, y = 0.5)
legend(x = end(blocks[["MEF2C_block"]]) - 80000,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_blocks$Tissue)),
       col = unique(BSseq_pos_blocks$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

#-------------------------------------------------------------------------------
# (B) RGS9 DMR
#

# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["RGS9"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(end(dmrs[["RGS9"]]) + 2500, end(dmrs[["RGS9"]]) + 3500),
      y = c(0.5, 0.5))
text("1 kb", x = end(dmrs[["RGS9"]]) + 3000, y = 0.55)
legend(x = start(dmrs[["RGS9"]]) - 4200,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

#-------------------------------------------------------------------------------
# (C) SEMA7A DMR
#

# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["SEMA7A"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(start(dmrs[["SEMA7A"]]) - 3000, start(dmrs[["SEMA7A"]]) - 2000),
      y = c(0.5, 0.5))
text("1 kb", x = start(dmrs[["SEMA7A"]]) - 2500, y = 0.55)
legend(x = end(dmrs[["SEMA7A"]]) - 4000,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

#-------------------------------------------------------------------------------
# (D) SATB2 DMR
#

# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["SATB2"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(start(dmrs[["SATB2"]]) - 3000, start(dmrs[["SATB2"]]) - 2000),
      y = c(0.5, 0.5))
text("1 kb", x = start(dmrs[["SATB2"]]) - 2500, y = 0.55)
legend(x = start(dmrs[["SATB2"]]) - 4000,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

#-------------------------------------------------------------------------------
# (E) MEF2C DMR
#

# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["MEF2C_small_DMR"]],
                   extend = 4000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(end(dmrs[["MEF2C_small_DMR"]]) + 1500,
            end(dmrs[["MEF2C_small_DMR"]]) + 2500),
      y = c(0.2, 0.2))
text("1 kb", x = end(dmrs[["MEF2C_small_DMR"]]) + 2000, y = 0.25)
legend(x = start(dmrs[["MEF2C_small_DMR"]]) - 4000,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

dev.off()

#-------------------------------------------------------------------------------
# (bonus) SIX3 DMR
#

pdf("../Graphs/SIX3-for-editing.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
# NOTE: Extend 5kb

# With ATAC-seq data
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1,
                   BSseq = BSseq_pos_dmrs,
                   region = dmrs[["SIX3"]],
                   extend = 25000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
lines(x = c(end(dmrs[["SIX3"]]) + 1500,
            end(dmrs[["SIX3"]]) + 2500),
      y = c(0.2, 0.2))
text("1 kb", x = end(dmrs["SIX3"]) + 2000, y = 0.25)
legend(x = start(dmrs[["SIX3"]]) - 4000,
       y = 0.2,
       legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
       col = unique(BSseq_pos_dmrs$Tissue_color),
       lwd = 2,
       bty = "n",
       ncol = 2)

dev.off()

#-------------------------------------------------------------------------------
# Enhancer-rich genes
#

pdf("../Graphs/enhancer_rich_genes.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
for (i in grep("ENS", names(dmrs))) {
  # With ATAC-seq data
  plotRegionWithATAC(cov = atac_cov_dmrs,
                     cov.region = NULL,
                     cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                     cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                     cov.lty = atac_cov_dmrs[[1]]$LTY,
                     cov.lwd = NULL,
                     cov.regionCol = alpha("red", 0.1),
                     k = 1001,
                     BSseq = BSseq_pos_dmrs,
                     region = dmrs[[i]],
                     extend = 40000,
                     main = "",
                     addRegions = pos_dmrs,
                     annoTrack = list("Linked" =
                                        unique(unlist(fantom5_enhancers_by_gene)),
                                      "H3K27ac" = H3K27ac_brain),
                     cex.anno = 1,
                     geneTrack = RefSeq_exons3_dmrs,
                     cex.gene = 1.5,
                     col = BSseq_pos_dmrs$Tissue_color,
                     lty = BSseq_pos_dmrs$lty,
                     lwd = NULL,
                     BSseqStat = NULL,
                     stat = "tstat.corrected",
                     stat.col = "black",
                     stat.lwd = 1,
                     stat.lty = 1,
                     stat.ylim = c(-8, 8),
                     mainWithWidth = TRUE,
                     regionCol = alpha("red", 0.1),
                     addTicks = FALSE,
                     addPoints = FALSE,
                     pointsMinCov = 5,
                     highlightMain = FALSE)
  # lines(x = c(start(dmrs[i]) + 10000,
  #             start(dmrs[i]) + 60000),
  #       y = c(0.8, 0.8))
  # text("50 kb", x = start(dmrs[i]) + 35000, y = 0.85)
  # legend(x = end(dmrs[i]) - 80000,
  #        y = 0.7,
  #        legend = sub("NA", "NAcc", unique(BSseq_pos_dmrs$Tissue)),
  #        col = unique(BSseq_pos_dmrs$Tissue_color),
  #        lwd = 2,
  #        bty = "n",
  #        ncol = 2)
}
dev.off()

#-------------------------------------------------------------------------------
# DLX6 (ENSG00000006377.9) locus
#

DLX6_locus <- reduce(unstrand(
  c(granges(unflattened_features$genes["ENSG00000006377.9"]),
    granges(fantom5_enhancers_by_gene[["ENSG00000006377.9"]]))))
DLX6_locus <- GRanges(seqnames(DLX6_locus[1]),
                      IRanges(min(start(DLX6_locus)),
                              max(end(DLX6_locus))))
DLX6_enhancers <- fantom5_enhancers_by_gene[["ENSG00000006377.9"]]

pdf("../Graphs/DLX6.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = DLX6_locus,
                   extend = 10000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("Permissive" =
                                      permissive_enhancers,
                                    "DLX6-linked" = DLX6_enhancers,
                                    "H3K27ac" = H3K27ac_brain),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
dev.off()

#-------------------------------------------------------------------------------
# MEF2C (ENSG00000081189.9) locus
#

MEF2C_locus <- reduce(unstrand(
  c(granges(unflattened_features$genes["ENSG00000081189.9"]),
    granges(fantom5_enhancers_by_gene[["ENSG00000081189.9"]]))))
MEF2C_locus <- GRanges(seqnames(MEF2C_locus[1]),
                      IRanges(min(start(MEF2C_locus)),
                              max(end(MEF2C_locus))))
MEF2C_enhancers <- fantom5_enhancers_by_gene[["ENSG00000081189.9"]]

pdf("../Graphs/MEF2C.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = MEF2C_locus,
                   extend = 10000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("Permissive" =
                                      permissive_enhancers,
                                    "MEF2C-linked" = MEF2C_enhancers,
                                    "H3K27ac" = H3K27ac_brain),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
dev.off()

#-------------------------------------------------------------------------------
# ZBTB16 (ENSG00000109906.9) locus
#

ZBTB16_locus <- reduce(unstrand(
  c(granges(unflattened_features$genes["ENSG00000109906.9"]),
    granges(fantom5_enhancers_by_gene[["ENSG00000109906.9"]]))))
ZBTB16_locus <- GRanges(seqnames(ZBTB16_locus[1]),
                       IRanges(min(start(ZBTB16_locus)),
                               max(end(ZBTB16_locus))))
ZBTB16_enhancers <- fantom5_enhancers_by_gene[["ENSG00000109906.9"]]

pdf("../Graphs/ZBTB16.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = ZBTB16_locus,
                   extend = 10000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("Permissive" =
                                      permissive_enhancers,
                                    "ZBTB16-linked" = ZBTB16_enhancers,
                                    "H3K27ac" = H3K27ac_brain),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
dev.off()

#-------------------------------------------------------------------------------
# ETS1 (ENSG00000134954.10) locus
#

ETS1_locus <- reduce(unstrand(
  c(granges(unflattened_features$genes["ENSG00000134954.10"]),
    granges(fantom5_enhancers_by_gene[["ENSG00000134954.10"]]))))
ETS1_locus <- GRanges(seqnames(ETS1_locus[1]),
                        IRanges(min(start(ETS1_locus)),
                                max(end(ETS1_locus))))
ETS1_enhancers <- fantom5_enhancers_by_gene[["ENSG00000134954.10"]]

pdf("../Graphs/ETS1.pdf",
    width = 7,
    height = 7,
    useDingbats = FALSE)
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = ETS1_locus,
                   extend = 10000,
                   main = "",
                   addRegions = pos_dmrs,
                   annoTrack = list("Permissive" =
                                      permissive_enhancers,
                                    "ETS1-linked" = ETS1_enhancers,
                                    "H3K27ac" = H3K27ac_brain),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)
dev.off()

#-------------------------------------------------------------------------------
# KLF5 (ENSG00000102554.9) locus
#

KLF5_locus <- reduce(unstrand(
  c(granges(unflattened_features$genes["ENSG00000102554.9"]),
    granges(fantom5_enhancers_by_gene[["ENSG00000102554.9"]]))))
KLF5_locus <- GRanges(seqnames(KLF5_locus[1]),
                      IRanges(min(start(KLF5_locus)),
                              max(end(KLF5_locus))))
KLF5_gene <- unflattened_features$genes["ENSG00000102554.9"]
KLF5_enhancers <- fantom5_enhancers_by_gene[["ENSG00000102554.9"]]

pdf("../Graphs/KLF5-for-editing-1.pdf",
    width = 5,
    height = 7,
    useDingbats = FALSE)
plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   # cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.col = scales::alpha("white", 0),
                   # cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lty = 0,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.3),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = KLF5_locus,
                   extend = 10000,
                   main = NULL,
                   addRegions = pos_dmrs,
                   annoTrack = list("cutout" = c(resize(KLF5_gene,
                                                        width(KLF5_gene) +
                                                          2 * 5000,
                                                        fix = "center"),
                                                 resize(KLF5_enhancers,
                                                        width(KLF5_enhancers) +
                                                          2 * 5000,
                                                        fix = "center"),
                                                 ignore.mcols = TRUE),
                                    "KLF5-linked" = KLF5_enhancers,
                                    "BrainEnh" = Brain_enh_dmrs,
                                    "DMRs" = pos_dmrs,
                                    "DAPs" = atac_NA_posvsBA9_pos_dmrs[
                                      abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1]),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   # col = BSseq_pos_dmrs$Tissue_color,
                   col = scales::alpha("white", 0),
                   # lty = BSseq_pos_dmrs$lty,
                   lty = 0,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.3),
                   addTicks = FALSE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)

op <- par(no.readonly = TRUE)
opar <- par(mar = c(0, 4.1, 0, 0), oma = c(5, 0, 4, 2), mfrow = c(1, 1))
layout(matrix(1:4, ncol = 1), heights = c(3, 3, 1, 1))
plot(x = start(KLF5_locus),
     y = 1,
     pch = NA_integer_,
     xlim = c(start(KLF5_locus) - 5000, end(KLF5_locus) + 5000),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
lines(x = c(start(KLF5_locus) + 10000,
            start(KLF5_locus) + 60000),
      y = c(0.8, 0.8))
text("50 kb", x = start(KLF5_locus) + 35000, y = 0.85)
par(op)

plotRegionWithATAC(cov = atac_cov_dmrs,
                   cov.region = NULL,
                   cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                     abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                   cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                   cov.lty = atac_cov_dmrs[[1]]$LTY,
                   cov.lwd = NULL,
                   cov.regionCol = alpha("red", 0.1),
                   k = 1001,
                   BSseq = BSseq_pos_dmrs,
                   region = KLF5_gene,
                   extend = 5000,
                   main = NULL,
                   addRegions = pos_dmrs,
                   annoTrack = list("KLF5-linked" = KLF5_enhancers,
                                    "BrainEnh" = Brain_enh_dmrs),
                   cex.anno = 1,
                   geneTrack = RefSeq_exons3_dmrs,
                   cex.gene = 1.5,
                   col = BSseq_pos_dmrs$Tissue_color,
                   lty = BSseq_pos_dmrs$lty,
                   lwd = NULL,
                   BSseqStat = NULL,
                   stat = "tstat.corrected",
                   stat.col = "black",
                   stat.lwd = 1,
                   stat.lty = 1,
                   stat.ylim = c(-8, 8),
                   mainWithWidth = TRUE,
                   regionCol = alpha("red", 0.1),
                   addTicks = TRUE,
                   addPoints = FALSE,
                   pointsMinCov = 5,
                   highlightMain = FALSE)

op <- par(no.readonly = TRUE)
opar <- par(mar = c(0, 4.1, 0, 0), oma = c(5, 0, 4, 2), mfrow = c(1, 1))
layout(matrix(1:4, ncol = 1), heights = c(3, 3, 1, 1))
plot(x = start(KLF5_gene),
     y = 1,
     pch = NA_integer_,
     xlim = c(start(KLF5_gene) - 5000, end(KLF5_gene) + 5000),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
lines(x = c(start(KLF5_gene) + 1000,
            start(KLF5_gene) + 6000),
      y = c(0.8, 0.8))
text("5 kb", x = start(KLF5_locus) + 3500, y = 0.85)
dev.off()
par(op)

pdf("../Graphs/KLF5-for-editing-2.pdf",
    width = 2.5,
    height = 5,
    useDingbats = FALSE)
for (i in seq_along(KLF5_enhancers)) {
  plotRegionWithATAC(cov = atac_cov_dmrs,
                     cov.region = NULL,
                     cov.addRegions = atac_NA_posvsBA9_pos_dmrs[
                       abs(atac_NA_posvsBA9_pos_dmrs$logFC) > 1],
                     cov.col = atac_cov_dmrs[[1]]$TISSUE_COLOR,
                     cov.lty = atac_cov_dmrs[[1]]$LTY,
                     cov.lwd = NULL,
                     cov.regionCol = alpha("red", 0.1),
                     k = 1001,
                     BSseq = BSseq_pos_dmrs,
                     region = KLF5_enhancers[i],
                     extend = 5000,
                     main = NULL,
                     addRegions = pos_dmrs,
                     annoTrack = list("KLF5-linked" = KLF5_enhancers,
                                      "BrainEnh" = Brain_enh_dmrs),
                     cex.anno = 1,
                     geneTrack = RefSeq_exons3_dmrs,
                     cex.gene = 1.5,
                     col = BSseq_pos_dmrs$Tissue_color,
                     lty = BSseq_pos_dmrs$lty,
                     lwd = NULL,
                     BSseqStat = NULL,
                     stat = "tstat.corrected",
                     stat.col = "black",
                     stat.lwd = 1,
                     stat.lty = 1,
                     stat.ylim = c(-8, 8),
                     mainWithWidth = TRUE,
                     regionCol = alpha("red", 0.1),
                     addTicks = TRUE,
                     addPoints = FALSE,
                     pointsMinCov = 5,
                     highlightMain = FALSE)
}

op <- par(no.readonly = TRUE)
opar <- par(mar = c(0, 4.1, 0, 0), oma = c(5, 0, 4, 2), mfrow = c(1, 1))
layout(matrix(1:4, ncol = 1), heights = c(3, 3, 1, 1))
plot(x = start(KLF5_enhancers[1]),
     y = 1,
     pch = NA_integer_,
     xlim = c(start(KLF5_enhancers[1]) - 5000, end(KLF5_enhancers[1]) + 5000),
     ylim = c(0, 1),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
lines(x = c(start(KLF5_enhancers[1]) + 2000,
            start(KLF5_enhancers[1]) + 3000),
      y = c(0.3, 0.3))
text("1 kb", x = start(KLF5_enhancers[1]) + 2500, y = 0.35)
dev.off()
par(op)

# Expression data
e <- filter(rna_atac_meth, gene_symbol == "KLF5")
e_df <- data_frame(
  sample = factor(ifelse(grepl("BA9", colnames(e$exp[[1]])), "BA9", "NAcc")),
  colour = ifelse(grepl("BA9", colnames(e$exp[[1]])), "deepskyblue", "chocolate1"),
  exp = as.vector(e$exp[[1]]))

pdf("../Graphs/KLF5-expression-for-editing.pdf", height = 5, width = 3)
stripchart(exp ~ sample,
           data = e_df,
           col = e_df$colour,
           vertical = TRUE,
           pch = 16,
           method = "jitter",
           yaxt = "n",
           xaxt = "n",
           ylab = "log2(exp)",
           frame.plot = FALSE)
axis(side = 2, at = c(3, 6))
axis(side = 1, at = c(1, 2), labels = c("BA9", "NAcc"))
lines(x = c(1.5, 1.5),
      y = c(mean(e$exp[[1]][, grepl("BA9", colnames(e$exp[[1]]))]),
            mean(e$exp[[1]][, grepl("BA9", colnames(e$exp[[1]]))]) +
              e$expLogFC))
text(round(e$expLogFC, 1),
     x = 1.7,
     y = mean(e$exp[[1]][, grepl("BA9", colnames(e$exp[[1]]))]) +
       e$expLogFC / 2)
dev.off()
