# Enhanced plotManyRegions() to include ATAC-seq data
# Peter Hickey
# 2016-10-13

library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(readr)

#-------------------------------------------------------------------------------
# Construct GPos objects spanning autosomes in hg19

# NOTE: GPos objects are currently limited to length == .Machine$integer.max so
#       have to split hg19 up by seqlevel
si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
sls <- paste0("chr", 1:22)
list_of_gpos <- lapply(sls, function(sl) {
  GPos(si[sl])
})
names(list_of_gpos) <- sls

#-------------------------------------------------------------------------------
# Load RleList of coverage vectors
fl <- list.files("ATAC-seq/extdata/flow-sorted-brain-atac/objects",
                 pattern = glob2rx("*rep*_cov.rds"),
                 full.names = TRUE)
list_of_cov <- lapply(fl,
                      readRDS)
sn <- gsub("_cov.rds", "", basename(fl))

#-------------------------------------------------------------------------------
# Load library sizes
list_of_ls <-
  readRDS("ATAC-seq/objects/library-size-flow-sorted-brain-ATAC-seq.rds")
list_of_ls_no_dup <- lapply(list_of_ls, "[[", "ls_no_dup")
ave_ls_no_dup <- mean(unlist(list_of_ls_no_dup))

#-------------------------------------------------------------------------------
# Construct SummarizedExperiment objects
# NOTE: assay is coverage scaled by library size scaled by 1,000,000 reads
# NOTE: list_of_se is a list of SummarizedExperiment object where the
#       @rowRanges slot of each SummarizedExperiment is a GPos object covering
#       a single seqlevel and the assay is a Rle for each sample of
#       (normalised) sequencing coverage

list_of_se <- lapply(sls, function(sl) {
  norm_cov <- DataFrame(Map(function(c, ls) {
    c / ls * 1000000
  }, c = lapply(list_of_cov, "[[", sl),
  ls = list_of_ls_no_dup))
  colnames(norm_cov) <- sn
  gpos <- list_of_gpos[[sl]]
  SummarizedExperiment(assays = SimpleList(norm_cov = norm_cov),
                       rowRanges = list_of_gpos[[sl]])
})
names(list_of_se) <- sls

# Add colData
cd <- readRDS("ATAC-seq/objects/colData-flow-sorted-brain-atac-seq.rds")
list_of_se <- lapply(list_of_se, function(se) {
  colData(se) <- cd
  se
})

saveRDS(list_of_se,
        "ATAC-seq/extdata/flow-sorted-brain-atac/objects/flow-sorted-brain-atac.cov.list_of_se.rds")

# Only keep pos samples
list_of_se <- lapply(list_of_se, function(se) {
  se[, se$NEUN == "pos"]
})

#-------------------------------------------------------------------------------
# Load regions

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/DMRs_BLOCKs_for_plotting.rda")
regions <- unlist(GRangesList(lapply(DMRs_BLOCKs_for_plotting, granges)))

#-------------------------------------------------------------------------------
# Load BSseq object

load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/BS.fit.small.sorted.somatic.all.rda")

# Only keep pos samples
BSseq <- BS.fit.small.sorted.somatic.all[, BS.fit.small.sorted.somatic.all$NeuN == "pos"]

#-------------------------------------------------------------------------------
# Load gene models
#

# TODO: Will likely change this object to use GENCODE gene models
load("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Objects_for_Paper/RefSeq_exons3.rds")

#-------------------------------------------------------------------------------
# Load brain-specific enhancers
#

load("Objects/Brain_enh.rda")

#-------------------------------------------------------------------------------
# Load ATAC-seq peaks
#

atac_peaks <- granges(readRDS("/dcl01/hansen/data/flow-sorted-brain-atac/objects/flow-sorted-brain-atac.overall.se.rds"))

#-------------------------------------------------------------------------------
# Load significantly differential ATAC-seq peaks
#

ATAC <- read_csv("/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/topTable.NA_posvsBA9_pos.ATAC-seq.csv")
colnames(ATAC) <- c("X","chr","start","end","width","strand","logFC","AveExpr","t","P.Value","adj.P.Val","B")
ATAC_gr <- GRanges(ATAC)
sig_ATAC_gr <- ATAC_gr[ATAC_gr$adj.P.Val <= 0.05 &
                         abs(ATAC_gr$logFC > 1)]

#-------------------------------------------------------------------------------
# Plotting functions

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

.plotCovData <- function(cov, region, extend, col, lty, lwd) {
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  cov <- subsetByOverlaps(cov, gr)

  ## Extract basic information
  sampleNames <- colnames(cov)
  names(sampleNames) <- sampleNames
  positions <- start(cov)
  coverage <- assay(cov)

  ## get col, lwd, lty
  colEtc <- .covGetCol(object = cov, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(min(coverage), max(coverage)), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "Normalised ATAC-seq coverage")
  axis(side = 2,
       at = (max(coverage) - min(coverage)) * c(0.2, 0.5, 0.8))

  sapply(1:ncol(cov), function(sampIdx) {
    .covPlotLines(positions, coverage[, sampIdx], col = colEtc$col[sampIdx],
                  lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx])
  })
}

plotGeneTrack <- bsseq:::plotGeneTrack
plotRegionWithATAC <- function(cov, cov.col = NULL, cov.lty = NULL,
                               cov.lwd = NULL,
                               BSseq, region = NULL, extend = 0, main = "",
                               addRegions = NULL,
                               annoTrack = NULL, cex.anno = 1,
                               geneTrack = NULL, cex.gene = 1.5,
                               col = NULL, lty = NULL, lwd = NULL,
                               BSseqStat = NULL, stat = "tstat.corrected",
                               stat.col = "black",
                               stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                               mainWithWidth = TRUE,
                               regionCol = alpha("red", 0.1), addTicks = TRUE,
                               addPoints = FALSE,
                               pointsMinCov = 5, highlightMain = FALSE) {

  opar <- par(mar = c(0,4.1,0,0), oma = c(5,0,4,2), mfrow = c(1,1))
  on.exit(par(opar))
  if (is.null(BSseqStat) && is.null(geneTrack)) {
    layout(matrix(1:3, ncol = 1), heights = c(2, 2,1))
  } else if (is.null(geneTrack)) {
    layout(matrix(1:4, ncol = 1), heights = c(2,2,2,1))
  } else {
    layout(matrix(1:5, ncol = 1), heights = c(2,2,2,1,0.3))
  }

  # Plot BSseq object
  bsseq:::.plotSmoothData(BSseq = BSseq, region = region, extend = extend,
                          addRegions = addRegions,
                          col = col, lty = lty, lwd = lwd, regionCol = regionCol,
                          addTicks = addTicks, addPoints = addPoints,
                          pointsMinCov = pointsMinCov,
                          highlightMain = highlightMain)
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)

  # Plot cov object
  .plotCovData(cov = cov, region = region, extend = extend, col = cov.col,
               lty = cov.lty, lwd = cov.lwd)

  if(!is.null(BSseqStat)) {
    BSseqStat <- subsetByOverlaps(BSseqStat, gr)
    if(is(BSseqStat, "BSseqTstat")) {
      stat.values <- getStats(BSseqStat)[, stat]
      stat.type <- "tstat"
    }
    if(is(BSseqStat, "BSseqStat")) {
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
         ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "", ylab = stat.type)
    axis(side = 2, at = c(-5,0,5))
    abline(h = 0, col = "grey60")
    .bsPlotLines(start(BSseqStat), stat.values, lty = stat.lty, col = stat.col, lwd = stat.lwd,
                 plotRange = c(start(gr), end(gr)))
  }

  if (!is.null(annoTrack)) {
    bsseq:::plotAnnoTrack(gr, annoTrack, cex.anno)
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

plotManyRegionsWithATAC <- function(list_of_se, cov.col = NULL, cov.lty = NULL,
                                    cov.lwd = NULL,
                                    BSseq, regions = NULL, extend = 0,
                                    main = "", addRegions = NULL,
                                    annoTrack = NULL, cex.anno = 1,
                                    geneTrack = NULL, cex.gene = 1.5, col = NULL,
                                    lty = NULL, lwd = NULL, BSseqStat = NULL,
                                    stat = "tstat.corrected",
                                    stat.col = "black", stat.lwd = 1,
                                    stat.lty = 1, stat.ylim = c(-8, 8),
                                    mainWithWidth = TRUE,
                                    regionCol = alpha("red", 0.1),
                                    addTicks = TRUE, addPoints = FALSE,
                                    pointsMinCov = 5, highlightMain = FALSE,
                                    verbose = TRUE) {

  cat("[plotManyRegionsWithATAC] preprocessing ...")
  if (!is.null(regions)) {
    if (is(regions, "data.frame")) {
      gr <- data.frame2GRanges(regions, keepColumns = FALSE)
    } else {
      gr <- regions
    }
    if (!is(gr, "GRanges")) {
      stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
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
    se <- subsetByOverlaps(list_of_se[[as.character(seqnames(gr[i]))]], gr[i])
    assay(se) <- do.call(cbind, lapply(assay(se), as.matrix))
    se
  }))
  if (!is.null(BSseqStat)) {
    BSseqStat <- subsetByOverlaps(BSseqStat, gr)
  }

  if (length(start(BSseq)) == 0) {
    stop("No overlap between BSseq data and regions")
  }
  if (!is.null(main) && length(main) != length(gr)) {
    main <- rep(main, length = length(gr))
  }
  cat("done\n")
  for (ii in seq_along(gr)) {
    if (verbose) {
      cat(sprintf("[plotManyRegionsWithATAC] plotting region %d (out of %d)\n",
                  ii, length(regions)))
    }
    plotRegionWithATAC(cov = cov, cov.col = cov.col, cov.lty = cov.lty,
                       cov.lwd = cov.lwd,
                       BSseq = BSseq, region = regions[ii, ],
                       extend = extend, col = col, lty = lty, lwd = lwd,
                       main = main[ii], BSseqStat = BSseqStat,
                       stat = stat, stat.col = stat.col, stat.lwd = stat.lwd,
                       stat.lty = stat.lty, stat.ylim = stat.ylim,
                       addRegions = addRegions, regionCol = regionCol,
                       mainWithWidth = mainWithWidth,
                       annoTrack = annoTrack, cex.anno = cex.anno,
                       geneTrack = geneTrack, cex.gene = cex.gene,
                       addTicks = addTicks, addPoints = addPoints,
                       pointsMinCov = pointsMinCov,
                       highlightMain = highlightMain)
  }
}

#-------------------------------------------------------------------------------
# Plot all NeuN+ samples

plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions[1],
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)

pdf("tmp-regions-pos.pdf")
plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions,
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Plot all NA-NeuN+ samples

plotManyRegionsWithATAC(list_of_se = list_of_se_NA_plot_regions,
                        cov.col = list_of_se_NA_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_NA_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_NA_plot_regions,
                        regions = regions[1],
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_NA_plot_regions$Tissue_color,
                        lty = BSseq_NA_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)

pdf("tmp-regions-NA-pos.pdf")
plotManyRegionsWithATAC(list_of_se = list_of_se_NA_plot_regions,
                        cov.col = list_of_se_NA_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_NA_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_NA_plot_regions,
                        regions = regions,
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_NA_plot_regions$Tissue_color,
                        lty = BSseq_NA_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Make objects to focus on just the regions to be used in figure
#

# NOTE: Arbitrarily adding 100kb either side of region
plot_regions <- resize(regions, 200000 + width(regions), fix = "center")

BSseq_plot_regions <- subsetByOverlaps(BSseq, plot_regions)
BSseq_NA_plot_regions <- BSseq_plot_regions[, BSseq_plot_regions$Tissue == "NA"]

list_of_se_plot_regions <- lapply(list_of_se, function(se) {
  subsetByOverlaps(se, plot_regions)
})
list_of_se_NA_plot_regions <- lapply(list_of_se_plot_regions, function(se) {
  se[, se$TISSUE == "NA"]
})
atac_peaks_plot_regions <- subsetByOverlaps(atac_peaks, plot_regions)
sig_ATAC_gr_plot_regions <- subsetByOverlaps(sig_ATAC_gr, plot_regions)
RefSeq_exons3_plot_regions <-
  RefSeq_exons3[queryHits(findOverlaps(makeGRangesFromDataFrame(RefSeq_exons3),
                                       plot_regions)), ]
Brain_enh_plot_regions <- subsetByOverlaps(Brain_enh, plot_regions)

save(regions,
     BSseq_plot_regions, BSseq_NA_plot_regions,
     list_of_se_plot_regions, list_of_se_NA_plot_regions,
     atac_peaks_plot_regions,
     sig_ATAC_gr_plot_regions,
     RefSeq_exons3_plot_regions,
     Brain_enh_plot_regions,
     file = "~/tmp-plotManyRegionsWithATAC.rda")

#-------------------------------------------------------------------------------
# Plotting experiments
#

# Transparency of cov
plotManyRegionsWithATAC(list_of_se = list_of_se_NA_plot_regions,
                        cov.col = scales::alpha(
                          list_of_se_NA_plot_regions[[1]]$TISSUE_COLOR, 1 / 6),
                        cov.lty = ifelse(
                          list_of_se_NA_plot_regions[[1]]$NEUN == "pos", 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_NA_plot_regions,
                        regions = regions[1],
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_NA_plot_regions$Tissue_color,
                        lty = BSseq_NA_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)

# Transparency of cov and overlap average cov
.plotCovData <- function(cov, region, extend, col, lty, lwd) {
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  cov <- subsetByOverlaps(cov, gr)

  ## Extract basic information
  sampleNames <- colnames(cov)
  names(sampleNames) <- sampleNames
  positions <- start(cov)
  coverage <- assay(cov)

  ## get col, lwd, lty
  colEtc <- .covGetCol(object = cov, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(min(coverage), max(coverage)), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "Normalised ATAC-seq coverage")
  axis(side = 2,
       at = (max(coverage) - min(coverage)) * c(0.2, 0.5, 0.8))

  sapply(1:ncol(cov), function(sampIdx) {
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

plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty = ifelse(
                          list_of_se_plot_regions[[1]]$NEUN == "pos", 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions[1],
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)

# Transparency of cov and overlap average cov with smoothing and fixed ylim
pdf("tmp-vary-k.pdf")
for (kk in c(1, 11, 101, 1001, 10001)) {
  .plotCovData <- function(cov, region, extend, col, lty, lwd, k = kk) {
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
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = c(0, 2), xlim = c(start(gr), end(gr)),
         xlab = "",
         ylab = "Normalised ATAC-seq coverage")
    axis(side = 2,
         at = c(0.4, 1, 1.6))

    sapply(1:ncol(cov), function(sampIdx) {
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

  plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                          cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                          cov.lty =
                            ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                   1, 2),
                          cov.lwd = NULL,
                          BSseq = BSseq_plot_regions,
                          regions = regions[width(regions) > 5000][1],
                          extend = 10000,
                          main = paste0("k = ", kk),
                          addRegions = regions,
                          annoTrack = list(ATAC = atac_peaks_plot_regions,
                                           sigATACwFC = sig_ATAC_gr_plot_regions,
                                           BrainEnh = Brain_enh_plot_regions),
                          cex.anno = 1,
                          geneTrack = RefSeq_exons3_plot_regions,
                          cex.gene = 1.5,
                          col = BSseq_plot_regions$Tissue_color,
                          lty = BSseq_plot_regions$lty,
                          lwd = NULL,
                          BSseqStat = NULL,
                          stat = "tstat.corrected",
                          stat.col = "black",
                          stat.lwd = 1,
                          stat.lty = 1,
                          stat.ylim = c(-8, 8),
                          mainWithWidth = TRUE,
                          regionCol = scales::alpha("red", 0.1),
                          addTicks = TRUE,
                          addPoints = FALSE,
                          pointsMinCov = 5,
                          highlightMain = FALSE,
                          verbose = TRUE)
}
dev.off()

# Transparency of cov and overlap average cov with fixed y-lim
.plotCovData <- function(cov, region, extend, col, lty, lwd) {
  gr <- bsseq:::.bsGetGr(BSseq, region, extend)
  cov <- subsetByOverlaps(cov, gr)

  ## Extract basic information
  sampleNames <- colnames(cov)
  names(sampleNames) <- sampleNames
  positions <- start(cov)
  coverage <- assay(cov)

  ## get col, lwd, lty
  colEtc <- .covGetCol(object = cov, col = col, lty = lty, lwd = lwd)

  ## The actual plotting
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 2), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "Normalised ATAC-seq coverage")
  axis(side = 2, at = c(0.4, 1, 1.6))

  sapply(1:ncol(cov), function(sampIdx) {
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

plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions[1],
                        extend = 10000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)

#-------------------------------------------------------------------------------
# Trying to figure out how to plot ATAC-seq in blocks
#

# Transparency of cov and overlap average cov without smoothing

pdf("tmp-blocks-k-1.pdf")
.plotCovData <- function(cov, region, extend, col, lty, lwd, k = 1) {
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
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 2), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "Normalised ATAC-seq coverage")
  axis(side = 2,
       at = c(0.4, 1, 1.6))

  sapply(1:ncol(cov), function(sampIdx) {
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
plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions[width(regions) > 5000],
                        extend = 10000,
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

# Transparency of cov and overlap average cov with runmed(k = 1001) smoothing
# and fixed ylim. Smoothing value based on experimentations in tmp-vary-k.pdf
# NOTE: Don't use loess smoothing because slow when the number of points is
#       large (> 1000)
pdf("tmp-blocks-k-1001.pdf")
.plotCovData <- function(cov, region, extend, col, lty, lwd, k = 1001) {
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
  plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 2), xlim = c(start(gr), end(gr)),
       xlab = "",
       ylab = "Normalised ATAC-seq coverage")
  axis(side = 2,
       at = c(0.4, 1, 1.6))

  sapply(1:ncol(cov), function(sampIdx) {
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
plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions[width(regions) > 5000],
                        extend = 10000,
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

#-------------------------------------------------------------------------------
# Plotting genes models
#

tx_plot_regions <- subsetByOverlaps(transcripts_by_gene, regions)

# > table(lengths(tx_plot_regions))
#
# 1  2  4  6 10 14 15 17 18
# 5  2  3  2  1  1  1  1  1

exons_plot_regions <- subsetByOverlaps(exons_by_transcript,
                                       resize(regions, width(regions) + 1000000,
                                              fix = "center"))

# NOTE: Abusing current plotting routine
# TODO: gene_id and gene_name are really tx_name and isoforms are set to 1
exons_df <- data.frame(chr = as.character(unlist(seqnames(exons_plot_regions),
                                                 use.names = FALSE)),
                       start = unlist(start(exons_plot_regions),
                                      use.names = FALSE),
                       end = unlist(end(exons_plot_regions), use.names = FALSE),
                       gene_id = rep(names(exons_plot_regions),
                                     lengths(exons_plot_regions)),
                       exon_number = unlist(
                         lapply(exons_plot_regions, function(x) x$exon_rank),
                         use.names = FALSE),
                       strand = as.character(unlist(strand(exons_plot_regions),
                                                    use.names = FALSE)),
                       gene_name = rep(names(exons_plot_regions),
                                       lengths(exons_plot_regions)),
                       isoforms = 1)

plotGeneTrack <- function(gr, geneTrack, cex) {
  geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
  ov <- findOverlaps(geneTrack_gr, gr)
  genes <- geneTrack[queryHits(ov), ]
  plot(start(gr), 1, type = "n", xaxt = "n",
       yaxt = "n", bty = "n",
       ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)), xlab = "",
       ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
  genes$offset <- match(genes$gene_id, unique(genes$gene_id))
  if (nrow(genes) > 0) {
    n_tx <- max(genes$offset)
    for (g in 1:nrow(genes)) {
      geneind2 <- which(geneTrack$gene_name == genes$gene_name[g])
      geneind2 <- geneind2[which(geneTrack$isoforms[geneind2] == 1)]
      direction <- unique(geneTrack$strand[geneind2])
      ES <- geneTrack$start[geneind2]
      EE <- geneTrack$end[geneind2]
      Exons <- cbind(ES, EE)
      lines(x = c(min(ES), max(EE)), y = c(genes$offset[g] / (n_tx + 2),
                                           genes$offset[g] / (n_tx + 2)))
      apply(Exons, 1, function(x) {
        polygon(c(x[1], x[2], x[2], x[1]),
                c(genes$offset[g] / (n_tx + 1), genes$offset[g] / (n_tx + 1),
                  genes$offset[g] / (n_tx + 3), genes$offset[g] / (n_tx + 3)),
                col = "darkgrey")
      })
      # text((max(start(gr), min(ES)) + min(end(gr), max(EE))) / 2,
      # 1.2, genes$gene_name[g], cex = cex)
    }
  }
}

pdf("tmp-gencode.pdf")
plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions,
                        extend = 20000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = exons_df,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

pdf("tmp-refseq.pdf")
plotGeneTrack <- bsseq:::plotGeneTrack
plotManyRegionsWithATAC(list_of_se = list_of_se_plot_regions,
                        cov.col = list_of_se_plot_regions[[1]]$TISSUE_COLOR,
                        cov.lty =
                          ifelse(list_of_se_plot_regions[[1]]$NEUN == "pos",
                                 1, 2),
                        cov.lwd = NULL,
                        BSseq = BSseq_plot_regions,
                        regions = regions,
                        extend = 20000,
                        main = "",
                        addRegions = regions,
                        annoTrack = list(ATAC = atac_peaks_plot_regions,
                                         sigATACwFC = sig_ATAC_gr_plot_regions,
                                         BrainEnh = Brain_enh_plot_regions),
                        cex.anno = 1,
                        geneTrack = RefSeq_exons3_plot_regions,
                        cex.gene = 1.5,
                        col = BSseq_plot_regions$Tissue_color,
                        lty = BSseq_plot_regions$lty,
                        lwd = NULL,
                        BSseqStat = NULL,
                        stat = "tstat.corrected",
                        stat.col = "black",
                        stat.lwd = 1,
                        stat.lty = 1,
                        stat.ylim = c(-8, 8),
                        mainWithWidth = TRUE,
                        regionCol = scales::alpha("red", 0.1),
                        addTicks = TRUE,
                        addPoints = FALSE,
                        pointsMinCov = 5,
                        highlightMain = FALSE,
                        verbose = TRUE)
dev.off()

# NOTE: Need to carefully check that gene is plotted correctly, e.g., RBFOX3
#       isn't plotted because while it contains a DMR, there are no exons
#       within `extend = 10000` bp of the DMR (but does work for
#       `extend = 20000`)
