# Convert all bulk ATAC-seq data from flow sorted brain GTEx reads to
# SummarizedExperiment objects. ATAC-seq reads are counted based on peaks
# called in bulk ATAC-seq data from various groupings of the data
# Also create a HDF5Array-backed SummarizedExperiment of sequencing coverage
# Peter Hickey
# 2017-06-30

library(dplyr)
library(purrr)
library(stringr)
library(SummarizedExperiment)
library(Matrix)
library(rtracklayer)
library(GenomicAlignments)
library(atacr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DelayedArray)

options("mc.cores" = 27)

extdir <- "../extdata"

# Load colData
ATAC_CD <- readRDS("../objects/colData-flow-sorted-brain-atac-seq.rds")

grs <- list.files(file.path(extdir, "flow-sorted-brain-atac", "objects"),
                  "^.*_gr\\.rds$",
                  full.names = TRUE)
# Only keep those samples from the main project
grs <- grs[sapply(rownames(ATAC_CD), grep, grs)]

# ------------------------------------------------------------------------------
# Load GRanges objects from make-GRanges.R and retain only non-duplicate reads
#

list_of_gr <- mclapply(grs, function(x) {
  gr <- readRDS(x)
  # NOTE: Drop mcols to reduce object size
  granges(gr)[!gr$isDuplicate]
})

# NOTE: Sometimes hit "long vectors" issue because some GRanges objects are so
#       long. So read those long ones in serially to avoid this error.
i <- which(!sapply(list_of_gr, function(x) is(x, "GRanges")))
while (!all(sapply(list_of_gr, function(x) is(x, "GRanges")))) {
  message("i = ", paste0(i, collapse = ", "))
  i <- which(!sapply(list_of_gr, function(x) is(x, "GRanges")))
  list_of_gr[i] <- lapply(grs[i], function(x) {
    gr <- readRDS(x)
    # NOTE: Drop mcols to reduce object size
    granges(gr)[!gr$isDuplicate]
  })
}
names(list_of_gr) <- gsub("_gr.rds", "", basename(grs))

# ------------------------------------------------------------------------------
# Import and construct regions
#

groups <- c("NA-pos", "NA-neg", "BA9-pos", "BA9-neg",
            "pos", "neg", "overall")
list_of_summits <- mclapply(groups, function(group) {
  summits <- import(file.path(extdir, "flow-sorted-brain-atac", "data",
                              "macs2",
                              paste0("flow-sorted-brain-atac.",
                                     group, "_summits.bed")),
                    genome = "hg19")
  summits <- sortSeqlevels(summits)
  summits
})
names(list_of_summits) <- groups

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
list_of_reduced_narrowPeaks <- mclapply(groups, function(group) {
  narrowPeaks <- import(file.path(extdir, "flow-sorted-brain-atac", "data",
                                  "macs2",
                                  paste0("flow-sorted-brain-atac.",
                                         group, "_peaks.narrowPeak")),
                        genome = "hg19",
                        extraCols = extraCols_narrowPeak,
                        format = "BED")
  reduce(sortSeqlevels(narrowPeaks))
})
names(list_of_reduced_narrowPeaks) <- paste0(groups, "_narrowPeak_reduced")
list_of_reduced_narrowPeaks <-
  c(list_of_reduced_narrowPeaks,
    list(union_narrowPeak_reduced =
           Reduce(union,
                  list_of_reduced_narrowPeaks[c("NA-pos_narrowPeak_reduced",
                                                "NA-neg_narrowPeak_reduced",
                                                "BA9-pos_narrowPeak_reduced",
                                                "BA9-neg_narrowPeak_reduced")]
                  ),
         union_pos_narrowPeak_reduced =
           Reduce(union,
                  list_of_reduced_narrowPeaks[c("NA-pos_narrowPeak_reduced",
                                                "BA9-pos_narrowPeak_reduced")])
         ))

list_of_regions <- c(list_of_reduced_narrowPeaks,
                     lapply(list_of_summits, function(x) {
                       reduce(resize(x, width = 500, fix = "center"))
                     }))

# ------------------------------------------------------------------------------
# Get known 'bad' regions
#

# Get the ENCODE mappability consensus blacklist (used to filter out spurious
# peaks by Buenrostro et al. 2015)
# NOTE: This should work but doesn't; reported to bioc-support
#     (https://support.bioconductor.org/p/80040/)
# encode_blacklist <- import(BEDFile("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
eb_dir <- tempdir()
eb_path_hg19 <-
  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
download.file(url = eb_path_hg19,
              destfile = file.path(eb_dir, basename(eb_path_hg19)))
encode_blacklist_hg19 <- import(file.path(eb_dir, basename(eb_path_hg19)),
                                genome = "hg19")

# Load the hg19 blacklist created by Jason Buenrostro
jdb_blacklist_path_hg19 <-
  "/dcl01/hansen/data/atac-seq-blacklists/JDB_blacklist.hg19.bed"
jdb_blacklist_hg19 <- import(jdb_blacklist_path_hg19, genome = "hg19")

# Combine into common blacklist
blacklist <- sort(reduce(c(granges(encode_blacklist_hg19),
                           jdb_blacklist_hg19)))

# ------------------------------------------------------------------------------
# Count reads to regions
#

# NOTE: This may need to be run 'manually' depending on the whims of the gods'
lapply(names(list_of_regions), function(name) {
  message(name)
  # NOTE: reduce() so that peaks are disjoint
  regions <- subsetByOverlaps(list_of_regions[[name]], blacklist, invert = TRUE)
  list_of_se <- mclapply(names(list_of_gr), function(sampleName) {
    gr <- list_of_gr[[sampleName]]
    se <- summarizeOverlaps(features = regions, reads = gr)
    colnames(se) <- sampleName
    se
  }, mc.cores = getOption("mc.cores"))
  se <- do.call(cbind, list_of_se)
  colData(se) <- as(combine(as.data.frame(colData(se)),
                            as.data.frame(ATAC_CD)),
                    "DataFrame")

  saveRDS(se, file.path(extdir, "flow-sorted-brain-atac", "objects",
                        paste0("flow-sorted-brain-atac.", name, ".se.rds")))
})

### ============================================================================
### HDF5Array-backed SummarizedExperiment(s) of ATAC-seq coverage
###

# NOTE: Cannot make a single genome-wide SummarizedExperiment due to long
#       vector issue. Instead, make a per-base SummarizedExperiment for each
#       chromosome with 3 assays:
#         1. raw coverage
#         2. cpm
#         3. log2-cpm
#       Only do this for chr1-22, chrX, and chrY
sample_names <- rownames(ATAC_CD)
names(sample_names) <- sample_names
list_of_cov <- mclapply(sample_names, function(this_sample) {
  readRDS(file.path(extdir, "flow-sorted-brain-atac", "objects",
                    paste0(this_sample, "_cov.rds")))
}, mc.cores = options("mc.cores"))
lib_sizes <- mclapply(list_of_cov, function(cov) {
  sum(as.numeric(sum(cov)))
}, mc.cores = options("mc.cores"))

seqlevels <- paste0("chr", c(1:22, c("X", "Y")))
mclapply(seqlevels, function(seqlevel) {
  # Raw data for this seqlevel
  list_of_raw <- lapply(list_of_cov, "[[", seqlevel)
  # Scaled data for this seqlevel
  list_of_cpm <- mapply(function(raw, ls) {
    values <- cpm(x = runValue(raw), lib.size = ls, log = FALSE)
    Rle(values, runLength(raw))
  }, raw = list_of_raw, ls = lib_sizes)
  list_of_log2cpm <- mapply(function(raw, ls) {
    values <- cpm(x = runValue(raw), lib.size = ls, log = TRUE)
    Rle(values, runLength(raw))
  }, raw = list_of_raw, ls = lib_sizes)

  # Construct assays: raw, cpm, log2cpm
  # NOTE: Have to construct assays as a DelayedMatrix with cbind()-ed RleMatrix
  #       columns instead of as a single RleMatrix. This is because in the
  #       latter case I run into issues when trying to index into the result.
  raw_rlearray <- do.call(
    cbind, lapply(names(list_of_raw), function(this_sample) {
      RleArray(rle = list_of_raw[[this_sample]],
               dim = c(length(list_of_raw[[this_sample]]), 1),
               dimnames = list(NULL, this_sample))
    }))
  cpm_rlearray <- do.call(
    cbind, lapply(names(list_of_cpm), function(this_sample) {
      RleArray(rle = list_of_cpm[[this_sample]],
               dim = c(length(list_of_cpm[[this_sample]]), 1),
               dimnames = list(NULL, this_sample))
    }))
  log2cpm_rlearray <- do.call(
    cbind, lapply(names(list_of_log2cpm), function(this_sample) {
      RleArray(rle = list_of_log2cpm[[this_sample]],
               dim = c(length(list_of_log2cpm[[this_sample]]), 1),
               dimnames = list(NULL, this_sample))
    }))

  # Construct SummarizedExperiment
  se <- SummarizedExperiment(
    rowRanges = GPos(seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevel]),
    assays = SimpleList(raw = raw_rlearray,
                        cpm = cpm_rlearray,
                        log2cpm = log2cpm_rlearray),
    colData = ATAC_CD)

  # Save as HDF5Array-backed SummarizedExperiment
  saveHDF5SummarizedExperiment(
    x = se,
    dir = file.path(extdir, "flow-sorted-brain-atac", "objects",
                    paste0("flow-sorted-brain-atac.", seqlevel,
                           ".coverage.SummarizedExperiment")),
    verbose = TRUE)
}, mc.cores = options("mc.cores"))

### ============================================================================
### HDF5Array-backed SummarizedExperiment of (scaled) coverage +/- 1500bp
### around midpoints of DMRs and DARs (or the full DMR/DAR, whichever is larger)
###

#-------------------------------------------------------------------------------
# Features (DMRs and DARs)
#

window <- 3000
load("../../Objects/All_Annotated_DMRs_GRanges.rda")
dmrs <- c(Annotated_NEG_DMRs_gr, Annotated_POS_DMRs_gr,
          Annotated_POSvNEG_DMRs_gr, Annotated_Unsorted_DMRs_gr,
          ignore.mcols = TRUE)
dmrs <- reduce(resize(dmrs, pmax(window, width(dmrs)), fix = "center"))
load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
daps <- lapply(atac[-1], function(x) granges(x[x$adj.P.Val < 0.05 &
                                                 abs(x$logFC) > 1, ]))
dars <- c(dars_pos, dars_pos_vs_neg, dars_neg, ignore.mcols = TRUE)
dars <- reduce(resize(dars, pmax(window, width(dars)), fix = "center"))

features <- reduce(c(dmrs, dars, ignore.mcols = TRUE))
list_of_features <- split(features, seqnames(features))

#-------------------------------------------------------------------------------
# A SummarizedExperiment of ATAC-seq cpm at features
#

seqlevels <- paste0("chr", 1:22)
names(seqlevels) <- seqlevels
paths_to_se <- file.path(
  extdir, "flow-sorted-brain-atac", "objects",
  paste0("flow-sorted-brain-atac.", seqlevels,
         ".coverage.SummarizedExperiment"))
names(paths_to_se) <- seqlevels

list_of_se <- lapply(names(list_of_features), function(sl) {
  path <- paths_to_se[[sl]]
  se <- loadHDF5SummarizedExperiment(path)
  # NOTE: Only retain cpm assay
  assays(se) <- assays(se)["cpm"]
  subsetByOverlaps(query = se,
                   subject = list_of_features[[sl]])
})
se <- do.call(rbind, list_of_se)
options("DelayedArray.block.size" = 45000000L)
saveHDF5SummarizedExperiment(
  x = se,
  dir = file.path(extdir, "flow-sorted-brain-atac", "objects",
                  "coverage_near_dmrs_and_dars.SummarizedExperiment"),
  verbose = TRUE)
