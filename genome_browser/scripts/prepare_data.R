# Prepare WGBS, ATAC-seq, and RNA-seq data for use in a genome browser
# Peter Hickey
# 2018-03-19

# ==============================================================================
# Setup
#

library(rtracklayer)
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DelayedArray)
library(edgeR)
library(GenomicAlignments)
library(readr)
library(matrixStats)

extdir <- "../extdata"

# NOTE: Don't increase mc.cores (export.bw() craps out)
options("mc.cores" = 8)

load(file.path("..", "..", "integrating-dmrs-dars-and-degs", "objects",
               "assays-and-features.rda"))
load(file.path("..", "..", "Objects", "Annotated_POSvNEG_BLOCKs_GRanges.rda"))
load("../../Objects/All_BLOCK_POS_DMRs_fwer50.rda")

### ============================================================================
### Construct bigWig files
###

# ==============================================================================
# mCG (small smooth)
#

CG_BSseq_unsorted <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.unsorted.fit.small.somatic.all"))

# Per-sample mCG
mclapply(colnames(CG_BSseq_unsorted), function(this_sample) {
  message("Converting ", this_sample)
  gr <- rowRanges(CG_BSseq_unsorted)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  gr$score <- as.vector(getMeth(CG_BSseq_unsorted[, this_sample]))
  this_sample <- gsub("NA", "NAcc", this_sample)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(this_sample, "_unsorted.small_smooth.mCG.bw")))
}, mc.cores = options("mc.cores"))

# Per-condition average mCG
conditions <- c("BA9", "BA24", "HC", "NA", "caudate")
mclapply(conditions, function(condition) {
  message("Converting ", condition)
  gr <- rowRanges(CG_BSseq_unsorted)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  # NOTE: Realizing data to speed up subsequent call to rowMeans2() at cost of
  #       larger memory footprint
  meth <- as.array(getMeth(
    CG_BSseq_unsorted[, grep(condition, colnames(CG_BSseq_unsorted))]))
  gr$score <- rowMeans2(meth, na.rm = TRUE)
  condition <- gsub("NA", "NAcc", condition)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(condition, "_unsorted.small_smooth.mCG.bw")))
})

CG_BSseq_small_smooth <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.small.sorted.somatic.all"))

# Per-sample mCG
mclapply(colnames(CG_BSseq_small_smooth), function(this_sample) {
  message("Converting ", this_sample)
  gr <- rowRanges(CG_BSseq_small_smooth)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  gr$score <- as.vector(getMeth(CG_BSseq_small_smooth[, this_sample]))
  this_sample <- gsub("NA", "NAcc", this_sample)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(this_sample, ".small_smooth.mCG.bw")))
}, mc.cores = options("mc.cores"))

# Per-condition average mCG
conditions <- apply(
  X = expand.grid(c("BA9", "BA24", "HC", "NA"), c("_neg", "_pos")),
  FUN = paste,
  MARGIN = 1,
  collapse = "")
mclapply(conditions, function(condition) {
  message("Converting ", condition)
  gr <- rowRanges(CG_BSseq_small_smooth)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  # NOTE: Realizing data to speed up subsequent call to rowMeans2() at cost of
  #       larger memory footprint
  meth <- as.array(getMeth(
    CG_BSseq_small_smooth[, grep(condition, colnames(CG_BSseq_small_smooth))]))
  gr$score <- rowMeans2(meth, na.rm = TRUE)
  condition <- gsub("NA", "NAcc", condition)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(condition, ".small_smooth.mCG.bw")))
}, mc.cores = options("mc.cores"))

# ==============================================================================
# mCG (large smooth)
#

# NOTE: No large smooth for unsorted data

CG_BSseq_large_smooth <- loadHDF5SummarizedExperiment(
  file.path(extdir, "flow-sorted-brain-wgbs", "objects",
            "BS.fit.large.sorted.somatic.all"))

# Per-sample mCG
mclapply(colnames(CG_BSseq_large_smooth), function(this_sample) {
  message("Converting ", this_sample)
  gr <- rowRanges(CG_BSseq_large_smooth)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  gr$score <- as.vector(getMeth(CG_BSseq_large_smooth[, this_sample]))
  this_sample <- gsub("NA", "NAcc", this_sample)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(this_sample, ".large_smooth.mCG.bw")))
}, mc.cores = options("mc.cores"))

# Per-condition average mCG
conditions <- apply(
  X = expand.grid(c("BA9", "BA24", "HC", "NA"), c("_neg", "_pos")),
  FUN = paste,
  MARGIN = 1,
  collapse = "")
mclapply(conditions, function(condition) {
  message("Converting ", condition)
  gr <- rowRanges(CG_BSseq_large_smooth)
  seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
  # NOTE: Realizing data to speed up subsequent call to rowMeans2() at cost of
  #       larger memory footprint
  meth <- as.array(getMeth(
    CG_BSseq_large_smooth[, grep(condition, colnames(CG_BSseq_large_smooth))]))
  gr$score <- rowMeans2(meth, na.rm = TRUE)
  condition <- gsub("NA", "NAcc", condition)
  export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                       paste0(condition, ".large_smooth.mCG.bw")))
})

# TODO: Test out 'MethylC' track from WashU browser

# ==============================================================================
# mCH (small smooth)
#

# NOTE: Only NeuN+ samples have mCH data
CH_BSseq_names <- c("pos_CA", "pos_CT", "neg_CA", "neg_CT")
names(CH_BSseq_names) <- CH_BSseq_names
list_of_CH_BSseq <- lapply(CH_BSseq_names, function(n) {
  BSseq <- loadHDF5SummarizedExperiment(
    dir = file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                    paste0(n, "_small-flow-sorted-brain-wgbs")))
  BSseq <- BSseq[, grep("pos", colnames(BSseq))]
  colnames(BSseq) <- gsub("NA", "NAcc", colnames(BSseq))
})

# Per-sample mCH for all strands and contexts
lapply(names(list_of_CH_BSseq), function(n) {
  BSseq <- list_of_CH_BSseq[[n]]
  mclapply(colnames(BSseq), function(this_sample) {
    message("Converting ", this_sample, " (", n, ")")
    gr <- rowRanges(BSseq)
    seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
    gr$score <- as.vector(getMeth(BSseq[, this_sample]))
    this_sample <- gsub("NA", "NAcc", this_sample)
    export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                         paste0(this_sample, ".m", n, ".bw")))
  }, mc.cores = getOption("mc.cores"))
})

# Per-condition average mCH for all strands and contexts
conditions <- paste0(c("BA9", "BA24", "HC", "NA"), "_pos")
mclapply(conditions, function(condition) {
  lapply(names(list_of_CH_BSseq), function(n) {
    BSseq <- list_of_CH_BSseq[[n]]
    message("Converting ", condition, " (", n, ")")
    gr <- rowRanges(BSseq)
    seqinfo(gr) <- intersect(seqinfo(gr), seqinfo(BSgenome.Hsapiens.UCSC.hg19))
    # NOTE: Realizing data in memory to speed up subsequent call to rowMeans2()
    #       at cost of larger memory footprint
    meth <- as.array(getMeth(BSseq[, grep(condition, colnames(BSseq))]))
    gr$score <- rowMeans2(meth, na.rm = TRUE)
    condition <- gsub("NA", "NAcc", condition)
    export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                         paste0(condition, ".m", n, ".bw")))
  })
}, mc.cores = getOption("mc.cores"))

# Per-condition average mCH for all strands and contexts where negative
# strand loci have their methylation level negated and data from the two
# strands are combined into a single bigWig file
mclapply(conditions, function(condition) {
  condition <- gsub("NA", "NAcc", condition)
  mclapply(c("CA", "CT"), function(context) {
    message(condition, " (", context, "): Loading pos strand")
    pos_strand <- import(
      file.path(extdir, "flow-sorted-brain-wgbs", "data",
                "bigWig", paste0(condition, ".pos_m", context, ".bw")))
    message(condition, " (", context, "): Loading neg strand")
    neg_strand <- import(
      file.path(extdir, "flow-sorted-brain-wgbs", "data",
                "bigWig", paste0(condition, ".neg_m", context, ".bw")))
    # Negate methylation levels on negative strand
    neg_strand$score <- -neg_strand$score

    message(condition, " (", context, "): Combining pos and neg strand")
    gr <- sort(c(pos_strand, neg_strand))
    message(condition, " (", context, "): Exporting to bigWig")
    export(gr, file.path(extdir, "flow-sorted-brain-wgbs", "data", "bigWig",
                         paste0(condition, ".m", context, ".bw")))
  }, mc.cores = 2)
}, mc.cores = 4)

# ==============================================================================
# ATAC-seq
#

ATAC_CD <-
  readRDS("../../ATAC-seq/objects/colData-flow-sorted-brain-atac-seq.rds")

# ------------------------------------------------------------------------------
# Rle coverage vector to bigWig(s): cpm
#

# NOTE: Only retain chr1-22, chrX, and chrY
seqlevels <- paste0("chr", c(1:22, "X", "Y"))
# Per-sample ATAC-seq coverage
mclapply(rownames(ATAC_CD), function(this_sample) {
  # Load data
  message("Loading ", this_sample)
  cov <- readRDS(file.path(extdir, "flow-sorted-brain-atac", "objects",
                           paste0(this_sample, "_cov.rds")))
  this_sample <- gsub("-", "_", gsub("NA", "NAcc", this_sample))
  cov <- cov[seqlevels]

  # Raw coverage
  message("Converting ", this_sample, " (raw)")
  gr <- as(cov, "GRanges")

  # Counts per-million
  message("Converting ", this_sample, " (cpm)")
  gr$score <- as.vector(cpm(gr$score, lib.size = sum(as.numeric(sum(cov)))))
  export(gr_cpm,
         file.path(extdir, "flow-sorted-brain-atac", "data", "bigWig",
                   paste0(this_sample, ".ATAC-seq.cpm.bw")))
}, mc.cores = getOption("mc.cores"))

# Per-condition average ATAC-seq coverage
# NOTE: Only using rep1 data
conditions <- c("BA9_neg", "BA9_pos", "NA_neg", "NA_pos")
sample_names <- rownames(ATAC_CD)
sample_names <- grep("rep1", sample_names, value = TRUE)
names(sample_names) <- sample_names
list_of_cov <- mclapply(sample_names, function(this_sample) {
  readRDS(file.path(extdir, "flow-sorted-brain-atac", "objects",
                    paste0(this_sample, "_cov.rds")))
}, mc.cores = options("mc.cores"))
lib_sizes <- mclapply(list_of_cov, function(cov) {
  sum(as.numeric(sum(cov)))
}, mc.cores = options("mc.cores"))

lapply(conditions, function(condition) {
  # NOTE: Bit of fussing with condition due to it being written differently in
  #       different file names
  condition <- gsub("_", "-", condition)
  list_of_cov <- list_of_cov[grep(condition, names(list_of_cov))]
  lib_sizes <- lib_sizes[grep(condition, names(lib_sizes))]
  condition <- gsub("NA", "NAcc", gsub("-", "_", condition))
  message("Converting ", condition)
  seqlevels <- paste0("chr", c(1:22, c("X", "Y")))

  # Loop over each seqlevel and construct GRanges object
  list_of_gr_cpm_by_group <- mclapply(seqlevels, function(seqlevel) {
    message(seqlevel, ": Getting raw covereage")
    # Raw data for this seqlevel
    list_of_raw <- lapply(list_of_cov, "[[", seqlevel)
    raw <- do.call(cbind, lapply(list_of_raw, as.vector))

    # cpm data for this seqlevel
    # NOTE: zero rows map to zero cpm, so no need to use these
    message(seqlevel, ": Computing cpm")
    all_zero <- matrixStats::rowAlls(raw, value = 0)
    cpm <- matrix(0, nrow = nrow(raw), ncol = 1)
    cpm[!all_zero] <- cpmByGroup(raw[!all_zero, ],
                                 group = rep(1, ncol(raw)))

    message(seqlevel, ": Constructing GRanges")
    gr_cpm_by_group <- as(setNames(RleList(Rle(cpm)), seqlevel),
                          "GRanges")

    # Return GRanges instances
    message(seqlevel, ": Finished!")
    gr_cpm_by_group
  }, mc.cores = options("mc.cores"))

  # Construct GRanges instances
  gr_cpm_by_group <- do.call(c, list_of_gr_cpm_by_group)
  seqinfo(gr_cpm_by_group) <- intersect(seqinfo(gr_cpm_by_group),
                                        seqinfo(BSgenome.Hsapiens.UCSC.hg19))

  # Export as bigWig files
  message("Exporting ", condition, " (cpm)")
  export(gr_cpm_by_group,
         file.path(extdir, "flow-sorted-brain-atac", "data", "bigWig",
                   paste0(condition, ".ATAC-seq.cpm.bw")))
})

# ------------------------------------------------------------------------------
# RNA-seq
#

# NOTE: Can't use Salmon files since no coverage info, so instead using HiSat2
#       alignments.
# TODO: Is using HiSat2 alignments kosher? Could instead re-run Salmon with the
#       (new) --writeMappings option
# TODO: HiSat2 bigWig look very noisy (because of unspliced mRNA?). Probably
#       also due to variable y-lim scaling (want to ensure ymax >> 1 cpm)
RNA_CD <- readRDS(file.path(extdir, "flow-sorted-brain-rna-seq", "objects",
                            "colData-flow-sorted-brain-rna-seq.rds"))
seqlevels <- paste0("chr", c(1:22, "X", "Y"))

mclapply(rownames(RNA_CD), function(this_sample) {
  message("Parsing BAM ", this_sample)
  bam <- BamFile(file = file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                                  "bam", paste0(this_sample, ".sorted.bam")))
  cov <- coverage(bam)[seqlevels]
  this_sample <- gsub("-", "_", gsub("NA", "NAcc", this_sample))

  # Raw coverage
  message("Converting ", this_sample, " (raw)")
  gr <- as(cov, "GRanges")

  # Counts per-million
  message("Converting ", this_sample, " (cpm)")
  gr$score <- as.vector(cpm(gr$score, lib.size = sum(as.numeric(sum(cov)))))
  export(gr,
         file.path(extdir, "flow-sorted-brain-rna-seq", "data", "bigWig",
                   paste0(this_sample, ".RNA-seq.cpm.bw")))
}, mc.cores = options("mc.cores"))

# Per-condition average RNA-seq coverage
conditions <- c("BA9_neg", "BA9_pos", "NA_neg", "NA_pos")
sample_names <- rownames(RNA_CD)
names(sample_names) <- sample_names
list_of_cov <- mclapply(sample_names, function(this_sample) {
  message("Parsing BAM ", this_sample)
  bam <- BamFile(file = file.path(extdir, "flow-sorted-brain-rna-seq", "data",
                                  "bam", paste0(this_sample, ".sorted.bam")))
  cov <- coverage(bam)[seqlevels]
}, mc.cores = options("mc.cores"))
lib_sizes <- mclapply(list_of_cov, function(cov) {
  sum(as.numeric(sum(cov)))
}, mc.cores = options("mc.cores"))

lapply(conditions, function(condition) {
  # NOTE: Bit of fussing with condition due to it being written differently in
  #       different file names
  condition <- gsub("_", "-", condition)
  list_of_cov <- list_of_cov[grep(condition, names(list_of_cov))]
  lib_sizes <- lib_sizes[grep(condition, names(lib_sizes))]
  condition <- gsub("NA", "NAcc", gsub("-", "_", condition))
  message("Converting ", condition)

  # Loop over each seqlevel and construct GRanges object
  list_of_gr_cpm_by_group <- mclapply(seqlevels, function(seqlevel) {
    message(seqlevel, ": Getting raw covereage")
    # Raw data for this seqlevel
    list_of_raw <- lapply(list_of_cov, "[[", seqlevel)
    raw <- do.call(cbind, lapply(list_of_raw, as.vector))

    # cpm data for this seqlevel
    # NOTE: zero rows map to zero cpm, so no need to use these
    message(seqlevel, ": Computing cpm")
    all_zero <- matrixStats::rowAlls(raw, value = 0)
    cpm <- matrix(0, nrow = nrow(raw), ncol = 1)
    cpm[!all_zero] <- cpmByGroup(raw[!all_zero, ],
                                 group = rep(1, ncol(raw)))

    message(seqlevel, ": Constructing GRanges")
    gr_cpm_by_group <- as(setNames(RleList(Rle(cpm)), seqlevel),
                          "GRanges")

    # Return GRanges instances
    message(seqlevel, ": Finished!")
    gr_cpm_by_group
  }, mc.cores = options("mc.cores"))

  # Construct GRanges instances
  gr_cpm_by_group <- do.call(c, list_of_gr_cpm_by_group)
  seqinfo(gr_cpm_by_group) <- intersect(seqinfo(gr_cpm_by_group),
                                        seqinfo(BSgenome.Hsapiens.UCSC.hg19))

  # Export as bigWig files
  message("Exporting ", condition, " (cpm)")
  export(gr_cpm_by_group,
         file.path(extdir, "flow-sorted-brain-rna-seq", "data", "bigWig",
                   paste0(condition, ".RNA-seq.cpm.bw")))
})

### ============================================================================
### Construct BED files: (DMRs, DARs, blocks)
###

# ==============================================================================
# chrom.sizes file needed for constructing a bigBED file from a BED file
#

si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
write.table(as.data.frame(si)[, 1, drop = FALSE],
            file.path(extdir, "hg19.chrom.sizes"),
            col.names = FALSE,
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")
# NOTE: BED to bigBED requies the chromosomes are sorted lexicographically
si_bed <- si[sort(seqnames(si))]

# ==============================================================================
# CG-DMRs
#

seqlevels(Annotated_Unsorted_DMRs_gr) <-
  seqlevels(intersect(si_bed, seqinfo(Annotated_Unsorted_DMRs_gr)))
seqinfo(Annotated_Unsorted_DMRs_gr) <-
  merge(seqinfo(Annotated_Unsorted_DMRs_gr), si_bed)
export(
  sort(Annotated_Unsorted_DMRs_gr),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-DMRs.unsorted.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-DMRs.unsorted.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-DMRs.unsorted.bb")))

seqlevels(Annotated_POSvNEG_DMRs_gr) <-
  seqlevels(intersect(si_bed, seqinfo(Annotated_POSvNEG_DMRs_gr)))
seqinfo(Annotated_POSvNEG_DMRs_gr) <-
  merge(seqinfo(Annotated_POSvNEG_DMRs_gr), si_bed)
export(
  sort(Annotated_POSvNEG_DMRs_gr),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-DMRs.pos_vs_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-DMRs.pos_vs_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-DMRs.pos_vs_neg.bb")))

seqlevels(Annotated_POS_DMRs_gr) <-
  seqlevels(intersect(si_bed, seqinfo(Annotated_POS_DMRs_gr)))
seqinfo(Annotated_POS_DMRs_gr) <-
  merge(seqinfo(Annotated_POS_DMRs_gr), si_bed)
export(
  sort(Annotated_POS_DMRs_gr),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-DMRs.pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-DMRs.pos.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-DMRs.pos.bb")))

seqlevels(Annotated_NEG_DMRs_gr) <-
  seqlevels(intersect(si_bed, seqinfo(Annotated_NEG_DMRs_gr)))
seqinfo(Annotated_NEG_DMRs_gr) <-
  merge(seqinfo(Annotated_NEG_DMRs_gr), si_bed)
export(
  sort(Annotated_NEG_DMRs_gr),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-DMRs.neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-DMRs.neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-DMRs.neg.bb")))

seqlevels(dmrs_NAvsBA9pos) <-
  seqlevels(intersect(si_bed, seqinfo(dmrs_NAvsBA9pos)))
seqinfo(dmrs_NAvsBA9pos) <-
  merge(seqinfo(dmrs_NAvsBA9pos), si_bed)
export(
  sort(dmrs_NAvsBA9pos),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-DMRs.NAcc_pos_vs_BA9_pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-DMRs.NAcc_pos_vs_BA9_pos.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-DMRs.NAcc_pos_vs_BA9_pos.bb")))

# ==============================================================================
# CG-blocks
#

seqlevels(Annotated_POSvNEG_BLOCKs_gr) <-
  seqlevels(intersect(si_bed, seqinfo(Annotated_POSvNEG_BLOCKs_gr)))
seqinfo(Annotated_POSvNEG_BLOCKs_gr) <-
  merge(seqinfo(Annotated_POSvNEG_BLOCKs_gr), si_bed)
export(
  sort(Annotated_POSvNEG_BLOCKs_gr),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-blocks.pos_vs_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-blocks.pos_vs_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-blocks.pos_vs_neg.bb")))

pos_blocks <- GRanges(sig_block_dmrs)
seqlevels(pos_blocks) <- seqlevels(intersect(si_bed, seqinfo(pos_blocks)))
seqinfo(pos_blocks) <- merge(seqinfo(pos_blocks), si_bed)
export(
  sort(pos_blocks),
  file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
            "CG-blocks.pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                        "CG-blocks.pos.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "CG-blocks.pos.bb")))

# NOTE: We never did NEG blocks

# ==============================================================================
# CH-DMRs
#

list_of_candidate_CH_DMRs <- readRDS(
  file.path("..", "..", "nonCG", "objects", "list_of_candidate_CH_DMRs.rds"))
ns <- setNames(names(list_of_candidate_CH_DMRs), names(list_of_candidate_CH_DMRs))
list_of_CH_DMRs <- lapply(ns, function(n) {
  x <- list_of_candidate_CH_DMRs[[n]]
  x <- x[(x$fwer / x$successful_permutations) <= 0.05]
  if (grepl("\\+", n)) {
    strand(x) <- "+"
  } else {
    strand(x) <- "-"
  }
  x
})

list_of_CH_DMRs <- list(
  "CA-DMRs" = unname(c(list_of_CH_DMRs[["mCA (+)"]], list_of_CH_DMRs[["mCA (-)"]])),
  "CT-DMRs" = unname(c(list_of_CH_DMRs[["mCT (+)"]], list_of_CH_DMRs[["mCT (-)"]])))

lapply(names(list_of_CH_DMRs), function(n) {
  x <- list_of_CH_DMRs[[n]]
  seqlevels(x) <- seqlevels(intersect(si_bed, seqinfo(x)))
  seqinfo(x) <- merge(seqinfo(x), si_bed)
  # NOTE: Have to ignore strand to sort BED file in order required by
  #       bedToBigBed
  export(sort(x, ignore.strand = TRUE),
         file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                   paste0(n, ".bed")))
  # NOTE: Have to
  system(paste0("bedToBigBed ",
                file.path(extdir, "flow-sorted-brain-wgbs", "data", "BED",
                          paste0(n, ".bed")),
                " ",
                file.path(extdir, "hg19.chrom.sizes"),
                " ",
                file.path(extdir, "bigBED", paste0(n, ".bb"))))
})

# ==============================================================================
# Open chromatin regions (OCRs)
#

# NOTE: It would be nice to use the narrowPeak files to display the OCRs.
#       However, there is no way to specify the colour of narrowPeak or
#       bigNarrowPeak tracks on the UCSC genome browser and IGV just shows them
#       as ordinary BED files.

writeLines('table hg19
"Open chromatin regions (OCRs)"
(\n
  string  chrom;      "Reference sequence chromosome or scaffold"
  uint    chromStart; "Start position of feature on chromosome"
  uint    chromEnd;   "End position of feature on chromosome"
  string  name;       "Name of gene"
  uint    score;      "Score"
  char[1] strand;     "+ or - for strand"
  uint    thickStart; "Coding region start"
  uint    thickEnd;   "Coding region end"
  uint    reserved;   "Brain region colour"
)',
           "OCRs.as")
seqlevels(ocrs_BA9_neg) <- seqlevels(intersect(si_bed, seqinfo(ocrs_BA9_neg)))
seqinfo(ocrs_BA9_neg) <- merge(seqinfo(ocrs_BA9_neg), si_bed)
ocrs_BA9_neg$itemRgb <- unique(ATAC_CD$TISSUE_COLOR[ATAC_CD$TISSUE == "BA9"])
export(sort(ocrs_BA9_neg),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "OCRs.BA9_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "OCRs.BA9_neg.bed"),
              " -as=OCRs.as",
              " -type=bed9",
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "OCRs.BA9_neg.bb")))

seqlevels(ocrs_BA9_pos) <- seqlevels(intersect(si_bed, seqinfo(ocrs_BA9_pos)))
seqinfo(ocrs_BA9_pos) <- merge(seqinfo(ocrs_BA9_pos), si_bed)
ocrs_BA9_pos$itemRgb <- unique(ATAC_CD$TISSUE_COLOR[ATAC_CD$TISSUE == "BA9"])
export(sort(ocrs_BA9_pos),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "OCRs.BA9_pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "OCRs.BA9_pos.bed"),
              " -as=OCRs.as",
              " -type=bed9",
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "OCRs.BA9_pos.bb")))

seqlevels(ocrs_NAcc_neg) <- seqlevels(intersect(si_bed, seqinfo(ocrs_NAcc_neg)))
seqinfo(ocrs_NAcc_neg) <- merge(seqinfo(ocrs_NAcc_neg), si_bed)
ocrs_NAcc_neg$itemRgb <- unique(ATAC_CD$TISSUE_COLOR[ATAC_CD$TISSUE == "NA"])
export(sort(ocrs_NAcc_neg),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "OCRs.NAcc_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "OCRs.NAcc_neg.bed"),
              " -as=OCRs.as",
              " -type=bed9",
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "OCRs.NAcc_neg.bb")))

seqlevels(ocrs_NAcc_pos) <- seqlevels(intersect(si_bed, seqinfo(ocrs_NAcc_pos)))
seqinfo(ocrs_NAcc_pos) <- merge(seqinfo(ocrs_NAcc_pos), si_bed)
ocrs_NAcc_pos$itemRgb <- unique(ATAC_CD$TISSUE_COLOR[ATAC_CD$TISSUE == "NA"])
export(sort(ocrs_NAcc_pos),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "OCRs.NAcc_pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "OCRs.NAcc_pos.bed"),
              " -as=OCRs.as",
              " -type=bed9",
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "OCRs.NAcc_pos.bb")))

seqlevels(ocrs_overall) <- seqlevels(intersect(si_bed, seqinfo(ocrs_overall)))
seqinfo(ocrs_overall) <- merge(seqinfo(ocrs_overall), si_bed)
ocrs_overall$itemRgb <- "black"
export(sort(ocrs_overall),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "OCRs.union.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "OCRs.union.bed"),
              " -as=OCRs.as",
              " -type=bed9",
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "OCRs.union.bb")))

# ==============================================================================
# DARs
#

dars_pos_vs_neg <- GRanges(read_csv(
  file.path("..", "..", "ATAC-seq", "extdata",
            "DARs.ave_pos_vs_ave_neg.ATAC-seq.csv.gz")))
seqlevels(dars_pos_vs_neg) <- seqlevels(intersect(si_bed,
                                                  seqinfo(dars_pos_vs_neg)))
seqinfo(dars_pos_vs_neg) <- merge(seqinfo(dars_pos_vs_neg), si_bed)
export(sort(dars_pos_vs_neg),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "DARs.pos_vs_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "DARs.pos_vs_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED", "DARs.pos_vs_neg.bb")))

dars_pos <- GRanges(read_csv(
  file.path("..", "..", "ATAC-seq", "extdata",
            "DARs.NA_posvsBA9_pos.ATAC-seq.csv.gz")))
seqlevels(dars_pos) <- seqlevels(intersect(si_bed,
                                           seqinfo(dars_pos)))
seqinfo(dars_pos) <- merge(seqinfo(dars_pos), si_bed)
export(sort(dars_pos),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "DARs.NAcc_pos_vs_BA9_pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "DARs.NAcc_pos_vs_BA9_pos.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED",
                        "DARs.NAcc_pos_vs_BA9_pos.bb")))

dars_neg <- GRanges(read_csv(
  file.path("..", "..", "ATAC-seq", "extdata",
            "DARs.NA_negvsBA9_neg.ATAC-seq.csv.gz")))
seqlevels(dars_neg) <- seqlevels(intersect(si_bed, seqinfo(dars_neg)))
seqinfo(dars_neg) <- merge(seqinfo(dars_neg), si_bed)
export(sort(dars_neg),
       file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                 "DARs.NAcc_neg_vs_BA9_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-atac", "data", "BED",
                        "DARs.NAcc_neg_vs_BA9_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED",
                        "DARs.NAcc_neg_vs_BA9_neg.bb")))

# ==============================================================================
# RNA-seq
#

# NOTE: Just a BED file of start/end of GENCODE gene model
# NOTE: No easy way to encode logFC in BED format

seqlevels(degs) <- seqlevels(intersect(si_bed, seqinfo(degs)))
seqinfo(degs) <- merge(seqinfo(degs), si_bed)
names(degs) <- degs$gene_id
export(
  sort(degs, ignore.strand = TRUE),
  file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
            "DEGs.NAcc_pos_vs_BA9_pos.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
                        "DEGs.NAcc_pos_vs_BA9_pos.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED",
                        "DEGs.NAcc_pos_vs_BA9_pos.bb")))

degs_pos_vs_neg_df <- read_csv(
  file.path("..", "..", "RNA-seq", "extdata",
            "DEGs.ave_pos_vs_ave_neg.RNA-seq.csv.gz"))
degs_pos_vs_neg <- unflattened_features$genes[
  degs_pos_vs_neg_df[degs_pos_vs_neg_df$adj.P.Val < 0.05, ][["X1"]]]
seqlevels(degs_pos_vs_neg) <- seqlevels(intersect(si_bed,
                                                  seqinfo(degs_pos_vs_neg)))
seqinfo(degs_pos_vs_neg) <- merge(seqinfo(degs_pos_vs_neg), si_bed)
export(
  sort(degs_pos_vs_neg, ignore.strand = TRUE),
  file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
            "DEGs.pos_vs_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
                        "DEGs.pos_vs_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED",
                        "DEGs.pos_vs_neg.bb")))

degs_neg_df <- read_csv(
  file.path("..", "..", "RNA-seq", "extdata",
            "DEGs.NA_negvsBA9_neg.RNA-seq.csv.gz"))
degs_neg <- unflattened_features$genes[
  degs_neg_df[degs_neg_df$adj.P.Val < 0.05, ][["X1"]]]
seqlevels(degs_neg) <- seqlevels(intersect(si_bed, seqinfo(degs_neg)))
seqinfo(degs_neg) <- merge(seqinfo(degs_neg), si_bed)
export(
  sort(degs_neg, ignore.strand = TRUE),
  file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
            "DEGs.NAcc_neg_vs_BA9_neg.bed"))
system(paste0("bedToBigBed ",
              file.path(extdir, "flow-sorted-brain-rna-seq", "data", "BED",
                        "DEGs.NAcc_neg_vs_BA9_neg.bed"),
              " ",
              file.path(extdir, "hg19.chrom.sizes"),
              " ",
              file.path(extdir, "bigBED",
                        "DEGs.NAcc_neg_vs_BA9_neg.bb")))

# ==============================================================================
# Brain enhancers?
#

# TODO: Are we providing these ourselves?
