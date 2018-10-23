# Create M and Cov matrices to upload to GEO
# Lindsay Rizzardi and Peter Hickey
# 2018-03-22

library(data.table)
library(bsseq)
library(HDF5Array)

outdir <- "/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects"
setDTthreads(20)

options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

### ============================================================================
### mCG
###

# ------------------------------------------------------------------------------
# Sorted Samples
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/BS.fit.small.sorted.somatic.all")

regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))
M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.sorted.CpG_unstranded.txt"),
  sep = "\t")
system(paste0("gzip ", file.path(outdir, "M_matrix.sorted.CpG_unstranded.txt")))

Cov <- as.matrix(getCoverage(BSseq, type = "Cov"))
Coverage <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Coverage,
  file = file.path(outdir, "Cov_matrix.sorted.CpG_unstranded.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "Cov_matrix.sorted.CpG_unstranded.txt")))

# ------------------------------------------------------------------------------
# Unsorted samples
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/BS.unsorted.fit.small.somatic.all")
# NOTE: Drop the caudate samples
BSseq <- BSseq[, BSseq$Tissue != "caudate"]

regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))
M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.unsorted.CpG_unstranded.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "M_matrix.unsorted.CpG_unstranded.txt")))

Cov <- as.matrix(getCoverage(BSseq,type = "Cov"))
Coverage <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Coverage,
  file = file.path(outdir, "Cov_matrix.unsorted.CpG_unstranded.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "Cov_matrix.unsorted.CpG_unstranded.txt")))

### ============================================================================
### mCA
###

# ------------------------------------------------------------------------------
# mCA (+)
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/pos_CA_small-flow-sorted-brain-wgbs")
strand(BSseq) <- "+"
regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))

M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.sorted.CpA_forward_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "M_matrix.sorted.CpA_forward_strand.txt")))

Cov <- as.matrix(getCoverage(BSseq, type = "Cov"))
Cov <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Cov,
  file = file.path(outdir, "Cov_matrix.sorted.CpA_forward_strand.txt"),
  sep = "\t")
system(
  paste0(
    "gzip ", file.path(outdir, "Cov_matrix.sorted.CpA_forward_strand.txt")))

# ------------------------------------------------------------------------------
# mCA (-)
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/neg_CA_small-flow-sorted-brain-wgbs")
strand(BSseq) <- "-"
regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))
M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.sorted.CpA_reverse_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "M_matrix.sorted.CpA_reverse_strand.txt")))

Cov <- as.matrix(getCoverage(BSseq, type = "Cov"))
Cov <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Cov,
  file = file.path(outdir, "Cov_matrix.sorted.CpA_reverse_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "Cov_matrix.sorted.CpA_reverse_strand.txt")))

### ============================================================================
### mCT
###

# ------------------------------------------------------------------------------
# mCT (+)
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/pos_CT_small-flow-sorted-brain-wgbs")
strand(BSseq) <- "+"
regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))
M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.sorted.CpT_forward_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "M_matrix.sorted.CpT_forward_strand.txt")))

Cov <- as.matrix(getCoverage(BSseq, type = "Cov"))
Cov <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Cov,
  file = file.path(outdir, "Cov_matrix.sorted.CpT_forward_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "Cov_matrix.sorted.CpT_forward_strand.txt")))

# ------------------------------------------------------------------------------
# mCT (-)
#

BSseq <- loadHDF5SummarizedExperiment(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/neg_CT_small-flow-sorted-brain-wgbs")
strand(BSseq) <- "-"
regions <- data.table(
  chr = as.character(seqnames(BSseq)),
  start = start(BSseq),
  strand = as.factor(strand(BSseq)))
M <- as.matrix(getCoverage(BSseq, type = "M"))
Meth <- cbind(regions, as.data.table(M))
fwrite(
  x = Meth,
  file = file.path(outdir, "M_matrix.sorted.CpT_reverse_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "M_matrix.sorted.CpT_reverse_strand.txt")))

Cov <- as.matrix(getCoverage(BSseq, type = "Cov"))
Cov <- cbind(regions, as.data.table(Cov))
fwrite(
  x = Cov,
  file = file.path(outdir, "Cov_matrix.sorted.CpT_reverse_strand.txt"),
  sep = "\t")
system(
  paste0("gzip ", file.path(outdir, "Cov_matrix.sorted.CpT_reverse_strand.txt")))

# TODO: Delete old M and Cov txt files
