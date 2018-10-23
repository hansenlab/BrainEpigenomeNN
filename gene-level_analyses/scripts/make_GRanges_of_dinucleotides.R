# Construct GRanges of all dinucleotides in hg19
# Peter Hickey
# 2018-02-16

### ----------------------------------------------------------------------------
### Setup
###

library(BSgenome.Hsapiens.UCSC.hg19)

extdir <- "../extdata"

bsp <- new("BSParams",
           X = BSgenome.Hsapiens.UCSC.hg19,
           FUN = matchPattern,
           exclude = c("M", "_"))

dinucs <- levels(interaction(DNA_BASES, DNA_BASES, sep = ""))
names(dinucs) <- dinucs

### ----------------------------------------------------------------------------
### Functions
###


dinucs_gr <- mclapply(dinucs, function(dinuc) {
  fwd <- bsapply(bsp, pattern = dinuc)
  fwd <- do.call(c, lapply(names(fwd), function(n) {
    GRanges(n, as(fwd[[n]], "IRanges"), strand = "+")
  }))
  fwd <- resize(fwd, 1, fix = "start")
  rc_dinuc <- as.character(reverseComplement(DNAString(dinuc)))
  rev <- bsapply(bsp, pattern = rc_dinuc)
  rev <- do.call(c, lapply(names(rev), function(n) {
    GRanges(n, as(rev[[n]], "IRanges"), strand = "-")
  }))
  rev <- resize(rev, 1, fix = "start")
  feature <- c(fwd, rev)
  saveRDS(feature, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                             paste0(dinuc, ".GRanges.rds")))
  feature
})

### ----------------------------------------------------------------------------
### Compute and save results
###

# NOTE: Results are saved within the function
list_of_dinucs_gr <- mclapply(dinucs, function(dinuc) {
  fwd <- bsapply(bsp, pattern = dinuc)
  fwd <- do.call(c, lapply(names(fwd), function(n) {
    GRanges(n, as(fwd[[n]], "IRanges"), strand = "+")
  }))
  fwd <- resize(fwd, 1, fix = "start")
  rc_dinuc <- as.character(reverseComplement(DNAString(dinuc)))
  rev <- bsapply(bsp, pattern = rc_dinuc)
  rev <- do.call(c, lapply(names(rev), function(n) {
    GRanges(n, as(rev[[n]], "IRanges"), strand = "-")
  }))
  rev <- resize(rev, 1, fix = "start")
  feature <- c(fwd, rev)
  saveRDS(feature, file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                             paste0(dinuc, ".GRanges.rds")))
  feature
})

