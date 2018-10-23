# FWER procedure for CH-DMRs
# Peter Hickey
# 2017-09-06

library(parallel)
library(GenomicRanges)

options("mc.cores" = 40)

extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)

### ============================================================================
### Functions
###

# NOTE: A parallelised version of bsseq:::getFWER.fstat
getFWER.fstat <- function(null, type = "blocks",
                          mc.cores = getOption("mc.cores")) {
  reference <- null[[1]]
  null <- null[-1]
  null <- null[!sapply(null, is.null)]
  better <- unlist(mclapply(seq_len(nrow(reference)), function(ii) {
    message(ii)
    areaStat <- abs(reference$areaStat[ii])
    width <- reference$width[ii]
    n <- reference$n[ii]
    if (type == "blocks") {
      out <- vapply(null, function(nulldist) {
        any(abs(nulldist$areaStat) >= areaStat & nulldist$width >=
              width)
      }, logical(1L))
    }
    if (type == "dmrs") {
      out <- vapply(null, function(nulldist) {
        any(abs(nulldist$areaStat) >= areaStat & nulldist$n >=
              n)
      }, logical(1L))
    }
    sum(out)
  }, mc.cores = mc.cores))
  better
}

### ============================================================================
### Compute FWER and meanDiff
###

list_of_val <- mapply(function(strand, context) {
  message("context = ", context, ", strand = ", strand)
  message("Loading p")
  p <- readRDS(paste0("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/",
                      "null_DMRs/", strand, "_", context, "_small.",
                      "permutationMatrix.rds"))
  files <- list.files(
    "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/null_DMRs",
    pattern = paste0("^", strand, "_", context,
                     "_small\\.permutation.*dmrs\\.rds$"),
    full.names = TRUE)
  names(files) <- as.numeric(
    gsub("\\.dmrs\\.rds", "",
         gsub(paste0(strand, "\\_", context,
                     "\\_small\\.permutation\\_"), "",
              basename(files))))
  message("Loading nulls")
  null_dmrs <- mclapply(
    X = files,
    FUN = function(xx) {
      y <- readRDS(xx)
      if (is.null(y)) {
        return(NULL)
      }
      makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
    }, mc.cores = getOption("mc.cores"))
  message(paste0(length(null_dmrs), "/", nrow(p), " permutations succeeded"))

  message("Loading DMRs")
  dmrs <- readRDS(paste0("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/",
                         "null_DMRs/", strand, "_", context, "_small.",
                         "dmrs.rds"))

  message("Computing FWER")
  dmrs$fwer <- getFWER.fstat(c(list(dmrs), null_dmrs), type = "dmrs")
  dmrs$successful_permutations <- length(null_dmrs)

  message("Constructing and saving objects")
  null_dmrs <- null_dmrs[!S4Vectors:::sapply_isNULL(null_dmrs)]
  null_dmrs <- GRangesList(null_dmrs)
  dmrs <- makeGRangesFromDataFrame(dmrs, keep.extra.columns = TRUE)
  val <- list(dmrs = dmrs, null_dmrs = null_dmrs)
  saveRDS(val,
          paste0("../objects/", strand, "_", context,
                 ".DMRs_and_null_DMRs.rds"))
  val
}, strand = strands, context = contexts, SIMPLIFY = FALSE)
list_of_candidate_CH_DMRs <- lapply(list_of_val, "[[", "dmrs")

### ============================================================================
### Mean mCG and mCH in CH-DMRs
###

list_of_features <- list_of_candidate_CH_DMRs
list_of_BSseq <- c(list_of_CG_BSseq, list_of_CH_BSseq)
list_of_meanMeth <- lapply(list_of_features, function(features) {
  message("N_features = ", length(features))
  lapply(list_of_BSseq, function(BSseq) {
    fac <- interaction(BSseq$Tissue, BSseq$NeuN, sep = "_")
    BSseq <- subsetByOverlaps(BSseq, features)
    meth <- as.matrix(getMeth(BSseq))
    ov <- findOverlaps(BSseq, features)
    out <- lapply(split(meth, subjectHits(ov)), matrix, ncol = ncol(meth))
    outMatrix <- matrix(NA_real_, ncol = nlevels(fac), nrow = length(features),
                        dimnames = list(NULL, levels(fac)))
    mean_meth <- lapply(out, matrixStats::colMeans2, na.rm = TRUE)
    group_mean_meth <- do.call(
      rbind,
      lapply(mean_meth, function(xx) tapply(xx, fac, mean)))
    outMatrix[as.integer(rownames(group_mean_meth)), ] <- group_mean_meth
    outMatrix
  })
})

saveRDS(list_of_meanMeth, "../objects/list_of_meanMeth.rds")

# ------------------------------------------------------------------------------
# Make SummarizedExperiment objects with features as rows and meanMeth as assays
#

fn <- c("mCA (+)" = "pos_CA_DMRs",
        "mCA (-)" = "neg_CA_DMRs",
        "mCT (+)" = "pos_CT_DMRs",
        "mCT (-)" = "neg_CT_DMRs")
list_of_SE <- mapply(function(features, meanMeth, n) {
  cn <- c("idxStart", "idxEnd", "cluster", "n", "invdensity", "areaStat",
          "maxStat", "fwer", "successful_permutations")
  rr <- features
  mcols(rr) <- mcols(rr)[colnames(mcols(rr)) %in% cn]
  se <- SummarizedExperiment(assays = meanMeth, rowRanges = rr)
  saveRDS(se,
          paste0("../objects/", n, ".rds"))
  se
}, features = list_of_features, meanMeth = list_of_meanMeth,
n = fn[names(list_of_features)])
saveRDS(list_of_SE, "../objects/list_of_SE.CH-DMRs.rds")

# ------------------------------------------------------------------------------
# Add meanMeth in each tissue to candidate DMRs (only for the relevant strand
# and context)
#

ns <- names(list_of_candidate_CH_DMRs)
names(ns) <- ns
list_of_candidate_CH_DMRs <- lapply(ns, function(n) {
  DMRs <- list_of_candidate_CH_DMRs[[n]]
  meanMeth <- list_of_meanMeth[[n]][[n]]
  mcols(DMRs) <- cbind(mcols(DMRs), as(meanMeth, "DataFrame"))
  DMRs
})
saveRDS(list_of_candidate_CH_DMRs,
        "../objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
