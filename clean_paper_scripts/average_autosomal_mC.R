# Autosomal average (raw) mC (mCA, mCC, mCG, mCT)
# Peter Hickey
# 2018-03-25

# NOTE: Requires an absurd amount of memory to process mCH in this way
#       (qrsh -l mem_free=480G,h_vmem=481G)

library(bsseq)
library(HDF5Array)
library(matrixStats)
library(BSgenome.Hsapiens.UCSC.hg19)

options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

### ============================================================================
### Construct GRanges with locations of CA, CC, and CT dinucleotides on each
### strand of hg19
###

Cs <- c("CA", "CC", "CT")
names(Cs) <- Cs
list_of_Cs_gr <- lapply(Cs, function(C) {
  si <- keepSeqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg19),
                      paste0("chr", 1:22))
  params <- new("BSParams",
                X = BSgenome.Hsapiens.UCSC.hg19,
                FUN = matchPattern,
                exclude = setdiff(seqnames(BSgenome.Hsapiens.UCSC.hg19),
                                  seqnames(si)))
  options(verbose = TRUE)
  C_pos <- IRangesList(endoapply(bsapply(params, pattern = C), as, "IRanges"))
  C_pos <- GRanges(seqnames = Rle(names(C_pos), lengths(C_pos)),
                  ranges = unlist(C_pos, use.names = FALSE),
                  strand = "+",
                  seqinfo = si)
  width(C_pos) <- 1L

  C_neg <- IRangesList(
    endoapply(bsapply(params, pattern = reverseComplement(DNAString(C))),
              as, "IRanges"))
  C_neg <- GRanges(seqnames = Rle(names(C_neg), lengths(C_neg)),
                   ranges = unlist(C_neg, use.names = FALSE),
                   strand = "-",
                   seqinfo = si)
  C_neg <- resize(C_neg, width = 1, fix = "start")
  options(verbose = FALSE)
  list("+" = C_pos,
       "-" = C_neg)
})

### ============================================================================
### Load data and compute average mC in each Tissue:NeuN combination
###

# ------------------------------------------------------------------------------
# mCH
#

# NOTE: This object was saved with an old version of
#       SummarizedExperiment/HDF5Array so can't simply use
#       loadHDF5SummarizedExperiment()
CH_BSseq <- readRDS(
  "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/CH-flow-sorted-brain-wgbs/se.rds")
assays(CH_BSseq) <- list(
  M =
    HDF5Array("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/CH-flow-sorted-brain-wgbs/assays.h5",
              name = "assay001"),
  Cov =
    HDF5Array("/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/CH-flow-sorted-brain-wgbs/assays.h5",
              name = "assay002"))
CH_BSseq <- keepSeqlevels(CH_BSseq, paste0("chr", 1:22),
                          pruning.mode = "coarse")

# NOTE: Data are too large to load into memory at once. So split it up by
#       NeuN+ and NeuN- samples
CH_BSseq_pos <- CH_BSseq[, CH_BSseq$NeuN == "pos"]
CH_BSseq_neg <- CH_BSseq[, CH_BSseq$NeuN == "neg"]
list_of_CH_BSseq <- list("NeuN+" = CH_BSseq_pos,
                         "NeuN-" = CH_BSseq_neg)
list_of_ave_mCH_df <- lapply(list_of_CH_BSseq, function(BSseq) {
  mCH <- SummarizedExperiment(
    assays = as.matrix(getMeth(BSseq, type = "raw")),
    rowRanges = rowRanges(BSseq),
    colData = colData(BSseq))
  tissues <- unique(BSseq$Tissue)
  names(tissues) <- tissues
  neun <- unique(BSseq$NeuN)
  names(neun) <- neun
  ave_mCH <- lapply(tissues, function(tissue) {
    lapply(neun, function(neun) {
      cols <- which(mCH$Tissue == tissue & mCH$NeuN == neun)
      Cs <- c("CA", "CC", "CT")
      names(Cs) <- Cs
      list_of_cm <- lapply(Cs, function(C) {
        Cs_gr <- c(list_of_Cs_gr[[C]][["+"]], list_of_Cs_gr[[C]][["-"]])
        ol <- findOverlaps(BSseq, Cs_gr, type = "equal")
        rows <- queryHits(ol)
        cm <- colMeans2(assay(mCH, withDimnames = FALSE),
                        rows = rows,
                        cols = cols,
                        na.rm = TRUE)
        names(cm) <- colnames(mCH)[cols]
        cm
      })
    })
  })
  context <- sapply(strsplit(names(unlist(ave_mCH)), "\\."), "[[", 3)
  sample_names <- sapply(strsplit(names(unlist(ave_mCH)), "\\."), "[[", 4)
  ave_mCH_df <- data.frame(Context = context,
                           Donor = BSseq$Individual[match(sample_names,
                                                          colnames(BSseq))],
                           Tissue = BSseq$Tissue[match(sample_names,
                                                       colnames(BSseq))],
                           NeuN = BSseq$NeuN[match(sample_names,
                                                   colnames(BSseq))],
                           ave_mC = unlist(ave_mCH, use.names = FALSE),
                           stringsAsFactors = FALSE)
  ave_mCH_df
})
ave_mCH_df <- do.call(rbind, list_of_ave_mCH_df)
row.names(ave_mCH_df) <- NULL

# ------------------------------------------------------------------------------
# mCG
#

CG_BSseq <- loadHDF5SummarizedExperiment(
  dir = "/dcl01/hansen/data/flow-sorted-brain-wgbs/objects/BS.fit.small.sorted.somatic.all")
CG_BSseq <- keepSeqlevels(CG_BSseq, paste0("chr", 1:22))

mCG <- SummarizedExperiment(
  assays = as.matrix(getMeth(CG_BSseq, type = "raw")),
  rowRanges = rowRanges(CG_BSseq),
  colData = colData(CG_BSseq))

tissues <- unique(mCG$Tissue)
names(tissues) <- tissues
neun <- unique(mCG$NeuN)
names(neun) <- neun
ave_mCG <- lapply(tissues, function(tissue) {
  lapply(neun, function(neun) {
    cols <- which(mCG$Tissue == tissue & mCG$NeuN == neun)
    cm <- colMeans2(assay(mCG, withDimnames = FALSE),
                    cols = cols,
                    na.rm = TRUE)
    names(cm) <- colnames(mCG)[cols]
    cm
  })
})
sample_names <- sapply(strsplit(names(unlist(ave_mCG)), "\\."), "[[", 3)
ave_mCG_df <- data.frame(Context = "CG",
                         Donor = mCG$Individual[match(sample_names,
                                                      colnames(mCG))],
                         Tissue = mCG$Tissue[match(sample_names,
                                                   colnames(mCG))],
                         NeuN = mCG$NeuN[match(sample_names,
                                               colnames(mCG))],
                         ave_mC = unlist(ave_mCG, use.names = FALSE),
                         stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Combine and save
#

ave_mC_df <- rbind(ave_mCG_df, ave_mCH_df)
saveRDS(ave_mC_df, "average_autosomal_mC.rds")

### ============================================================================
### Plots
###

df <- readRDS("average_autosomal_mC.rds")

g <- ggplot(aes(x = Tissue, y = 100 * ave_mC, colour = Tissue),
            data = df) +
  facet_grid(Context ~ NeuN, scales = "free_y") +
  geom_jitter(width = 0.1) +
  scale_color_manual(
    values = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"))
save_plot("Global_C_Meth_autosomes.pdf", g, base_height = 5)

g1 <- ggplot(aes(x = Tissue, y = 100 * ave_mC, colour = Tissue),
             data = filter(df, Context == "CG")) +
  facet_grid(Context ~ NeuN) +
  geom_jitter(width = 0.1) +
  scale_color_manual(
    values = c("deeppink", "deepskyblue", "darkgrey", "chocolate1")) +
  ylim(c(75, 90))
save_plot("Global_CG_Meth_autosomes.pdf", g1, base_height = 5)

g2 <- ggplot(aes(x = Tissue, y = 100 * ave_mC, colour = Tissue),
             data = filter(df, Context != "CG")) +
  facet_grid(Context ~ NeuN) +
  geom_jitter(width = 0.1) +
  ylim(c(0, 13)) +
  scale_color_manual(
    values = c("deeppink", "deepskyblue", "darkgrey", "chocolate1"))
save_plot("Global_CH_Meth_autosomes.pdf", g2, base_height = 5)


df <- df %>%
  mutate(
    Tissue = factor(Tissue, levels = c("BA9", "BA24", "HC", "NA")),
    Context = factor(Context, levels = c("CG", "CA", "CC", "CT")))

df %>%
  group_by(Context, NeuN, Tissue) %>%
  summarise(mean = mean(ave_mC),
            sd = sd(ave_mC))

# Significant: NA-BA9
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CA", NeuN == "pos")))
# Significant: NA-BA9
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CC", NeuN == "pos")))
# Significant: NA-BA9, HC-BA24, NA-HC
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CG", NeuN == "pos")))
# Significant: NA-BA9, NA-BA24
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CT", NeuN == "pos")))

# Significant: NA-BA9
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CA", NeuN == "neg")))
# Significant: None
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CC", NeuN == "neg")))
# Significant: None
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CG", NeuN == "neg")))
# Significant: None
TukeyHSD(
  aov(
    ave_mC ~ Donor + Tissue,
    data = filter(df, Context == "CT", NeuN == "neg")))

lm_fit <- lm(
  ave_mC ~ Donor + Tissue * Context * NeuN,
  data = df)
anova(lm_fit)
aov_fit <- aov(
  ave_mC ~ Donor + Tissue * Context * NeuN,
  data = df)
TukeyHSD(aov_fit)

lm_fit <- lm(
  ave_mC ~ Donor + Tissue + Context * NeuN,
  data = df)
aov_fit <- aov(
  ave_mC ~ Donor + Tissue + Context * NeuN,
  data = df)

df %>%
  group_by(Tissue, Context, NeuN)

# Remember, these two are equivalent
t.test(
  ave_mC ~ Tissue,
  data = filter(
    df,
    Tissue %in% c("BA9", "NA"),
    NeuN == "pos",
    Context == "CG"),
  paired = TRUE)
lm_fit <- lm(
  ave_mC ~ Donor + Tissue,
  data = filter(df, Context == "CG", Tissue %in% c("BA9", "NA"), NeuN == "pos"))

lm_fit <- lm(
  ave_mC ~ Donor + Tissue * NeuN,
  data = filter(df, Context == "CG"))
aov_fit <- aov(
  ave_mC ~ Error(Donor) + Tissue * NeuN,
  data = filter(df, Context == "CG"))
lmer_fit <- lmer(
  ave_mC ~ (1 | Donor) + Tissue * NeuN,
  data = filter(df, Context == "CG"))
