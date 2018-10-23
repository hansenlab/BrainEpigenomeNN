# Count proportion of dinucleotides around gene bodies
# Peter Hickey
# 2018-02-11

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(matrixStats)

options("mc.cores" = 16)
extdir <- "../extdata"

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")

# Compute proportions
computeProportions <- function(genes, feature, ss = 0.01, ngb = 2) {
  stopifnot(ss > 0, ss < 1)
  w <- width(genes)
  # Resize the genes to include number of gene bodies (ngb) upstream and
  # downstream of the gene
  genes <- resize(resize(genes, w + ngb * w, fix = "start"),
                  w + 2 * ngb * w,
                  fix = "end")
  w <- ss * width(genes)
  p <- seq(ss, 1, ss)
  ol_tbl <- sapply(p, function(pp) {
    on_plus <- which(strand(genes) == "+" | strand(genes) == "*")
    on_plus_TSS <- start(genes)[on_plus]
    end(genes)[on_plus] <- on_plus_TSS + pp * width(genes)[on_plus] - 1L
    on_minus <- which(strand(genes) == "-")
    on_minus_TSS <- end(genes)[on_minus]
    start(genes)[on_minus] <- on_minus_TSS - pp * width(genes)[on_minus] + 1L
    genes <- resize(genes, w, fix = "end")
    countOverlaps(genes, feature) / width(genes)
  })
  perc <- 100 * colMeans2(ol_tbl, na.rm = TRUE)
  return(list(x = p, y = perc))
}

plotProportions <- function(val, n, ngb) {
  x <- 100 * val$x
  inc <- 2 * ngb + 1
  first_cut <- length(x) * ngb / inc
  second_cut <- length(x) * (ngb + 1) / inc
  x[seq(1, first_cut)] <-
    -rev(seq(0, ngb * 100 - 1, length.out = first_cut))
  x[seq(first_cut, second_cut)] <-
    seq(0, 100, length.out = second_cut - first_cut + 1)
  x[seq(second_cut, length(x))] <-
    seq(100, (ngb + 1) * 100, length.out = length(x) - second_cut + 1)
  plot(x,
       val$y,
       xlab = "% along gene body",
       ylab = "CpX density",
       main = n,
       type = "s",
       ylim = c(0, 10))
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
}

bsp <- new("BSParams",
           X = BSgenome.Hsapiens.UCSC.hg19,
           FUN = matchPattern,
           exclude = c("M", "_"))

dinucs <- levels(interaction(DNA_BASES, DNA_BASES, sep = ""))
names(dinucs) <- dinucs

# ------------------------------------------------------------------------------
# Stratify genes and transcripts be CGI-promoter status and PC/lncRNA status
#

genes <- gencode_features$genes
genes_idx <- interaction(
  ifelse(names(genes) %in% pc_genes, "PC", "lncRNA"),
  ifelse(names(genes) %in% unlist(cgi_promoters$GENEID),
         "with_CGI_promoter",
         "without_CGI_promoter"),
  sep = "_")
list_of_genes <- split(genes, genes_idx)

tx <- unlist(gencode_features$transcripts_by_gene)
tx$gene_name <- names(tx)
tx <- unname(tx)
tx_idx <- interaction(
  ifelse(tx$gene_name %in% pc_genes, "PC", "lncRNA"),
  ifelse(tx$tx_name %in% cgi_promoters$TXNAME,
         "with_CGI_promoter",
         "without_CGI_promtoer"),
  sep = "_")
list_of_tx <- split(tx, tx_idx)

# ------------------------------------------------------------------------------
# Construct GRanges of all dinucleotides
#

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

# ------------------------------------------------------------------------------
# Gene-level analysis
#

gene_props <- mclapply(dinucs, function(dinuc, gene_models, ngb = 2,
                                        ss = 0.01 / (1 + 2 * ngb)) {
  message("Loading ", dinuc)
  dinuc_gr <- readRDS(file.path(extdir, "flow-sorted-brain-wgbs", "objects",
                                paste0(dinuc, ".GRanges.rds")))
  message("Original strand ", dinuc)
  val <- lapply(gene_models,
                computeProportions,
                feature = dinuc_gr,
                ss = ss,
                ngb = ngb)
  message("Opposite strand ", dinuc)
  opposite_strand_val <- lapply(gene_models,
                                computeProportions,
                                feature = invertStrand(dinuc_gr),
                                ss = ss,
                                ngb = ngb)
  message("Saving ", dinuc)
  list(val = val,
       opposite_strand_val = opposite_strand_val)
}, gene_models = list_of_genes)

# TODO: Update plotProportions to handle all strands?
pdf("~/tmp.pdf", 7, 7)
par(mfrow = c(4, 4))
lapply(dinucs, function(dinuc) {
  plotProportions(gene_props[[dinuc]][["val"]][["PC_without_CGI_promoter"]],
                  dinuc, 2)
})
dev.off()
