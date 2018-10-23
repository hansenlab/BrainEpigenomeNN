# Count proportion of CpXs around gene bodies
# Peter Hickey
# 2018-02-06

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(matrixStats)

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")

f <- function(genes, feature, ss = 0.01, ngb = 2) {
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

f2 <- function(val, n, ngb) {
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

list_of_genes <- list(
  with_CGI_promoter =
    gencode_features$genes[intersect(rna_atac_meth$gene,
                                     unlist(cgi_promoters$GENEID))],
  without_CGI_promoter =
    gencode_features$genes[setdiff(rna_atac_meth$gene,
                                   unlist(cgi_promoters$GENEID))])

### ============================================================================
### CpGs
###

# NOTE: Use all CpGs in the reference genome rather than just those tested for
#       differential methylation
cpgs <- bsapply(bsp, pattern = "CG")
cpgs <- do.call(c, lapply(names(cpgs), function(n) {
  GRanges(n, as(cpgs[[n]], "IRanges"), strand = "*")
}))

cpgs_val <- lapply(list_of_genes,
                   f,
                   feature = cpgs,
                   ss = 0.01 / 5,
                   ngb = 2)
saveRDS(cpgs_val, "../objects/cpgs_val.rds")

pdf("../figures/CG_around_genes.pdf", width = 7, height = 7)
par(mfrow = c(1, 2))
f2(cpgs_val$with_CGI_promoter,
   "CpG:\nwith CGI promoter",
   2)
f2(cpgs_val$without_CGI_promoter,
   "CpG:\nwithout CGI promoter",
   2)
dev.off()

### ============================================================================
### CpAs
###

# NOTE: Get CpAs on each strand, retaining strand information
cpas <- bsapply(bsp, pattern = "CA")
cpas <- do.call(c, lapply(names(cpas), function(n) {
  GRanges(n, as(cpas[[n]], "IRanges"), strand = "+")
}))
cpas <- resize(cpas, 1, fix = "start")
tpgs <- bsapply(bsp, pattern = "TG")
tpgs <- do.call(c, lapply(names(tpgs), function(n) {
  GRanges(n, as(tpgs[[n]], "IRanges"), strand = "-")
}))
tpgs <- resize(tpgs, 1, fix = "start")
cpas_tpgs <- c(cpas, tpgs)
saveRDS(cpas_tpgs, "../objects/cpas_tpgs.rds")

cpas_tpgs_val <- lapply(list_of_genes,
                   f,
                   feature = cpas_tpgs,
                   ss = 0.01 / 5,
                   ngb = 2)
saveRDS(cpas_tpgs_val, "../objects/cpas_tpgs_val.rds")
cpas_tpgs_opposite_strand_val <- lapply(list_of_genes,
                                        f,
                                        feature = invertStrand(cpas_tpgs),
                                        ss = 0.01 / 5,
                                        ngb = 2)
saveRDS(cpas_tpgs_opposite_strand_val,
        "../objects/cpas_tpgs_opposite_strand_val")

pdf("../figures/CA_around_genes.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))
f2(cpas_tpgs_val$with_CGI_promoter,
   "CpA/TpG:\nwith CGI promoter",
   2)
f2(cpas_tpgs_val$without_CGI_promoter,
   "CpA/TpG:\nwithout CGI promoter",
   2)
f2(cpas_tpgs_opposite_strand_val$with_CGI_promoter,
   "CpA/TpG (opp. strand):\nwith CGI promoter",
   2)
f2(cpas_tpgs_opposite_strand_val$without_CGI_promoter,
   "CpA/TpG (opp. strand):\nwithout CGI promoter",
   2)
dev.off()
