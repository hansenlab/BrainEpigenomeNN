# Correlation of DiffEpi with gene expression around genes
# Peter Hickey
# 2018-02-07

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(cowplot)

load("../../integrating-dmrs-dars-and-degs/objects/assays-and-features.rda")
big_dars_pos <- dars_pos[abs(dars_pos$logFC) > 1]
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_candidate_CH_DMRs <- lapply(
  X = list_of_candidate_CH_DMRs,
  FUN = function(x) x[x$fwer <= 50])
list_of_candidate_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x$meanDiff <- x$NA_pos - x$BA9_pos
  x
})
strand(list_of_candidate_CH_DMRs[["mCA (+)"]]) <- "+"
strand(list_of_candidate_CH_DMRs[["mCA (-)"]]) <- "-"
strand(list_of_candidate_CH_DMRs[["mCT (+)"]]) <- "+"
strand(list_of_candidate_CH_DMRs[["mCT (-)"]]) <- "-"

f <- function(genes, feature, yname, ss = 0.01, ngb = 2) {
  stopifnot(ss > 0, ss < 1)
  w <- width(genes)
  # Resize the genes to include number of gene bodies (ngb) upstream and
  # downstream of the gene
  x <- resize(resize(genes, w + ngb * w, fix = "start"),
              w + 2 * ngb * w,
              fix = "end")
  upstream <- 0
  step_width <- ss * width(x)
  p <- seq(ss, 1, ss)
  cor_tbl <- lapply(p, function(pp) {
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_TSS <- start(x)[on_plus]
    start(x)[on_plus] <- on_plus_TSS - upstream
    end(x)[on_plus] <- on_plus_TSS + pp * width(x)[on_plus] - 1L
    on_minus <- which(strand(x) == "-")
    on_minus_TSS <- end(x)[on_minus]
    end(x)[on_minus] <- on_minus_TSS + upstream
    start(x)[on_minus] <- on_minus_TSS - pp * width(x)[on_minus] + 1L
    x <- resize(x, step_width, fix = "end")
    ol <- findOverlaps(feature, x)
    tmp <- data.frame(x = mcols(feature)[[yname]][queryHits(ol)],
                      y = x$logFC[subjectHits(ol)],
                      gene_id = names(x[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp), ])
    cor <- cor.test(tmp$x, tmp$y)
    data_frame(estimate = cor$estimate,
               lower = cor$conf.int[[1]],
               upper = cor$conf.int[[2]])
  })
  bind_cols(
    data_frame(x = p),
    bind_rows(cor_tbl))
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
       val$estimate,
       xlab = "% along gene body",
       ylab = "Correlation",
       main = n,
       type = "s",
       ylim = c(-1, 1))
  abline(v = 0, lty = 2)
  abline(v = 100, lty = 2)
}

genes <- gencode_features$genes[rna_atac_meth$gene]
genes$logFC <- inner_join(data_frame(gene = names(genes)),
                          dplyr::select(rna_atac_meth, gene, expLogFC)) %>%
  pull(expLogFC)
list_of_genes <- split(
  genes,
  ifelse(names(genes) %in% unlist(cgi_promoters$GENEID),
         "with_CGI_promoter",
         "without_CGI_promoter"))

### ============================================================================
### CG-DMRs
###

cg_dmr_val <- lapply(list_of_genes,
                     f,
                     feature = dmrs_NAvsBA9pos,
                     yname = "meanDiff",
                     ss = 0.01 / 5,
                     ngb = 2)
cg_dmr_df <- bind_rows(lapply(names(cg_dmr_val), function(n) {
  mutate(cg_dmr_val[[n]], state = n)
}))
par(mfrow = c(1, 2))
f2(cg_dmr_val[[1]], "with CGI", 2)
f2(cg_dmr_val[[2]], "without CGI", 2)
ggplot(cg_dmr_df, aes(x = x, y = estimate, col = state)) +
  geom_line()

### ============================================================================
### CA-DMRs
###

# ------------------------------------------------------------------------------
# pos strand
#

original_strand_ca_pos_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = list_of_candidate_CH_DMRs[["mCA (+)"]],
  yname = "meanDiff",
  ss = 0.01 / 5,
  ngb = 2)
original_strand_ca_pos_dmr_df <- bind_rows(
  lapply(names(original_strand_ca_pos_dmr_val), function(n) {
    mutate(original_strand_ca_pos_dmr_val[[n]], state = n)
  }))

unstranded_ca_pos_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = unstrand(list_of_candidate_CH_DMRs[["mCA (+)"]]),
                     yname = "meanDiff",
                     ss = 0.01 / 5,
                     ngb = 2)
unstranded_ca_pos_dmr_df <- bind_rows(
  lapply(names(unstranded_ca_pos_dmr_val), function(n) {
    mutate(unstranded_ca_pos_dmr_val[[n]], state = n)
  }))
inverted_strand_ca_pos_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = invertStrand(list_of_candidate_CH_DMRs[["mCA (+)"]]),
  yname = "meanDiff",
  ss = 0.01 / 5,
  ngb = 2)
inverted_strand_ca_pos_dmr_df <- bind_rows(
  lapply(names(inverted_strand_ca_pos_dmr_val), function(n) {
    mutate(inverted_strand_ca_pos_dmr_val[[n]], state = n)
  }))
par(mfrow = c(3 ,2))
f2(ca_pos_dmr_val[["with_CGI_promoter"]], "with CGI", 2)
f2(ca_pos_dmr_val[["without_CGI_promoter"]], "without CGI", 2)
f2(unstranded_ca_pos_dmr_val[["with_CGI_promoter"]], "US with CGI", 2)
f2(unstranded_ca_pos_dmr_val[["without_CGI_promoter"]], "US without CGI", 2)
f2(inverted_strand_ca_pos_dmr_val[["with_CGI_promoter"]], "IS with CGI", 2)
f2(inverted_strand_ca_pos_dmr_val[["without_CGI_promoter"]], "IS without CGI", 2)

ca_pos_dmr_list <- list(original = original_strand_ca_pos_dmr_df,
                      unstranded = unstranded_ca_pos_dmr_df,
                      inverted = inverted_strand_ca_pos_dmr_df)
ca_pos_dmr_df <- bind_rows(lapply(names(ca_pos_dmr_list), function(n) {
  mutate(ca_pos_dmr_list[[n]],
         strand = n)
}))
ggplot(ca_pos_dmr_df, aes(x = x, y = estimate, col = strand)) +
  geom_line() +
    facet_grid(~ state)

# ------------------------------------------------------------------------------
# neg strand
#

original_strand_ca_neg_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = list_of_candidate_CH_DMRs[["mCA (-)"]],
  yname = "meanDiff",
  ss = 0.01 / 5,
  ngb = 2)
original_strand_ca_neg_dmr_df <- bind_rows(
  lapply(names(original_strand_ca_neg_dmr_val), function(n) {
    mutate(original_strand_ca_neg_dmr_val[[n]], state = n)
  }))

unstranded_ca_neg_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = unstrand(list_of_candidate_CH_DMRs[["mCA (-)"]]),
  yname = "meanDiff",
  ss = 0.01 / 5,
  ngb = 2)
unstranded_ca_neg_dmr_df <- bind_rows(
  lapply(names(unstranded_ca_neg_dmr_val), function(n) {
    mutate(unstranded_ca_neg_dmr_val[[n]], state = n)
  }))
inverted_strand_ca_neg_dmr_val <- lapply(
  list_of_genes,
  f,
  feature = invertStrand(list_of_candidate_CH_DMRs[["mCA (-)"]]),
  yname = "meanDiff",
  ss = 0.01 / 5,
  ngb = 2)
inverted_strand_ca_neg_dmr_df <- bind_rows(
  lapply(names(inverted_strand_ca_neg_dmr_val), function(n) {
    mutate(inverted_strand_ca_neg_dmr_val[[n]], state = n)
  }))

par(mfrow = c(3 ,2))
f2(ca_neg_dmr_val[["with_CGI_promoter"]], "with CGI", 2)
f2(ca_neg_dmr_val[["without_CGI_promoter"]], "without CGI", 2)
f2(unstranded_ca_neg_dmr_val[["with_CGI_promoter"]], "US with CGI", 2)
f2(unstranded_ca_neg_dmr_val[["without_CGI_promoter"]], "US without CGI", 2)
f2(inverted_strand_ca_neg_dmr_val[["with_CGI_promoter"]], "IS with CGI", 2)
f2(inverted_strand_ca_neg_dmr_val[["without_CGI_promoter"]], "IS without CGI", 2)

ca_neg_dmr_list <- list(original = original_strand_ca_neg_dmr_df,
                        unstranded = unstranded_ca_neg_dmr_df,
                        inverted = inverted_strand_ca_neg_dmr_df)
ca_neg_dmr_df <- bind_rows(lapply(names(ca_neg_dmr_list), function(n) {
  mutate(ca_neg_dmr_list[[n]],
         strand = n)
}))
ggplot(ca_neg_dmr_df, aes(x = x, y = estimate, col = strand)) +
  geom_line() +
  facet_grid(~ state)

### ============================================================================
### DARs
###

# UP TO HERE: ggplot versions

dars_pos_val <- lapply(
  list_of_genes,
  f,
  feature = dars_pos,
  yname = "logFC",
  ss = 0.01 / 5,
  ngb = 2)
dars_pos_df <- bind_rows(lapply(names(dars_pos_val), function(n) {
  mutate(dars_pos_val[[n]],
         state = n)
}))
par(mfrow = c(1 ,2))
f2(dars_pos_val[["with_CGI_promoter"]], "with CGI", 2)
f2(dars_pos_val[["without_CGI_promoter"]], "without CGI", 2)
ggplot(dars_pos_df, aes(x = x, y = estimate, col = state)) +
  geom_line()

### ============================================================================
### DARs
###

big_dars_pos_val <- lapply(
  list_of_genes,
  f,
  feature = big_dars_pos,
  yname = "logFC",
  ss = 0.01 / 5,
  ngb = 2)
big_dars_pos_df <- bind_rows(lapply(names(big_dars_pos_val), function(n) {
  mutate(big_dars_pos_val[[n]],
         state = n)
}))
par(mfrow = c(1 ,2))
f2(big_dars_pos_dmr_val[["with_CGI_promoter"]], "with CGI", 2)
f2(big_dars_pos_dmr_val[["without_CGI_promoter"]], "without CGI", 2)
ggplot(big_dars_pos_df, aes(x = x, y = estimate, col = state)) +
  geom_line()
