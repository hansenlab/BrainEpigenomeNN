# Compute and plot correlation of DiffEpi marks with DiffExp around TSSs of
# protein coding genes
# Peter Hickey
# 2018-02-12

library(GenomicRanges)
library(purrr)
library(data.table)
library(dplyr)
library(scales)

###---------------------------------------------------------------------------
### Load data
###

load("../../genomic-features/objects/unflattened-GENCODE-v19-features.rda")
load("../objects/assays-and-features.rda")
big_dars_pos <- dars_pos[abs(dars_pos$logFC) > 1]
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_candidate_CH_DMRs <- lapply(
  X = list_of_candidate_CH_DMRs,
  FUN = function(x) x[x$fwer <= 50])
CA_pos_DMRs <- list_of_candidate_CH_DMRs[["mCA (+)"]]

pc_transcripts <- unflattened_features_pc_transcripts$transcripts
# NOTE: names() used to identify transcripts from the same gene
names(pc_transcripts) <-
  unflattened_features_pc_transcripts$transcripts$GENEID

lnc_transcripts <- unflattened_features_lnc_transcripts$transcripts
# NOTE: names() used to identify transcripts from the same gene
names(lnc_transcripts) <-
  unflattened_features_lnc_transcripts$transcripts$GENEID

### ---------------------------------------------------------------------------
### Correlations
###

# NOTE: Not the most efficient way to do this, but I at least understand what
#       it's doing
f <- function(s, tx) {
  # NOTE: Don't compute s == 0 for downstream (computed in upstream to avoid
  #       zero-width ranges)

  # ----------------------------------------------------------------------------
  # CG-DMRs
  #
  mCG_downstream <- map_df(s[s != 0], function(ss) {
    p <- resize(promoters(tx,
                          downstream = ss,
                          upstream = 0),
                width = 1,
                fix = "end")
    ol <- findOverlaps(dmrs_NAvsBA9pos, p)
    tmp <- data.table(x = dmrs_NAvsBA9pos$meanDiff[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp)])
    # NOTE: Meth and expression inversely correlated, so take negative meth
    y <- cor.test(-tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = s[s != 0])

  mCG_upstream <- map_df(s, function(ss) {
    p <- resize(promoters(tx,
                          downstream = 0,
                          upstream = ss),
                width = 1,
                fix = "start")
    ol <- findOverlaps(dmrs_NAvsBA9pos, p)
    tmp <- data.table(x = dmrs_NAvsBA9pos$meanDiff[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp)])
    # NOTE: Meth and expression inversely correlated, so take negative meth
    y <- cor.test(-tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = -s)
  mCG <- bind_rows(mCG_downstream, mCG_upstream) %>%
    arrange(x)

  # ----------------------------------------------------------------------------
  # CH-DMRs
  #

  # NOTE: Using mCA (+) to represent mCH
  mCH_downstream <- map_df(s[s != 0], function(ss) {
    p <- resize(promoters(tx,
                          downstream = ss,
                          upstream = 0),
                width = 1,
                fix = "end")
    ol <- findOverlaps(CA_pos_DMRs, p)
    tmp <- data.table(x = CA_pos_DMRs$NA_pos[queryHits(ol)] -
                        CA_pos_DMRs$BA9_pos[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp)])
    # NOTE: Meth and expression inversely correlated, so take negative meth
    y <- cor.test(-tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = s[s != 0])

  mCH_upstream <- map_df(s, function(ss) {
    p <- resize(promoters(tx,
                          downstream = 0,
                          upstream = ss),
                width = 1,
                fix = "start")
    ol <- findOverlaps(CA_pos_DMRs, p)
    tmp <- data.table(x = CA_pos_DMRs$NA_pos[queryHits(ol)] -
                        CA_pos_DMRs$BA9_pos[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp)])
    # NOTE: Meth and expression inversely correlated, so take negative meth
    y <- cor.test(-tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = -s)
  mCH <- bind_rows(mCH_downstream, mCH_upstream) %>%
    arrange(x)

  # ----------------------------------------------------------------------------
  # DARs
  #

  atac_downstream <- map_df(s[s != 0], function(ss) {
    p <- resize(promoters(tx,
                          downstream = ss,
                          upstream = 0),
                width = 1,
                fix = "end")
    ol <- findOverlaps(dars_pos, p)
    tmp <- data.table(x = dars_pos$logFC[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp), ])
    y <- cor.test(tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = s[s != 0])

  atac_upstream <- map_df(s, function(ss) {
    p <- resize(promoters(tx,
                          downstream = 0,
                          upstream = ss),
                width = 1,
                fix = "start")
    ol <- findOverlaps(dars_pos, p)
    tmp <- data.table(x = dars_pos$logFC[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp), ])
    y <- cor.test(tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = y$conf.int[1],
               upper = y$conf.int[2])
  }) %>%
    mutate(x = -s)

  atac <- bind_rows(atac_downstream, atac_upstream) %>%
    arrange(x)

  # ----------------------------------------------------------------------------
  # bigDARs
  #

  big_atac_downstream <- map_df(s[s != 0], function(ss) {
    p <- resize(promoters(tx,
                          downstream = ss,
                          upstream = 0),
                width = 1,
                fix = "end")
    ol <- findOverlaps(big_dars_pos, p)
    tmp <- data.table(x = big_dars_pos$logFC[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp), ])
    y <- cor.test(tmp$x, tmp$y, use = "complete.obs")
    data_frame(cor = y$estimate,
               lower = ifelse(is.null(y$conf.int[1]), NA_real_,
                              is.null(y$conf.int[1])),
               upper = ifelse(is.null(y$conf.int[2]), NA_real_,
                              is.null(y$conf.int[2])))
  }) %>%
    mutate(x = s[s != 0])

  big_atac_upstream <- map_df(s, function(ss) {
    p <- resize(promoters(tx,
                          downstream = 0,
                          upstream = ss),
                width = 1,
                fix = "start")
    ol <- findOverlaps(big_dars_pos, p)
    tmp <- data.table(x = big_dars_pos$logFC[queryHits(ol)],
                      y = rna_seq_de_pos$gene_level$logFC[
                        match(names(p[subjectHits(ol)]),
                              rna_seq_de_pos$gene_level$gene_id)],
                      gene_id = names(p[subjectHits(ol)]))
    tmp <- unique(tmp[complete.cases(tmp), ])
    y <- try(cor.test(tmp$x, tmp$y, use = "complete.obs"),
             silent = TRUE)
    if (is(y, "try-error")) {
      return(data_frame(cor = NA_real_,
                        lower = NA_real_,
                        upper = NA_real_))
    }
    data_frame(cor = y$estimate,
               lower = ifelse(is.null(y$conf.int[1]), NA_real_,
                              is.null(y$conf.int[1])),
               upper = ifelse(is.null(y$conf.int[2]), NA_real_,
                              is.null(y$conf.int[2])))
  }) %>%
    mutate(x = -s)

  big_atac <- bind_rows(big_atac_downstream, big_atac_upstream) %>%
    arrange(x)

  # ----------------------------------------------------------------------------
  # Plot
  #

  # NOTE: Colours are from Dark2 palette
  #       (http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=4)
  cols <- c("mCG" = "#1b9e77",
            "mCH" = "#d95f02",
            "atac" = "#7570b3",
            "big_atac" = "#e7298a")

  plot(cor ~ x, data = mCG, type = "l", ylim = c(0, 1), col = cols["mCG"],
       lwd = 3,
       ylab = "|cor|",
       xlab = "Distance to TSS")
  polygon(x = c(mCG$x, rev(mCG$x)),
          y = c(mCG$lower, rev(mCG$upper)),
          col = alpha(cols["mCG"], 0.3),
          lty = 0)
  lines(cor ~ x, data = mCH, type = "l", ylim = c(0, 1), col = cols["mCH"],
        lwd = 3)
  polygon(x = c(mCG$x, rev(mCH$x)),
          y = c(mCH$lower, rev(mCH$upper)),
          col = alpha(cols["mCH"], 0.3),
          lty = 0)
  lines(cor ~ x, data = atac, type = "l", ylim = c(0, 1), col = cols["atac"],
        lwd = 3)
  polygon(x = c(atac$x, rev(atac$x)),
          y = c(atac$lower, rev(atac$upper)),
          col = alpha(cols["atac"], 0.3),
          lty = 0)
  lines(cor ~ x, data = big_atac, type = "l", ylim = c(0, 1),
        col = cols["big_atac"],
        lwd = 3)
  polygon(x = c(big_atac$x, rev(big_atac$x)),
          y = c(big_atac$lower, rev(big_atac$upper)),
          col = alpha(cols["big_atac"], 0.3),
          lty = 0)

  legend("topleft",
         legend = c("CG-DMR", "CH-DMR", "DAR", "bigDAR"),
         col = cols[c("mCG", "mCH", "atac", "big_atac")],
         lty = c(1, 1),
         lwd = c(2, 2))

  # ----------------------------------------------------------------------------
  # Return correlations
  #

  list(mCG = mCG, mCH = mCH, atac = atac, big_atac = big_atac)
}

###---------------------------------------------------------------------------
### Compute and plot for different sized windows around TSSs
###

pdf("../figures/expression_epi_correlation_from_TSS.pdf", height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = pc_transcripts)
saveRDS(val, "../objects/expression_epi_correlation_from_TSS.rds")
dev.off()

# TODO: Stratify genes or transcripts?

# Transcripts with CGI promoters
pdf("../figures/expression_epi_correlation_from_TSS.with_CGI_promoter.pdf",
    height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = pc_transcripts[pc_transcripts$TXNAME %in%
                                          unlist(cgi_promoters$TXNAME)])
saveRDS(val, "../objects/expression_epi_correlation_from_TSS.with_CGI_promoter.rds")
dev.off()

# Transcripts without CGI promoters
pdf("../figures/expression_epi_correlation_from_TSS.without_CGI_promoter.pdf",
    height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = pc_transcripts[!pc_transcripts$TXNAME %in%
                                          unlist(cgi_promoters$TXNAME)])
saveRDS(val, "../objects/expression_epi_correlation_from_TSS.without_CGI_promoter.rds")
dev.off()

###---------------------------------------------------------------------------
### Compute and plot for different sized windows around TSSs (lncRNAs)
###

pdf("../figures/expression_epi_correlation_from_TSS.lncRNA.pdf",
    height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = lnc_transcripts)
saveRDS(val, "../objects/expression_epi_correlation_from_TSS.lncRNA.rds")
dev.off()

# TODO: Stratify genes or transcripts?

# Transcripts with CGI promoters
pdf("../figures/expression_epi_correlation_from_TSS.with_CGI_promoter.lncRNA.pdf",
    height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = lnc_transcripts[lnc_transcripts$TXNAME %in%
                                           unlist(cgi_promoters$TXNAME)])
saveRDS(val,
        "../objects/expression_epi_correlation_from_TSS.with_CGI_promoter.lncRNA.rds")
dev.off()

# Transcripts without CGI promoters
pdf("../figures/expression_epi_correlation_from_TSS.without_CGI_promoter.lncRNA.pdf",
    height = 6, width = 6)
s <- list("0-1000000:1000" = seq(0, 1000000, 1000),
          "0-100000:100" = seq(0, 100000, 100),
          "0-10000:10" = seq(0, 10000, 10))
val <- sapply(s, f, tx = lnc_transcripts[!lnc_transcripts$TXNAME %in%
                                          unlist(cgi_promoters$TXNAME)])
saveRDS(val,
        "../objects/expression_epi_correlation_from_TSS.without_CGI_promoter.lncRNA.rds")
dev.off()
