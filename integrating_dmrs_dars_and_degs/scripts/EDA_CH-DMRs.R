library(GenomicRanges)
load("../objects/assays-and-features.rda")
list_of_candidate_CH_DMRs <-
  readRDS("../../nonCG/objects/list_of_candidate_CH_DMRs.with_meanMeth.rds")
list_of_candidate_CH_DMRs <- lapply(list_of_candidate_CH_DMRs, function(x) {
  x[x$fwer / x$successful_permutations <= 0.05]
})
pos_CA_DMRs <- list_of_candidate_CH_DMRs[[1]]
pos_CA_DMRs$meanDiff <- pos_CA_DMRs$NA_pos - pos_CA_DMRs$BA9_pos
CH_DMRs <- Reduce(union, list_of_candidate_CH_DMRs)
big_dars_pos <- dars_pos[abs(dars_pos$logFC) > 1]

cg_ol <- findOverlaps(dmrs_NAvsBA9pos, degs, ignore.strand = TRUE)
ch_ol <- findOverlaps(pos_CA_DMRs, degs, ignore.strand = TRUE)

# CG-DMRs are strongly hypermethylated, CH-DMRs are hypomethylated, and DEGs
# are evenly split
# This doesn't change much when restricting to DMRs that overlap DEGs
par(mfrow = c(3, 1))
plot(density(dmrs_NAvsBA9pos$meanDiff), main = "CG-DMRs", xlab = "meanDiff")
lines(density(dmrs_NAvsBA9pos$meanDiff[queryHits(cg_ol)]),
      col = "red")
legend("topright", legend = c("All", "ol DEG"), col = c("black", "red"),
       lty = c(1, 1))
abline(v = 0, lty = 2)
plot(density(pos_CA_DMRs$meanDiff), main = "CH-DMRs", xlab = "meanDiff")
lines(density(pos_CA_DMRs$meanDiff[queryHits(ch_ol)]), col = "red")
legend("topright", legend = c("All", "ol DEG"), col = c("black", "red"),
       lty = c(1, 1))
abline(v = 0, lty = 2)
plot(density(degs$logFC), main = "DEGs", ylab = "logFC")
abline(v = 0, lty = 2)

# DEGs with and without CH-DMRs
degs_with_ch_dmr <- subsetByOverlaps(degs, CH_DMRs, ignore.strand = TRUE)
degs_without_ch_dmr <- subsetByOverlaps(degs, CH_DMRs, invert = TRUE,
                                        ignore.strand = TRUE)
# Slight association for presence/absence of CGI promoters
fisher.test(table(overlapsAny(promoters_by_gene[degs$gene_id], CH_DMRs, ignore.strand = TRUE),
                  overlapsAny(promoters_by_gene[degs$gene_id],
                              cgi_promoters_by_gene,
                              ignore.strand = TRUE)))
# Strong association for presence/absence of promoter CG-DMR
fisher.test(table(overlapsAny(promoters_by_gene[degs$gene_id], CH_DMRs, ignore.strand = TRUE),
                  overlapsAny(promoters_by_gene[degs$gene_id],
                              dmrs_NAvsBA9pos,
                              ignore.strand = TRUE)))
# Strong association for presence/absence of promoter DAR
fisher.test(table(overlapsAny(promoters_by_gene[degs$gene_id], CH_DMRs, ignore.strand = TRUE),
                  overlapsAny(promoters_by_gene[degs$gene_id],
                              dars_pos,
                              ignore.strand = TRUE)))
# Strong association for presence/absence of promoter bigDAR
fisher.test(table(overlapsAny(promoters_by_gene[degs$gene_id], CH_DMRs, ignore.strand = TRUE),
                  overlapsAny(promoters_by_gene[degs$gene_id],
                              big_dars_pos,
                              ignore.strand = TRUE)))
# DEGs with CH-DMR are more common in NAcc (but this just reflects bias of
# CH-DMRs being hypomethylated)
# This means that DEGs with CH-DMRs are more likely to be NAcc-specific genes
par(mfrow = c(1, 1))
boxplot(subsetByOverlaps(degs, CH_DMRs,
                         ignore.strand = TRUE)$logFC,
        subsetByOverlaps(degs, CH_DMRs,
                         ignore.strand = TRUE, invert = TRUE)$logFC,
        names = c("CH-DMR", "No CH-DMR"),
        ylab = "logFC")
# DEGs with CH-DMR have larger absolute logFC
par(mfrow = c(1, 1))
boxplot(abs(subsetByOverlaps(degs, CH_DMRs,
                             ignore.strand = TRUE)$logFC),
        abs(subsetByOverlaps(degs, CH_DMRs,
                             ignore.strand = TRUE, invert = TRUE)$logFC),
        names = c("CH-DMR", "No CH-DMR"),
        ylab = "|logFC|")
t.test(abs(subsetByOverlaps(degs, CH_DMRs,
                            ignore.strand = TRUE)$logFC),
       abs(subsetByOverlaps(degs, CH_DMRs,
                            ignore.strand = TRUE, invert = TRUE)$logFC))
pol <- sapply(degs, function(x) {
  sum(width(intersect(x, CH_DMRs, ignore.strand = TRUE))) / width(x)
}) * 100
# The more of the gene covered by a CH-DMR, the greater the |logFC|
par(mfrow = c(2, 1))
scatter.smooth(pol, abs(degs$logFC),
               xlab = "% CH-DMR", ylab = "|logFC|",
               lpars = list(col = "red", lwd = 2))
scatter.smooth(pol[pol > 0], abs(degs$logFC[pol > 0]),
               xlab = "% CH-DMR (> 0)", ylab = "|logFC|",
               lpars = list(col = "red", lwd = 2))
# The more of the gene covered by a CH-DMR, the shorter it is likely to be
par(mfrow = c(2, 1))
scatter.smooth(pol, log10(width(degs)),
               xlab = "% CH-DMR", ylab = "log10(width)",
               lpars = list(col = "red", lwd = 2))
scatter.smooth(pol[pol > 0], log10(width(degs))[pol > 0],
               xlab = "% CH-DMR (> 0)", ylab = "log10(width)",
               lpars = list(col = "red", lwd = 2))

# DEGs with CH-DMRs are longer than those without
par(mfrow = c(1, 1))
boxplot(log10(width(degs_with_ch_dmr)),
        log10(width(degs_without_ch_dmr)),
        names = c("CH-DMR", "No CH-DMR"),
        ylab = "log10(width)",
        main = "DEGs")
