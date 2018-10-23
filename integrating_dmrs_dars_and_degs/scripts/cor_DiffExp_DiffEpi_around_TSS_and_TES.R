# Compute and plot correlation of DiffEpi marks with DiffExp around TSS and
# TESs of protein coding genes
# Peter Hickey
# 2018-02-01

tes <- readRDS("../objects/expression_epi_correlation_from_TES.rds")
tss <- readRDS("../objects/expression_epi_correlation_from_TSS.rds")

pdf("../figures/expression_epi_correlation_from_TSS_and_TES.pdf")
lapply(seq_len(ncol(tss)), function(j) {
  par(mfrow = c(2, 2))
  lapply(seq_len(nrow(tss)), function(i) {
    tss_ <- tss[i, j][[1]]
    tes_ <- tes[i, j][[1]]
    plot(cor ~ x,
         data = tss_,
         type = "l",
         xlab = "Distance from TSS/TES",
         main = rownames(tss)[i],
         ylim = c(0, 1))
    lines(cor ~ x, data = tes_, col = "red")
    legend("topleft",
           legend = c("TSS", "TES"),
           col = c("black", "red"),
           lty = 1)
    abline(v = 0, lty = 2)
  })
})
dev.off()
