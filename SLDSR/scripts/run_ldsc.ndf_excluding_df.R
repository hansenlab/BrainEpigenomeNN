# Run LDSC with  'non differential features (ndf) excluding differential
# features (df)' categories (made in make-annot.R)
# Peter Hickey
# 2018-06-11

# NOTE: On JHPCE, requires `module load python/2.7.6`

library(parallel)
options("mc.cores" = 40)

categories <- readRDS("../objects/categories.rds")
categories <- categories[c("POS_CG-DMRs", "POSvsNEG_CG-DMRs",
                           "CH-DMRs",
                           "POS_DARs", "POSvsNEG_DARs",
                           "POS_bigDARs", "POSvsNEG_bigDARs")]
gwasss <- list.files("../extdata/munge_sumstats/Phase1", full.names = TRUE,
                     pattern = glob2rx("*.sumstats.gz"))

seqlevels <- 1:22

### ============================================================================
### Adjusting for baseline + setdiff({CNS, chromHMM_union, H3K27ac},
### differential feature)
###

mclapply(names(categories), function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    cmd <- paste0("python ",
                  "/users/phickey/software/ldsc/ldsc.py ",
                  "--h2 ", x, " ",
                  "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                  "--ref-ld-chr ../extdata/Phase1/baseline/baseline.,",
                  paste0("../output/ldsc/CNS_excluding_",  cn,
                         ".Phase1.,"),
                  paste0("../output/ldsc/chromHMM_union_excluding_",  cn,
                         ".Phase1.,"),
                  paste0("../output/ldsc/H3K27ac_excluding_",  cn,
                         ".Phase1.,"),
                  "../output/ldsc/", cn, ".Phase1. ",
                  "--overlap-annot ",
                  "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                  "--out ../output/ldsc/", cn,
                  ".adjusting_for_ndf_excluding_df.",
                  bn,
                  ".Phase1 ",
                  "--print-coefficients")
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
}, mc.cores = 4)
