# Run LD score estimation with functional categories (made in make-annot.R)
# Peter Hickey
# 2018-06-11

# NOTE: On JHPCE, requires `module load python/2.7.6`

library(parallel)
options("mc.cores" = 20)

args <- commandArgs(TRUE)
i <- as.integer(args[1])
message("i = ", i)

categories <- readRDS("../objects/complex_set_op_features.rds")

seqlevels <- 1:22

message("category = ", names(categories)[i])

### ============================================================================
### LD Score estimation
###

lapply(names(categories)[i], function(cn) {
  message(cn)
  mclapply(seqlevels, function(sl) {
    cmd <- paste0("python ",
                  "/users/phickey/software/ldsc/ldsc.py ",
                  "--l2 ",
                  "--bfile ../extdata/Phase1/1000G_plinkfiles/1000G.mac5eur.",
                  sl, " ",
                  "--ld-wind-cm 1 ",
                  "--annot ../output/ldsc/", cn, ".Phase1.", sl, ".annot.gz ",
                  "--out ../output/ldsc/", cn, ".Phase1.", sl, " ",
                  "--print-snps ../extdata/Phase1/hapmap3_snps/hm.", sl, ".snp")
    print(cmd)
    system(cmd)
  }, mc.cores = getOption("mc.cores"))
})

