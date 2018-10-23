# Run LDSC with custom functional categories (made in make-annot.R +
# run_LDScore_estimation.R)
# Peter Hickey
# 2018-06-11

# NOTE: On JHPCE, requires `module load python/2.7.6`

library(parallel)
options("mc.cores" = 10)

args <- commandArgs(TRUE)
i <- as.integer(args[1])
message("i = ", i)

categories <- readRDS("../objects/categories.rds")
gwasss <- list.files("../extdata/munge_sumstats/Phase1", full.names = TRUE,
                     pattern = glob2rx("*.sumstats.gz"))

seqlevels <- 1:22

message("category = ", names(categories)[i])

### ============================================================================
### Adjusting for baseline
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (cn == "CNS") {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
                    "../extdata/Phase1/baseline/baseline. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/CNS.", bn, ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/", cn, ".Phase1.,",
                    "../extdata/Phase1/baseline/baseline. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, ".", bn, ".Phase1 ",
                    "--print-coefficients")
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

### ============================================================================
### Marginal analysis (feature + no baseline features except base)
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (cn == "CNS") {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/base.Phase1.,",
                    "../extdata/Phase1/cell_type_groups/CNS. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, "_no_baseline_except_base.", bn,
                    ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/base.Phase1.,",
                    "../output/ldsc/", cn, ".Phase1. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, "_no_baseline_except_base.", bn,
                    ".Phase1 ",
                    "--print-coefficients")
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

### ============================================================================
### Adjusting for baseline + CNS + chromHMM_union + H3K27ac
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (!cn %in% c("CNS", "chromHMM_union", "H3K27ac")) {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../extdata/Phase1/baseline/baseline.,",
                    "../extdata/Phase1/cell_type_groups/CNS.,",
                    "../output/ldsc/chromHMM_union.Phase1.,",
                    "../output/ldsc/H3K27ac.Phase1.,",
                    "../output/ldsc/", cn, ".Phase1. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn,
                    ".adjusting_for_baseline_and_CNS_and_chromHMM_union_and_H3K27ac.",
                    bn,
                    ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("echo Nothing to do for ", cn)
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

### ============================================================================
### Complex set ops features in schizophrenia
###

# CG_only aka CG_3: baseline + {setdiff(X, union(Y, Z))} + {union(X, Y, Z) - {setdiff(X, union(Y, Z))}} + setdiff(POS_CG-DMRs, union(X, Y, Z))
# where X, Y, Z = chromHMM, CNS, H3K27ac
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/baseline/baseline.,",
              "../output/ldsc/chromHMM_only.Phase1.,",
              "../output/ldsc/H3K27ac_only.Phase1.,",
              "../output/ldsc/CNS_only.Phase1.,",
              "../output/ldsc/two_plus.Phase1.,",
              "../output/ldsc/CG_only.Phase1. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CG_only.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# CG_shared aka CG_9: baseline + {setdiff(X, union(Y, Z))} + {union(X, Y, Z) - {setdiff(X, union(Y, Z))}} + intersect(POS_CG-DMRs, union(X, Y, Z))
# TODO: Rename `two_plus_no_CG_shared` files `two_plus_no_CG`
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/baseline/baseline.,",
              "../output/ldsc/chromHMM_only_no_CG.Phase1.,",
              "../output/ldsc/H3K27ac_only_no_CG.Phase1.,",
              "../output/ldsc/CNS_only_no_CG.Phase1.,",
              "../output/ldsc/two_plus_no_CG_shared.Phase1.,",
              "../output/ldsc/CG_only.Phase1.,",
              "../output/ldsc/CG_shared.Phase1. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CG_shared.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

### ============================================================================
### Focus on top traits and see what we can find
### (Schizophrenia, IQ, BMI, ADHD, Neuroticism, Years_of_education,
### Ever_smoked all have Z > 1.96; Depressive_symptoms, Epilepsy,
### College_attainment all have Z > 1.8)
###

# ------------------------------------------------------------------------------
# Schizophrenia-specific stuff
#

# CNS + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_POS_CG-DMRs.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# H3K27ac + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/H3K27ac_and_POS_CG-DMRs.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# POS_CG-DMRs vs. POSvsNEG_CG-DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../output/ldsc/POSvsNEG_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/POS_CG-DMRs_and_POSvsNEG_CG-DMRs.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# CNS + H3K27ac + chromHMM_union + overall_OCRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/overall_OCRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_overall_OCRss.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# POS_DARs + overall_OCRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../output/ldsc/POS_DARs.Phase1.,",
              "../output/ldsc/overall_OCRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/POS_DARs_and_OCRs.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# overall_OCRs + POS_CG-DMRs + POSvsNEG_CG-DMRs + POS_CG-blocks +
# POSvsNEG_CG-blocks + POS_DARs + POSvsNEG_DARs + CH-DMRs + H3K27ac +
# chromHMM_union + Luo_small_CG-DMRs + CNS
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../output/ldsc/overall_OCRs.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../output/ldsc/POSvsNEG_CG-DMRs.Phase1.,",
              "../output/ldsc/POS_CG-blocks.Phase1.,",
              "../output/ldsc/POSvsNEG_CG-blocks.Phase1.,",
              "../output/ldsc/POS_DARs.Phase1.,",
              "../output/ldsc/POSvsNEG_DARs.Phase1.,",
              "../output/ldsc/CH-DMRs.Phase1.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/Luo_small_CG-DMRs.Phase1.,",
              "../extdata/Phase1/cell_type_groups/CNS.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/the_works.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# POS_CG-DMRs no baseline features except base
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../output/ldsc/base.Phase1. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/POS_CG-DMRs_no_baseline.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# CNS no baseline features except base
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Schizophrenia.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/base.Phase1. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_no_baseline.Schizophrenia.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# IQ-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/IQ.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.IQ.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# BMI-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/BMI.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.BMI.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# ADHD-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/ADHD.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.ADHD.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# Neuroticism-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Neuroticism.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.Neuroticism.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# Years_of_education-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Years_of_education.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.Years_of_education.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# ------------------------------------------------------------------------------
# Depressive_symptoms-specific stuff
#

# CNS + H3K27ac + chromHMM_union + POS_CG_DMRs
cmd <- paste0("python ",
              "/users/phickey/software/ldsc/ldsc.py ",
              "--h2 ../output/munge_sumstats/Phase1/Depressive_symptoms.sumstats.gz ",
              "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
              "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
              "../output/ldsc/H3K27ac.Phase1.,",
              "../output/ldsc/chromHMM_union.Phase1.,",
              "../output/ldsc/POS_CG-DMRs.Phase1.,",
              "../extdata/Phase1/baseline/baseline. ",
              "--overlap-annot ",
              "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
              "--out ../output/ldsc/CNS_and_H3K27ac_and_chromHMM_union_and_POS_CG-DMRs.Depressive_symptoms.Phase1 ",
              "--print-coefficients")
print(cmd)
system(cmd)

# TODO: Choose final categories
