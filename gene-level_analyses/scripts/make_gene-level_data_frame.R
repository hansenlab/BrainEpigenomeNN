# Construct a data frame with all gene level measurements
# Peter Hickey
# 2018-02-12

### ============================================================================
### Aim
###

# For every GENCODE v19 gene (n = 33,351) want data to compute the following:
#
# [x] Correlation of DiffEpi with gene expression logFC around TSS +/- distance
#     - distance = seq(0, 1000, 1)
#     - distance = seq(0, 10000, 10)
#     - distance = seq(0, 100000, 100)
#     - distance = seq(0, 1000000, 1000)
# [x] Correlation of DiffEpi with gene expression logFC around TES +/- distance
#     - distance = seq(0, 1000, 1)
#     - distance = seq(0, 10000, 10)
#     - distance = seq(0, 100000, 100)
#     - distance = seq(0, 1000000, 1000)
#     - Remember, for CA-DMRs, want opposite strand as well as original strand
# [x] Dinucleotide proportion in 1% bins along gene body +/- 2 GBEs
# [x] Presence of DiffEpi in 1% bins along gene body +/- 2 GBEs
# [x] Absolute mCA (stranded) & mCG (unstranded) in 1% bins along gene body +/- 2 GBEs
#     - NA_pos
#     - BA9_pos
# [x] Correlation with logFC in 1% bins along gene body +/- 2 GBEs
# [x] logFC (NA_pos vs BA9_pos)
# [x] Expression estimate
#     - NA_pos
#     - BA9_pos
# [x] Percentage of gene body covered by DiffEpi
# [x] Average mCA and mCG over gene body
#     - This be computed from binned measurements using a weighted mean, with
#       weights given my number of CX in bin
#
# NOTE: GBE = gene body equivalents, i.e. the length of the gene body

# NOTE: This script should just join up the results of other scripts
