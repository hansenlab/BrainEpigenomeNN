# Call peaks in bulk ATAC-seq data from flow sorted brain GTEx project
# Peter Hickey
# 2016-09-12

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=8G
#$ -l h_vmem=9G
#$ -m n
#$ -l h_fsize=500G
#$ -pe local 7

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load macs/2.1.0

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-atac/data/ATAC-seq-samples.txt"

OUT_DIR="../extdata/flow-sorted-brain-atac/data/macs2"
mkdir -p ${OUT_DIR}
# Directory for logging files
LOG_DIR="../extdata/flow-sorted-brain-atac/logs"
mkdir -p ${LOG_DIR}/macs2
# Directory containing filtered BAM files
BAM_DIR="../extdata/flow-sorted-brain-atac/data/bam"

# Create bash arrays with the various subests of BAM files
mapfile -t PREFIX < <( tail -n +2 ${SAMPLES_TBL} | cut -f1 )
# Append and prepend to bash array
# (http://web.archive.org/web/20101114051536/http://codesnippets.joyent.com/posts/show/1826)
BAMS=( "${PREFIX[@]/%/.filtered.bam}" )
BAMS=( "${BAMS[@]/#/${BAM_DIR}/}")
POS=($(grep "pos" <(printf "%s\n" "${BAMS[@]}")))
NEG=($(grep "neg" <(printf "%s\n" "${BAMS[@]}")))
NA_POS=($(grep "NA" <(printf "%s\n" "${POS[@]}")))
NA_NEG=($(grep "NA" <(printf "%s\n" "${NEG[@]}")))
BA9_POS=($(grep "BA9" <(printf "%s\n" "${POS[@]}")))
BA9_NEG=($(grep "BA9" <(printf "%s\n" "${NEG[@]}")))

### =========================================================================
### Run MACS2 on aggregated data since all samples are from LCLs
### -------------------------------------------------------------------------
###

# NOTE: Using macs2 options given in Buenrostro, J. D. et al. Single-cell
#       chromatin accessibility reveals principles of regulatory variation.
#       Nature 523, 486–490 (2015) for **scaATAC-seq**. The original ATAC-seq
#       paper, Buenrostro, J. D., Giresi, P. G., Zaba, L. C., Chang, H. Y. &
#       Greenleaf, W. J. Transposition of native chromatin for fast and sensitive
#       epigenomic profiling of open chromatin, DNA-binding proteins and
#       nucleosome position. Nat. Methods 10, 1213–1218 (2013), used ZINBA to
#       call peaks
# NOTE: Run each process in the background and wait to collect results

# Call peaks using all samples (overall)
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.overall \
  -t ${BAMS[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.overall.callpeaks.log &

# Call peaks separately in NeuN+
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.pos \
  -t ${POS[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.pos.callpeaks.log &

# Call peaks separately in NeuN+
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.neg \
  -t ${NEG[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.neg.callpeaks.log &

# Call peaks separately in NA-pos
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.NA-pos \
  -t ${NA_POS[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.NA-pos.callpeaks.log &

# Call peaks separately in NA-neg
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.NA-neg \
  -t ${NA_NEG[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.NA-neg.callpeaks.log &

# Call peaks separately in BA9-pos
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.BA9-pos \
  -t ${BA9_POS[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.BA9-pos.callpeaks.log &

# Call peaks separately in BA9-neg
macs2 callpeak \
  --nomodel \
  --nolambda \
  --keep-dup all \
  --call-summits \
  --outdir ${OUT_DIR} \
  --name flow-sorted-brain-atac.BA9-neg \
  -t ${BA9_NEG[@]} \
  &> ${LOG_DIR}/macs2/flow-sorted-brain-atac.BA9-neg.callpeaks.log &

wait
