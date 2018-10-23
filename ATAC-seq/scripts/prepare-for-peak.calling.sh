# Prepare bulk ATAC-seq data from flow sorted brain GTEx project for peak
# calling with MACS2
# Peter Hickey
# 2016-09-12

# We want to create an aggregated BAM file using data from all samples. Each
# sample's data is filtered to remove the following:
# 1. chrM reads (don't care about peaks on chrM. Because so many reads come
# from chrM, removing these greatly reduces the size of the final file).
# 2. Remove duplicate reads (this is done in the original publication prior to
# peak calling).
# 3. Require mapQ >= 30 (done by Buenrostro et al. in scATAC-seq paper)

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=4G
#$ -l h_vmem=5G
#$ -m n
#$ -l h_fsize=500G

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load samtools

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-atac/data/ATAC-seq-samples.txt"
# Directory containing sorted BAM files
BAM_DIR="../extdata/flow-sorted-brain-atac/data/bam"

### =========================================================================
### Create filtered individual BAM files in parallel and then merge
### -------------------------------------------------------------------------
###

# TODO: Generalise to work on MARCC and JHPCE
LINE=$(sed -n "$((${SGE_TASK_ID} + 1))p" ${SAMPLES_TBL})

PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

samtools view -b \
              -o ${BAM_DIR}/${PREFIX}.filtered.bam \
              -F 1024 \
              -q 30 \
              ${BAM_DIR}/${PREFIX}.markdup.bam \
              chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 \
              chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
              chr21 chr22 chrX chrY
