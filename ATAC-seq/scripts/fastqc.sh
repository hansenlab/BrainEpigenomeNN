#!/bin/bash
set -e -u
# Run FastQC on bulk ATAC-seq data from flow sorted brain GTEx project
# Peter Hickey
# 2016-09-10

# NOTE: Requires edits to work on JHPCE


### =========================================================================
### SGE and SLURM variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -m n
#$ -pe local 16

#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --n-tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:0:0

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

# NOTE: FastQC is not available as a module on JHPCE or MARCC. Instead, this
#       scripts assumes that FastQC is installed and available as `fastqc`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-atac/data/ATAC-seq-samples.txt"
# Directory for FastQC output files
FASTQC_DIR="../extdata/flow-sorted-brain-atac/data/fastqc"
mkdir -p ${FASTQC_DIR}

# TODO: Generalise TMPDIR to work on both JHPCE and MARCC
TMPDIR="/scratch/groups/khanse10/tmp/"

### =========================================================================
### Run FastQC on each pair of FASTQs
### -------------------------------------------------------------------------
###

LINE=$(sed -n "$((${SLURM_ARRAY_TASK_ID} + 1))p" ${SAMPLES_TBL})
PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

fastqc --outdir ${FASTQC_DIR} \
       --threads 16 \
       --noextract \
       ${R1s} \
       ${R2s}

### =========================================================================
### Run FastQC on all pairs of FASTQs from each replicate
### -------------------------------------------------------------------------
###

# NOTE: FastQC doesn't work with FIFOs, so need to create temporary files
cat ${R1s} > ${TMPDIR}/${PREFIX}.R1s.fastq.gz
cat ${R2s} > ${TMPDIR}/${PREFIX}.R2s.fastq.gz

fastqc --outdir ${FASTQC_DIR} \
       --threads 16 \
       --noextract \
       ${TMPDIR}/${PREFIX}.R1s.fastq.gz \
       ${TMPDIR}/${PREFIX}.R2s.fastq.gz

rm ${TMPDIR}/${PREFIX}.R1s.fastq.gz ${TMPDIR}/${PREFIX}.R2s.fastq.gz
