#!/bin/bash
set -e -u
# Run FastQC on bulk RNA-seq data from flow sorted brain GTEx project
# Peter Hickey
# 2016-09-10

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -m n
#$ -pe local 8

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

# NOTE: FastQC is not available as a module on JHPCE. Instead, this scripts
#       assumes that FastQC is installed and available as `fastqc`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-rna-seq/data/RNA-seq-samples.txt"
# Directory for FastQC output files
FASTQC_DIR="../extdata/flow-sorted-brain-rna-seq/data/fastqc"
mkdir -p ${FASTQC_DIR}

### =========================================================================
### Run FastQC on each pair of FASTQs
### -------------------------------------------------------------------------
###

LINE=$(sed -n "$((${SGE_TASK_ID} + 1))p" ${SAMPLES_TBL})
PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

fastqc --outdir ${FASTQC_DIR} \
       --threads 8 \
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
       --threads 8 \
       --noextract \
       ${TMPDIR}/${PREFIX}.R1s.fastq.gz \
       ${TMPDIR}/${PREFIX}.R2s.fastq.gz

rm ${TMPDIR}/${PREFIX}.R1s.fastq.gz ${TMPDIR}/${PREFIX}.R2s.fastq.gz
