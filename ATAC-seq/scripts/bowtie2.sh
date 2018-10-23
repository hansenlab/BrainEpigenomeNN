#!/bin/bash
set -e -u
# Align bulk ATAC-seq data from flow sorted brain GTEx project data
# Peter Hickey
# 2016-09-10

# NOTE: Requires edits to work on JHPCE

### =========================================================================
### SGE and SLURM variables
### -------------------------------------------------------------------------
###

# NOTE: Not much gained by using more than 8 cores for bowtie2
# NOTE: Request number of CPUs equal to the number required by bowtie2 plus
#       6 (4 for samtools, 2 for trimadap).
# NOTE: The memory requested is probably higher than needs be, but I frequently
#       run into trouble when I try to descrease this for reasons I don't
#       fully understand
# NOTE: The --time argument on MARCC is hard to estimate

#$ -l mem_free=2G
#$ -l h_vmem=2.1G
#$ -m n
#$ -pe local 14
#$ -l h_fsize=500G

# TODO: This are apparently ignored by sbatch; WTF? Have to specify manually,
#       e.g.,
#       sbatch --array=1-27 --exclusive --mem=16G --ntasks=14 \
#              --cpus-per-task=1 --time=24:0:0 bowtie2.sh

#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --ntasks=14
#SBATCH --cpus-per-task=1
#SBATCH --time=24:0:0

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load bowtie2/2.2.5
module load samtools

# NOTE: By loading the bowtie2/2.2.5 module on JHPCE the shell variable
#       $BOWTIE2_INDEXES is created and which points to a directory containing
#       bowtei2 indices for the hg18, hg19, mm9, and mm10 reference genomes.
#       This does **not** happen on MARCC, so will need to specifiy index

# NOTE: trimadap is not available as a module on JHPCE/MARCC. Instead, this
#       scripts assumes that trimadap is installed and available as
#        `trimadap-mt`

# NOTE: Picard is not available as a module on JHPCE/MARCC. Instead, this
#       scripts assumes that trimadap is installed and available as
#        `picard`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-atac/data/ATAC-seq-samples.txt"

# TODO: Make work on both JHPCE and MARCC
# Path of reference genome index. In this case, we want Bowtie2 index for the
# human (hg19) reference genome
# IDX=${BOWTIE2_INDEXES}/hg19
IDX="../extdata/hg19/hg19"

# Directory for logging files
LOG_DIR="../extdata/flow-sorted-brain-atac/logs"
mkdir -p ${LOG_DIR}/bowtie2 ${LOG_DIR}/markdup
# Directory for sorted BAM files, their indices, and flagstat output
BAM_DIR="../extdata/flow-sorted-brain-atac/data/bam"
mkdir -p ${BAM_DIR}

# TODO: Generalise TMPDIR to work on both JHPCE and MARCC
TMPDIR="/scratch/groups/khanse10/tmp/"

### =========================================================================
### Useful info for debugging
### -------------------------------------------------------------------------
###

echo "hostname"
hostname
echo "TMPDIR = ${TMPDIR}"
echo "df -h TMPDIR"
df -h ${TMPDIR}
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### =========================================================================
### Trim adapters from reads, align with bowtie2, and create a sorted, indexed
### BAM file with read groups
### -------------------------------------------------------------------------
###
###

# TODO: Generalise to work on MARCC and JHPCE
LINE=$(sed -n "$((${SLURM_ARRAY_TASK_ID} + 1))p" ${SAMPLES_TBL})

PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

# TODO: Would be good to reduce the number of intermediate files created by
#       samtools sort
# NOTE: This uses named pipes (a.k.a. FIFOs) to trim reads and pass them to
#       bowtie2 on-the-fly
# NOTE: trimadap fills trimmed bases with 'X' characters. These are
#       subsequently trimmed by bowtie2 and these bases are excluded from the
#       output
# NOTE: Since we're aligning concatenated FASTQs, and each FASTQ has its own
#       PU, this is arbitrarily set to '1'. If needed, we can recover the PU
#       from the start of the QNAME.
(bowtie2 -X 2000 \
         --local \
         --threads 8 \
         --time \
         --dovetail \
         --rg-id ${PREFIX} \
         --rg LB:${SAMPLE}.${REPLICATE} \
         --rg PL:Illumina \
         --rg PU:1 \
         --rg SM:${SAMPLE} \
         -x ${IDX} \
         -1 <( trimadap-mt -3 CTGTCTCTTATACACATCTCCGAGCCCACGAGA <( cat ${R1s} ) ) \
         -2 <( trimadap-mt -3 CTGTCTCTTATACACATCTGACGCTGCCGACGA <( cat ${R2s} ) ) \
         | samtools view -u - \
         | samtools sort -o ${BAM_DIR}/${PREFIX}.sorted.bam \
                         -T ${TMPDIR}/${PREFIX} \
                         -O bam \
                         -@ 3 \
                         -m 2G \
                         - ) \
         2> ${LOG_DIR}/bowtie2/${PREFIX}.bowtie2.log
samtools index ${BAM_DIR}/${PREFIX}.sorted.bam

# Mark potential PCR duplicates
# TODO: Make INPUT the piped output of samtools sort
picard MarkDuplicates \
  INPUT=${BAM_DIR}/${PREFIX}.sorted.bam \
  OUTPUT=${BAM_DIR}/${PREFIX}.markdup.bam \
  METRICS_FILE=${LOG_DIR}/markdup/${PREFIX}.markdup.metrics \
  REMOVE_DUPLICATES=false \
  TMP_DIR=${TMPDIR} \
  CREATE_INDEX=true &> \
    ${LOG_DIR}/markdup/${PREFIX}.markdup.log

# Flagstat
samtools flagstat ${BAM_DIR}/${PREFIX}.markdup.bam \
 > ${BAM_DIR}/${PREFIX}.markdup.bam.flagstat
