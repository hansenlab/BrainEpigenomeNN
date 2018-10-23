#!/bin/bash
set -e -u
# Align bulk RNA-seq data from flow sorted brain GTEx project data using HISAT2
# Peter Hickey
# 2016-09-11

# NOTE: Requires edits to work on MARCC

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###


# NOTE: Not much gained by using more than 8 cores for HISAT2
# NOTE: Request number of CPUs equal to the number required by bowtie2 plus
#       4 (2 for samtools, 2 for seqtk).
# NOTE: The memory requested is probably higher than needs be, but I frequently
#       run into trouble when I try to descrease this for reasons I don't
#       fully understand
# NOTE: The --time argument on MARCC is hard to estimate

#$ -l mem_free=10G
#$ -l h_vmem=10.1G
#$ -pe local 12

#SBATCH --exclusive
#SBATCH --mem=120G
#SBATCH --n-tasks=12
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load samtools

# NOTE: HISAT2 is not available as a module on JHPCE and the version on MARCC
#       is somewhat outdated. Instead, this script assumes that HISAT2 is
#       installed and available as `hisat2`

# NOTE: seqtk is not available as a module on JHPCE or MARCC. Instead, this
#       scripts assumes that seqtk is installed and available as `seqtk`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-rna-seq/data/RNA-seq-samples.txt"
# Path of reference genome index. In this case, we want HISAT2 index for the
# human (hg19) reference genome
IDX="../extdata/GRCH37/grch37_snp_tran/genome_snp_tran"
# Directory for logging files
LOG_DIR="../extdata/flow-sorted-brain-rna-seq/logs"
mkdir -p ${LOG_DIR}/hisat2
# Directory for sorted BAM files, their indices, and flagstat output
BAM_DIR="../extdata/flow-sorted-brain-rna-seq/data/bam"
mkdir -p ${BAM_DIR}

# TODO: Generalise TMPDIR to work on both JHPCE and MARCC
# TMPDIR="/scratch/groups/khanse10/tmp/"

### =========================================================================
### Trim first 3bp from 5' end of R1, align with HISAT2, and create an
### sorted, indexed BAM file with read groups
### -------------------------------------------------------------------------
###
###

# TODO: Generalise to work on MARCC and JHPCE
LINE=$(sed -n "$((${SGE_TASK_ID} + 1))p" ${SAMPLES_TBL})

PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

# TODO: Would be good to reduce the number of intermediate files created by
#       samtools sort
# NOTE: This uses a named pipe (a.k.a. FIFO) to trim reads and pass them to
#       HISAT2 on-the-fly
# NOTE: Since we're aligning concatenated FASTQs, and each FASTQ has its own
#       PU, this is arbitrarily set to '1'. If needed, we can recover the PU
#       from the start of the QNAME.
(hisat2 --novel-splicesite-outfile ${BAM_DIR}/${PREFIX}.novel-splicesite.txt \
        --threads 8 \
        --time \
        --rg-id ${PREFIX} \
        --rg LB:${SAMPLE}.${REPLICATE} \
        --rg PL:Illumina \
        --rg PU:1 \
        --rg SM:${SAMPLE} \
        --add-chrname \
        -x ${IDX} \
        -1 <( seqtk trimfq -b 3 <( zcat ${R1s} ) ) \
        -2 <( zcat ${R2s} ) \
        | samtools view -u - \
        | samtools sort -o ${BAM_DIR}/${PREFIX}.sorted.bam \
                        -T ${TMPDIR}/${PREFIX} \
                        -O bam \
                        -m 2G \
                        - ) \
       2> ${LOG_DIR}/hisat2/${PREFIX}.hisat2.log

# Index and flagstat
samtools index ${BAM_DIR}/${PREFIX}.sorted.bam
samtools flagstat ${BAM_DIR}/${PREFIX}.sorted.bam \
 > ${BAM_DIR}/${PREFIX}.sorted.bam.flagstat
