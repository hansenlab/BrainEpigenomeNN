#!/bin/bash
set -e -u
# Quasi-map and quantify bulk RNA-seq data from flow sorted brain GTEx project
# data using Salmon and GENCODE v19 transcript database
# Peter Hickey
# 2015-09-09

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

# NOTE: Request number of CPUs equal to the number required by Salmon plus
#       3 (3 for zcat-ing and seqtk-ing the input FASTQs).

#$ -l mem_free=2G
#$ -l h_vmem=3G
#$ -m n
#$ -pe local 18
#$ -l h_fsize=500G

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

# NOTE: Salmon is not available as a module on JHPCE. Instead, this scripts
#       assumes that Picard is installed and available as `salmon`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

# The table relating sample IDs to FASTQ files and other metadata
SAMPLES_TBL="../extdata/flow-sorted-brain-rna-seq/data/RNA-seq-samples.txt"
# Path of transcript database index. In this case, we want Salmon index for the
# GENCODE v19
IDX="../extdata/flow-sorted-brain-rna-seq/data/GENCODE_v19/gencode.v19.pc_transcripts.lncRNA_transcripts.quasi_index"
# Path to the GTF file accompanying the transcript database. In this case, we
# want the GENCODE v19 GTF
GENEMAP="../extdata/flow-sorted-brain-rna-seq/data/GENCODE_v19/gencode.v19.annotation.gtf.gz"
# Directory for Salmon output
SALMON_DIR="../extdata/flow-sorted-brain-rna-seq/data/salmon"
mkdir -p ${SALMON_DIR}

### =========================================================================
### Trim first 3bp from 5' end of R1, quasi-map reads and quantify with Salmon
### -------------------------------------------------------------------------
###

LINE=$(sed -n "$((${SGE_TASK_ID} + 1))p" ${SAMPLES_TBL})
PREFIX=$(echo "${LINE}" | awk -F "\t"  '{print $1}')
R1s=$(echo "${LINE}" | awk -F "\t"  '{print $2}')
R2s=$(echo "${LINE}" | awk -F "\t"  '{print $3}')
SAMPLE=$(echo "${LINE}" | awk -F "\t"  '{print $4}')
REPLICATE=$(echo "${LINE}" | awk -F "\t"  '{print $5}')

### =========================================================================
### Useful info for debugging
### -------------------------------------------------------------------------
###

echo "hostname"
hostname
echo "SGE_TASK_ID = ${SGE_TASK_ID}"
echo "R1s = ${R1s}"
echo "R2s = ${R1s}"

# NOTE: This uses a named pipe (a.k.a. FIFO) to trim reads and pass them to
#       Salmon on-the-fly
# TODO: --writeMappings?
# TODO: --geneMap argument isn't working with gencode database
salmon quant -i ${IDX} \
				     -l A \
				     --seqBias \
				     --gcBias \
				     --threads 15 \
				     --geneMap ${GENEMAP} \
				     -1 <( seqtk trimfq -b 3 <( zcat ${R1s} ) ) \
				     -2 <( zcat ${R2s} ) \
				     -o ${SALMON_DIR}/${PREFIX}.transcripts_quant
