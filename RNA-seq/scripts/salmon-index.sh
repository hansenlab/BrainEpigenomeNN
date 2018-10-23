#!/bin/bash
set -e -u
# Download GENCODE v19 FASTAs and GTFs and
# Peter Hickey
# 2016-09-05

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=16G
#$ -l h_vmem=17G
#$ -m n
#$ -pe local 2

### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

# NOTE: Salmon is not available as a module on JHPCE. Instead, this scripts
#       assumes that Salmon is installed and available as `salmon`

### =========================================================================
### Key parameters as shell variables
### -------------------------------------------------------------------------
###

GENCODE_DIR="../extdata/flow-sorted-brain-rna-seq/data/GENCODE_v19"
mkdir -p ${GENCODE_DIR}

### =========================================================================
### Download FASTAs and GTFs and process them with Salmon index
### -------------------------------------------------------------------------
###

# Download 'Nucleotide sequences of coding transcripts on the reference
# chromosomes' FASTA
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
mv gencode.v19.pc_transcripts.fa.gz ${GENCODE_DIR}/
# Download 'Nucleotide sequences of long non-coding RNA transcripts on the
# reference chromosomes' FASTA
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
mv gencode.v19.lncRNA_transcripts.fa.gz ${GENCODE_DIR}

# Merge together the protein-coding and lncRNA FASTAs
cat ${GENCODE_DIR}/gencode.v19.pc_transcripts.fa.gz \
${GENCODE_DIR}/gencode.v19.lncRNA_transcripts.fa.gz > \
${GENCODE_DIR}/gencode.v19.pc_transcripts.lncRNA_transcripts.fa.gz

# Download 'Comprehensive gene annotation' GTF
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
mv gencode.v19.annotation.gtf.gz ${GENCODE_DIR}/
# Download 'Long non-coding RNA gene annotation' GTF
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz
mv gencode.v19.long_noncoding_RNAs.gtf.gz ${GENCODE_DIR}/

# Build a Salmon index of the merged protein-coding and lncRNAs FASTA
salmon index -t ${GENCODE_DIR}/gencode.v19.pc_transcripts.lncRNA_transcripts.fa.gz \
				-i ${GENCODE_DIR}/gencode.v19.pc_transcripts.lncRNA_transcripts.quasi_index \
				--gencode \
				--p 2 \
				--type quasi \
			 	-k 31
