#!/bin/bash

INPUT_DIR=$1

declare -a R1_fastq
declare -a R2_fastq
for R1 in $(ls ${INPUT_DIR}/*R1*)
do
    R1_fastq[${#R1_fastq[@]}+1]=$(echo "${R1}")
done
for R2 in $(ls ${INPUT_DIR}/*R2*)
do
    R2_fastq[${#R2_fastq[@]}+1]=$(echo "${R2}")
done

name=`basename ${R1_fastq[1]}`
samp=`basename ${name} | cut -d "_" -f1`

mkdir /dcl01/feinberg/data/gtex/gtex/GTEX_conversion/${samp}
cd /dcl01/feinberg/data/gtex/gtex/GTEX_conversion/${samp}

for (( i=1; i<=${#R1_fastq[@]}; i++ ))
do
    qsub -l cegs=TRUE,mem_free=5G,h_vmem=8G,h_fsize=100G -cwd -N bismarkalign_BSconversionRate_${samp} /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/GTEX_conversion.sh ${R1_fastq[i]} ${R2_fastq[i]} ${samp}
done