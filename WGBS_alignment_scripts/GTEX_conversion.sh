#!/bin/bash

R1_fastq=$1
R2_fastq=$2
samp=$3

/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark --bam --bowtie2 --path_to_bowtie /home/jhmi/sramazan/gtex/bowtie2-2.2.5/ -X 1000 \
    /amber3/feinbergLab/personal/lrizzard/TAB_genome/ --temp_dir /dcl01/feinberg/data/gtex/gtex/GTEX_conversion/tmp \
    -1 ${R1_fastq} -2 ${R2_fastq} \
    -o /dcl01/feinberg/data/gtex/gtex/GTEX_conversion/${samp} &>/dcl01/feinberg/data/gtex/gtex/GTEX_conversion/${samp}/${samp}_bismarklog