#!/bin/bash

samp=$1
output=$2

name=`basename ${samp} | cut -d "_" -f1`

samtools merge -@8 ${output}/${samp}/${name}.merged.bam *bismark_bt2_pe.bam
wait
samtools sort -n ${output}/${samp}/${name}.merged.bam ${output}/${samp}/${name}.merged.nsorted