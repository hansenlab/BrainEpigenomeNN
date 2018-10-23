#!/bin/bash

samp=$1
output=$2
tmp=$3

/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark_methylation_extractor -p --bedGraph --counts --cytosine_report --comprehensive --multicore 8 --gzip \
    --ignore 5 --ignore_r2 10 --genome_folder /amber3/feinbergLab/personal/lrizzard/genomes/hg19 \
    ${output}/${samp}/*.merged.nsorted.bam -o ${output}/${samp} --no_header >>${output}/${samp}/${samp}_extractlog

wait

rm bismarkalign_${samp}* bismarkmerge_${samp}* bismarkextract_${samp}*
rm -rf ${tmp}/${samp}
rm ${output}/${samp}/*fq.gz_bismark_bt2_pe.bam
rm *.o* *.p* *.e* *.po* *.pe*
#rm *_galorepair_bismarklog