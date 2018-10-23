#!/bin/bash

samp=$1
output=$2

qsub -l cegs2=TRUE,mem_free=12G,h_vmem=14G,h_fsize=200G -pe local 8 -cwd /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/methyl_extract_ignore.sh ${samp} ${output}