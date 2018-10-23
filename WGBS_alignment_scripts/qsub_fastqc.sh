#!/bin/bash

INPUT_DIR=$1
samp=`basename ${INPUT_DIR}`

qsub -l mem_free=20G,h_vmem=24G,h_fsize=100G -N fastqc_${samp} -cwd /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/fastqc.sh ${INPUT_DIR}