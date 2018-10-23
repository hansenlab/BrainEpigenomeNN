#!/bin/bash

INPUT_DIR=$1
OUTPUT_DIR=$2
TMP=$3

declare -a samples

for S in $(ls ${INPUT_DIR}/)
do 
    samples[${#samples[@]}+1]=$(echo "${S}")
    
done

for (( i=1; i<=${#samples[@]}; i++ ));do
    /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/run_align_SHARED.sh -d ${INPUT_DIR}/${samples[i]} -o ${OUTPUT_DIR} -t ${TMP}
done
