#!/bin/bash

INPUT_DIR=$1

declare -a samples
for S in $(ls ${INPUT_DIR}/)
do 
    samples[${#samples[@]}+1]=$(echo "${S}")
done

for (( i=1; i<=${#samples[@]}; i++ ));do
    /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/qsub_BSconRate.sh ${INPUT_DIR}/${samples[i]}
done
