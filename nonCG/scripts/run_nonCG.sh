#!/bin/bash



INPUT_DIR=/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG


declare -a CHH

for CHH in $(ls ${INPUT_DIR}/CHH_context_*)
do
	CHH[${CHH[@]}+1]=$(echo "${CHH}")
done

declare -a CHG

for CHG in $(ls ${INPUT_DIR}/CHG_context_*)
do
	CHG[${CHG[@]}+1]=$(echo "${CHG}")
done
samp=`basename ${INPUT_DIR}`

for (( i=1; i<=${#CHH[@]}; i++ ))
do
   qsub -l cegs2=TRUE,mem_free=20G,h_vmem=24G,h_fsize=500G -pe local 8 -N nonCG_${samp} -cwd nonCG.sh ${CHH[i]} ${CHG[i]} ${samp}
done
