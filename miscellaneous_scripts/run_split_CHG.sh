#!/bin/bash

declare -a CHG

for CHG in $(ls /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/*CHG_report.txt)
do
        CHG[${#CHG[@]}+1]=$(echo "${CHG}")
done

for (( i=1; i<=${#CHG[@]}; i++ ))
do
file=`basename ${CHG[i]}`
   qsub -l mem_free=40G,h_vmem=44G,h_fsize=200G -N nonCG_${file} -cwd split_CHG_by_chr.sh ${CHG[i]} 
  
done 

