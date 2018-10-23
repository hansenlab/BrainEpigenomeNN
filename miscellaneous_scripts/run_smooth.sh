#!/bin/bash

rawdir=/dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/Sorted_Individual_Objects/

for file in $(ls ${rawdir}*.txt)
do
    filename=`basename ${file}`
    #echo ${filename}
    #echo ${file}
    qsub -l mem_free=80G,h_vmem=81G,h_fsize=200G -N smooth_${filename} -cwd RunRscript.sh ${filename} ${rawdir}

done
