#!/bin/bash

#echo ${1}
#echo ${2}
#echo ${3}

#should look like:
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark2bedGraph --ample_memory --CX -o ${3}.CHH ${1}
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark2bedGraph --ample_memory --CX -o ${3}.CHG ${2}
wait
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/coverage2cytosine --ample_memory --CX --genome_folder /amber3/feinbergLab/personal/lrizzard/genomes/hg19 -o ${3}.CHH_report.txt ${3}.CHH.bismark.cov
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/coverage2cytosine --ample_memory --CX --genome_folder /amber3/feinbergLab/personal/lrizzard/genomes/hg19 -o ${3}.CHG_report.txt ${3}.CHG.bismark.cov


