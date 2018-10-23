#!/bin/bash
#$ -cwd

set -e
set -u
#echo=${1}
   name1=`basename ${1} | cut -d "_" -f1`
   name2=`basename ${1} | cut -d "_" -f2`
   name3=`basename ${1} | cut -d "_" -f3`
   
  
   samp=${name1}_${name2}_${name3}
   grep "chr1" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr1_CHG_report.txt
   grep "chr2" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr2_CHG_report.txt
   grep "chr3" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr3_CHG_report.txt
   grep "chr4" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr4_CHG_report.txt
   grep "chr5" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr5_CHG_report.txt
   grep "chr6" ${1]} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr6_CHG_report.txt
   grep "chr7" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr7_CHG_report.txt
   grep "chr8" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr8_CHG_report.txt
   grep "chr9" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr9_CHG_report.txt
   grep "chr10" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr10_CHG_report.txt
   grep "chr11" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr11_CHG_report.txt
   grep "chr12" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr12_CHG_report.txt
   grep "chr13" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr13_CHG_report.txt
   grep "chr14" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr14_CHG_report.txt
   grep "chr15" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr15_CHG_report.txt
   grep "chr16" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr16_CHG_report.txt
   grep "chr17" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr17_CHG_report.txt
   grep "chr18" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr18_CHG_report.txt
   grep "chr19" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr19_CHG_report.txt
   grep "chr20" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr20_CHG_report.txt
   grep "chr21" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr21_CHG_report.txt
   grep "chr22" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr22_CHG_report.txt
   grep "chrX" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chrX_CHG_report.txt
   grep "chrY" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chrY_CHG_report.txt


