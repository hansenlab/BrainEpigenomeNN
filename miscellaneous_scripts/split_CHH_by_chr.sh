
#!/bin/bash
#$ -cwd

set -e
set -u
#echo=${1}
   name1=`basename ${1} | cut -d "_" -f1`
   name2=`basename ${1} | cut -d "_" -f2`
   name3=`basename ${1} | cut -d "_" -f3`

  samp=${name1}_${name2}_${name3}
   
   grep "chr1" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr1_CHH_report.txt
   grep "chr2" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr2_CHH_report.txt
   grep "chr3" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr3_CHH_report.txt
   grep "chr4" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr4_CHH_report.txt
   grep "chr5" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr5_CHH_report.txt
   grep "chr6" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr6_CHH_report.txt
   grep "chr7" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr7_CHH_report.txt
   grep "chr8" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr8_CHH_report.txt
   grep "chr9" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr9_CHH_report.txt
   grep "chr10" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr10_CHH_report.txt
   grep "chr11" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr11_CHH_report.txt
   grep "chr12" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr12_CHH_report.txt
   grep "chr13" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr13_CHH_report.txt
   grep "chr14" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr14_CHH_report.txt
   grep "chr15" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr15_CHH_report.txt
   grep "chr16" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr16_CHH_report.txt
   grep "chr17" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr17_CHH_report.txt
   grep "chr18" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr18_CHH_report.txt
   grep "chr19" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr19_CHH_report.txt
   grep "chr20" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr20_CHH_report.txt
   grep "chr21" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr21_CHH_report.txt
   grep "chr22" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chr22_CHH_report.txt
   grep "chrX" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chrX_CHH_report.txt
   grep "chrY" ${1} >> /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/${samp}_chrY_CHH_report.txt

