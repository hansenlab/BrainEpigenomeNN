---
title: "nonCG Notes"
author: "Lindsay Rizzardi"
date: "October 2, 2015"
output: html_document
---

First need to perform bismark2bedGraph with --CX flag then run coverage2cytosine also with --CX flag using as input all three files below:
     CHG_context_5248_BA9_neg.all.nsorted.txt.gz
     CHH_context_5248_BA9_neg.all.nsorted.txt.gz

Test with one sample first:
     qrsh -l mem_free=22G,h_vmem=24G,h_fsize=500G -pe local 8
In the Methylation Folder:
     mkdir nonCG
Run on one file at a time:
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark2bedGraph --ample_memory -o 5248_BA9_neg.all.nsorted.CHG --CX CHG_context_5248_BA9_neg.all.nsorted.txt.gz

THEN:
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/coverage2cytosine --genome_folder /amber3/feinbergLab/personal/lrizzard/genomes/hg19 -o 5248_BA9_neg --CX 5248_BA9_neg.all.nsorted.CHG.bismark.cov

create scripts to run qsub on each Sample:

perl /amber3/feinbergLab/personal/lrizzard/BSseq_Amy/qsub_sge.pl -cwd HiSeq140_bismark2bedgraph.sh -l cegs2=TRUE,mem_free=20G,h_vmem=24G,h_fsize=500G -pe 8

move all files into new folder:
     cd /dcl01/feinberg/data/gtex/flow_sorted_nuclei/
     rsync -avc */*/*.nsorted.txt.gz /dcl01/feinberg/data/gtex/flow_sorted_nuclei/methylation/nonCG/fornonCG/
then run the run_nonCG.sh from the nonCG directory

Core dumped : try one CHH file with max memory and file size, try this next time (10/2/15) mf=200G vm=50G RAM and file size 800G -pe 4
qsub -l cegs,mf=80G,h_vmem=85G,h_fsize=800G   /amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark2bedGraph --CX --ample_memory -o CHH_context_5628_HC_pos CHH_context_5628_HC_pos.all.nsorted.txt.gz

so ran with more:Problem was way I was using Xin’s script! format as below! was only using 10G of RAM before

perl qsub_sge.pl align_nonCG.sh --resource cegs,mf=40G,h_vmem=45G,h_fsize=800G

four didn’t finish so run individually

perl qsub_sge.pl last4.sh --resource cegs,mf=80G,h_vmem=85G,h_fsize=800G

Go ahead and run cov2cyt without those last 4

perl qsub_sge.pl cov2cyt_nonCG.sh --resource cegs,mf=40G,h_vmem=45G,h_fsize=800G

When making/working with the bsseq objects make sure you request a TON of ram!!! each file is 32G
will make an object for each sample and then merge them together.



