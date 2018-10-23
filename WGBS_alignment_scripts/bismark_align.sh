#!/bin/bash
## 1 is R1 read, 2 is R2 read, 3 is input dir, 4 output dir, and 5 tmp dir

R1=$1
R2=$2
INPUT_DIR=$3
OUTPUT_DIR=$4
TMP_DIR=$5
samp=$6

#This gets you the sample name, index, lane, and fastq
name=`basename ${R1} | cut -d "_" -f1`
random=`basename ${R1} | cut -d "_" -f2`
lane=`basename ${R1} | cut -d "_" -f3`
num=`basename ${R1} | cut -d "_" -f5 | cut -d "." -f1`

if [[ -s $R1 && -s $R2 ]]
then
    lread1=`ls ${R1}`
    lread2=`ls ${R2}`
    echo "Trimgalore $samp $lread1 $lread2"
else
	echo "FAIL"
    rm -rf $TMP_DIR/$samp
    exit 1
fi

#Run trimgalore
/home/jhmi/sramazan/gtex/trim_galore_v0.4.0/trim_galore -q 25 --paired $lread1 $lread2 \
 --path_to_cutadapt /home/jhmi/sramazan/gtex/cutadapt-1.8.1/bin/cutadapt \
 -o $TMP_DIR/$samp &>$OUTPUT_DIR/$samp/${samp}_${lane}_${num}_trimlog

wait

#Run Bismark
/amber3/feinbergLab/personal/lrizzard/bin/bismark_v0.14.3/bismark --bam --bowtie2 \
    --path_to_bowtie /home/jhmi/sramazan/gtex/bowtie2-2.2.5/ -X 1000 --multicore 4 \
    /amber3/feinbergLab/personal/lrizzard/genomes/hg19/ --temp_dir $TMP_DIR \
    -1 $TMP_DIR/$samp/${name}_${random}_${lane}_R1_${num}_val_1.fq.gz \
    -2 $TMP_DIR/${samp}/${name}_${random}_${lane}_R2_${num}_val_2.fq.gz \
    -o $OUTPUT_DIR/$samp &>$OUTPUT_DIR/$samp/${samp}_${lane}_${num}_galorepair_bismarklog