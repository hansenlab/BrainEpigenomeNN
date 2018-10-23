#!/bin/bash
set -e
set -u

PROGRAM=$(basename $0)

echo -e "\nYou have initialized $PROGRAM, beginning to run...\n"

usage() {

    echo "$PROGRAM USAGE:"
    echo -e "$PROGRAM -d input_dir -o output_dir -t tmp_dir\n"
    echo -e "\t\t-d input_dir\tThe input directory full path"
    echo -e "\t\t-o output_dir\tThe output directory full path"
    echo -e "\t\t-t tmp_dir\tThe tmp directory full path"
}

if [[ $# -lt 3 ]]; then
    echo "ERROR: Not enough parameters"
    usage
    exit -1
fi

INPUT_DIR=""
OUTPUT_DIR=""
TMP_DIR=""

while getopts ":d:o:t:" optname
  do
  case "$optname" in
      "d")
          INPUT_DIR=$OPTARG
          ;;
      "o")
          OUTPUT_DIR=$OPTARG
          ;;
      "t")
	  TMP_DIR=$OPTARG
	  ;;
      "?")
          echo "ERROR: Unknown option $OPTARG"
          exit -1
          ;;
      ":")
	  echo "ERROR: No argument value for option $OPTARG"
          exit -1
          ;;
      *)
          # Should not occur
          echo"ERROR: Unknown error while processing options"
          exit -1
          ;;
  esac    
done

echo -e "Your input parameters are:"
echo -e "INPUT_DIR is:\t$INPUT_DIR"
echo -e "OUTPUT_DIR is:\t$OUTPUT_DIR"
echo -e "TMP_DIR is:\t$TMP_DIR"

declare -a R1_fastq
declare -a R2_fastq
for R1 in $(ls ${INPUT_DIR}/*R1*)
do 
    R1_fastq[${#R1_fastq[@]}+1]=$(echo "${R1}")
done
for R2 in $(ls ${INPUT_DIR}/*R2*)
do
    R2_fastq[${#R2_fastq[@]}+1]=$(echo "${R2}")
done

samp=`basename ${INPUT_DIR}`

if [[ $samp =~ ^Sample_* ]]
then
    samp=`basename ${INPUT_DIR} | cut -d "_" -f2`
fi

mkdir $OUTPUT_DIR/$samp $TMP_DIR/$samp
cd $OUTPUT_DIR/$samp

for (( i=1; i<=${#R1_fastq[@]}; i++ ))
do
   qsub -l mem_free=20G,h_vmem=24G,h_fsize=100G -pe local 5 -N bismarkalign_${samp} -cwd /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/bismark_align.sh ${R1_fastq[i]} ${R2_fastq[i]} ${INPUT_DIR} ${OUTPUT_DIR} ${TMP_DIR} ${samp}
done

qsub -l mem_free=16G,h_vmem=18G,h_fsize=100G -pe local 8 -hold_jid bismarkalign_${samp} -cwd -N bismarkmerge_${samp} -cwd /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/bam_merge.sh ${samp} ${OUTPUT_DIR}

qsub -l mem_free=16G,h_vmem=18G,h_fsize=200G -pe local 8 -hold_jid bismarkmerge_${samp} -cwd -N bismarkextract_${samp} -cwd /amber3/feinbergLab/personal/sramazan/bash/scripts/gtex/methyl_extract_ignore.sh ${samp} ${OUTPUT_DIR} ${TMP_DIR}