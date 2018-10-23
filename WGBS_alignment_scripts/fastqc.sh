#!/bin/bash

OUT_DIR=$1
/amber3/feinbergLab/personal/sramazan/bin/FastQC/fastqc --outdir ${OUT_DIR} --extract --format bam --dir /dcl01/feinberg/data/gtex/gtex/GTEX_conversion/tmp ${OUT_DIR}/*merged.nsorted.bam