#!/bin/bash
#$ -N smooth
#$ -cwd
#$ -m e
#$ -M lindsay.rizzardi@jhmi.edu

filename=$1
rawdir=$2
Rscript SmallSmooth.R ${filename} ${rawdir} ${filename}.Rout

