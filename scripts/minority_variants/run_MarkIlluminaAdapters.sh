#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/reverted.bam
  do
    bsub.py 8 MarkIlluminaAdapters -o ${folder} picard MarkIlluminaAdapters I=${bam_file} O=${folder}/markilluminaadapters.bam M=${folder}/markilluminaadapters_metrics.txt 
  done
done
