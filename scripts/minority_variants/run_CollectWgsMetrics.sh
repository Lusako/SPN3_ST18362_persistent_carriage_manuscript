#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/marked_duplicates.bam 
  do
    bsub.py 8 CollectWgsMetrics -o ${folder} picard CollectWgsMetrics I=${bam_file} O=${folder}/collect_wgs_metrics.txt R=GCA_015839035.1_ASM1583903v1_genomic.fna
  done
done
