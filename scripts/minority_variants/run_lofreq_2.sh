#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/indel_qualities_recalibrated.bam
  do
    bsub.py 8 lofreq_2 -o ${folder} lofreq call -f GCA_015839035.1_ASM1583903v1_genomic.fna --min-cov 20 --call-indels -o ${folder}/lofreq_final.vcf ${bam_file}
  done
done
