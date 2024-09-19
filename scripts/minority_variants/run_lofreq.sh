#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/*_clean.bam
  do
    bsub.py 8 lofreq -o ${folder} lofreq call -f GCA_015839035.1_ASM1583903v1_genomic.fna -o ${folder}/lofreq.vcf --min-cov 20 ${bam_file}
  done
done
