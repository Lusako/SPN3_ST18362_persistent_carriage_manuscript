#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/marked_duplicates_group.bam
  do
    gatk BaseRecalibrator -I ${bam_file} -R GCA_015839035.1_ASM1583903v1_genomic.fna --known-sites ${folder}/lofreq.vcf -O ${folder}/recalibration.table
  done
done
