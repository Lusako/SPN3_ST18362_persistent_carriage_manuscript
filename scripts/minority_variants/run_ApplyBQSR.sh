#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/marked_duplicates_group.bam
  do
    gatk ApplyBQSR -R GCA_015839035.1_ASM1583903v1_genomic.fna -I ${bam_file} --bqsr-recal-file ${folder}/recalibration.table -O ${folder}/recalibrated.bam
  done
done
