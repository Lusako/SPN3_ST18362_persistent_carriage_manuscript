#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/recalibrated.bam 
  do
    bsub.py 8 lofreq_indel -o ${folder} -e ${folder} lofreq indelqual --dindel --ref GCA_015839035.1_ASM1583903v1_genomic.fna --out ${folder}/indel_qualities_recalibrated.bam ${bam_file}
  done
done
