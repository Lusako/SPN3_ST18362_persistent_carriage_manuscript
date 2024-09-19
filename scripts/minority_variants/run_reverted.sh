#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/*_clean.bam 
  do
    bsub.py 8 reverted -o ${folder} picard RevertSam I=${bam_file} O=${folder}/reverted.bam
  done
done
