#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/45* 
  do
    bsub.py 8 CleanSam -o ${folder} picard CleanSam I=${bam_file} O=${folder}/${folder}_clean.bam
  done
done
