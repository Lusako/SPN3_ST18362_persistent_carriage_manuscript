#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/marked_duplicates.bam 
  do
    bsub.py 8 picard_AddOrReplaceReadGroups -o ${folder} picard AddOrReplaceReadGroups I=${bam_file} O=${folder}/marked_duplicates_group.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
  done
done
