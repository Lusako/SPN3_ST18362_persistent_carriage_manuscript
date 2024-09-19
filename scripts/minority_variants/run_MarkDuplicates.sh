#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/mergebamalignment.bam 
  do
    bsub.py 8 MarkDuplicates -o ${folder} picard MarkDuplicates I=${bam_file} O=${folder}/marked_duplicates.bam M=${folder}/marked_dup_metrics.txt REMOVE_SEQUENCING_DUPLICATES=true
  done
done
