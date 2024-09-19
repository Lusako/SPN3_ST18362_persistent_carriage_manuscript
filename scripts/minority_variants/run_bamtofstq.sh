#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/markilluminaadapters.bam
  do
    bsub.py 8 BamtoFatq -o ${folder} picard SamToFastq I=${bam_file} FASTQ=${folder}/samtofastq_interleaved.fq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true 
  done
done
