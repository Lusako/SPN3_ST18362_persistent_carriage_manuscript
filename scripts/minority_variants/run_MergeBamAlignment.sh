#!/bin/bash

for folder in 45*
do
  for bam_file in ${folder}/reverted.bam
  do
    bsub.py 8 MergeBamAlignment -o ${folder} picard MergeBamAlignment R=GCA_015839035.1_ASM1583903v1_genomic.fna UNMAPPED_BAM=${bam_file} ALIGNED_BAM=${folder}/bwa_mem.sam O=${folder}/mergebamalignment.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS
  done
done
