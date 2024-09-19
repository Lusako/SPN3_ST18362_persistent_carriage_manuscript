#!/bin/bash
for folder in 45*
do
  bsub.py 8 clean_gatk -o ${folder} picard SamToFastq I=${folder}/markilluminaadapters.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | bwa mem -M -t 7 -p 22027_1#37.contigs_velvet.fa /dev/stdin | picard MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=${folder}/reverted.bam OUTPUT=${folder}/Clean_gatk_piped.bam R=22027_1#37.contigs_velvet.fa CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS 
done
