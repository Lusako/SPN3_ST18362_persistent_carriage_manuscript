#!/bin/bash

# Define your arrays
obams=( "bam_file")
reffas=("GCA_015839035.1_ASM1583903v1_genomic.fna")
fq1s=("1.fastq_file")
fq2s=("2.fastq_file")

# Loop through each set of files
for i in "${!obams[@]}"; do
  obam="${obams[$i]}"
  reffa="${reffas[$i]}"
  fq1="${fq1s[$i]}"
  fq2="${fq2s[$i]}"

  echo "Processing $obam with $reffa, $fq1, and $fq2"

  bwa mem -M -t 7 "$reffa" "$fq1" "$fq2" | \
    samtools fixmate -m - - | \
    lofreq viterbi -f "$reffa" - | \
    samtools sort - | \
    samtools view -S -b -o "$obam"
done
