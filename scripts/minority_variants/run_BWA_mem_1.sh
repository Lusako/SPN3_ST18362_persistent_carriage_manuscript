#!/bin/bash

# Define your arrays
obams=("45907_2#336.bam" "45907_2#293.bam" "45897_1#42.bam" "45897_2#220.bam" "45897_1#211.bam" "45907_2#306.bam" "45907_2#14.bam" "45907_2#100.bam" "45907_1#169.bam" "45907_2#91.bam" "45907_2#155.bam" "45907_1#195.bam" "45907_2#174.bam")
reffas=("GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna" "GCA_015839035.1_ASM1583903v1_genomic.fna")
fq1s=("45907_2#336_1.fastq.gz" "45907_2#293_1.fastq.gz" "45897_1#42_1.fastq.gz" "45897_2#220_1.fastq.gz" "45897_1#211_1.fastq.gz" "45907_2#306_1.fastq.gz" "45907_2#14_1.fastq.gz" "45907_2#100_1.fastq.gz" "45907_1#169_1.fastq.gz" "45907_2#91_1.fastq.gz" "45907_2#155_1.fastq.gz" "45907_1#195_1.fastq.gz" "45907_2#174_1.fastq.gz")
fq2s=("45907_2#336_2.fastq.gz" "45907_2#293_2.fastq.gz" "45897_1#42_2.fastq.gz" "45897_2#220_2.fastq.gz" "45897_1#211_2.fastq.gz" "45907_2#306_2.fastq.gz" "45907_2#14_2.fastq.gz" "45907_2#100_2.fastq.gz" "45907_1#169_2.fastq.gz" "45907_2#91_2.fastq.gz" "45907_2#155_2.fastq.gz" "45907_1#195_2.fastq.gz" "45907_2#174_2.fastq.gz")

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
