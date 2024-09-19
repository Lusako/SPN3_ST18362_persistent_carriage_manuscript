#!/bin/bash

for folder in 45*
do
  for fastq_file in ${folder}/samtofastq_interleaved.fq
  do
    bwa mem -M -t 7 -p GCA_015839035.1_ASM1583903v1_genomic.fna ${fastq_file} > ${folder}/bwa_mem.sam
  done
done
