#!/bin/bash

for folder in 45*
do
  for vcf_file_ann in ${folder}/lofreq_final.ann.vcf 
  do
    cat ${vcf_file_ann} | ../scripts/vcfEffOnePerLine.pl | java -jar ../SnpSift.jar extractFields - CHROM POS REF ALT DP AF SB "ANN[*].GENE" "ANN[*].EFFECT" > ${folder}/output.csv
  done
done
