#!/bin/bash

for folder in 45*
do
  for vcf_file in ${folder}/*_lofreq_final.vcf 
  do
    java -Xmx8g -jar ../snpEff.jar CP046355 ${vcf_file} > ${folder}/lofreq_final.ann.vcf
  done
done
