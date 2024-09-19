#!/bin/bash

for folder in 45*
do
  for vcf_file in ${folder}/lofreq.vcf 
  do
    bsub.py 8 IndexFeatureFile -o ${folder} gatk IndexFeatureFile -I ${vcf_file}
  done
done
