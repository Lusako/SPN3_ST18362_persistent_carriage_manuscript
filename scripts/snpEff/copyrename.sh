#!/usr/bin/env bash

##: Do not expand a literal glob * if there are no files/directories
shopt -s nullglob

##: Save the whole path and files in an array.
files=(/Users/lusako/Library/CloudStorage/OneDrive-Malawi-LiverpoolWellcomeResearchProgramme/snpEff/lofreq_GPSC10_ref/45*/output.csv)

for file in "${files[@]}"; do
  IFS='/' read -ra path <<< "$file" ##: split via / and save in an array
  tmp=${path[@]:(-2)}  ##: Remain only the last 2, e.g. 1 input.txt
  new_file=${tmp// /_}   ##: Remove the space, becomes 1input.txt
  cp -v "$file" /Users/lusako/Library/CloudStorage/OneDrive-Malawi-LiverpoolWellcomeResearchProgramme/snpEff/lofreq_GPSC10_ref/"$new_file"
done
