module load ..
bwa/0.7.17-r1188
samtools/1.9--h91753b0_8
picard/2.22.2--0
gatk/4.1.4.1
conda
conda activate lofreq

bwa index -a bwtsw GCA_015839035.1_ASM1583903v1_genomic.fna

Creating the fasta index file
samtools faidx GCA_015839035.1_ASM1583903v1_genomic.fna

./run_BWA_mem_1.sh

./run_makedir_and_movefiles.sh

Creating the FASTA sequence dictionary file
picard CreateSequenceDictionary R=ref.fna O=ref.dict

./run_CleanSam.sh

./run_lofreq.sh

The average genome coverage values aren’t that high. It’s usually recommended to use an ad-hoc minimum coverage
of 10*1/minAF (below is a quote from Adam Lauring’s review paper). The minAF is usually 0.03 (3%), but in this case 
maybe you can increase it to 0.02 (2%), which would give you a min coverage depth for lofreq as 20.
"When the viral input is adequate, the precision of frequency estimates generally scales with depth of coverage. 
A good rule of thumb is that the coverage should be 10 times the reciprocal of a variant's frequency (i.e., 200-fold 
coverage to reliably estimate variants present at ≥5% frequency)." 
https://www.annualreviews.org/content/journals/10.1146/annurev-virology-010320-061642

GATK’s best practice protocol, i.e. that you mark duplicates (not for very high coverage data though), realign indels 
and recalibrate base qualities with GATK (BQSR).
https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
http://broadinstitute.github.io/picard/command-line-overview.html#IlluminaBasecallsToSam
https://gatk.broadinstitute.org/hc/en-us/articles/13832708374939-BaseRecalibrator
https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

./run_reverted.sh

./run_MarkIlluminaAdapters.sh

./run_bamtofstq.sh

./run_bwa_mem.sh

Check total number of threads of the server 
getconf _NPROCESSORS_ONLN - sanger = 64

./run_MergeBamAlignment.sh

warning: refre to - https://github.com/broadinstitute/picard/issues/1712

https://gatk.broadinstitute.org/hc/en-us/articles/360036729891-MarkDuplicates-Picard-
./run_MarkDuplicates.sh

https://gatk.broadinstitute.org/hc/en-us/articles/360036804671-CollectWgsMetrics-Picard-
./run_CollectWgsMetrics.sh

Index for gatk recalibration
./run_IndexFeatureFile.sh

#Create Index for BAM File - Already created picard MergeBamAlignment with CREATE_INDEX=true
#Create an index for the deduplicated BAM file.
#samtools index -bc marked_duplicates.bam

samtools view -H marked_duplicates.bam | grep '@RG'

./run_picard_AddOrReplaceReadGroups.sh

samtools view -H marked_duplicates_group.bam | grep '@RG'

./run_BaseRecalibrator.sh

./run_ApplyBQSR.sh

Please check for indels with: samtools view recalibrated.bam | grep 'B[ID]:Z:' 

./run_lofreq_indel.sh

Please check for indels with: samtools view indel_qualities_recalibrated.bam | grep 'B[ID]:Z:'

#alnqual: Insert base and indel alignment qualities
#Rarely needed because computed on the fly in lofreq call

./run_lofreq_2.sh

note: WARNING(lofreq_call.c|main_call): 8 indel calls (before filtering) were made without indel alignment-quality! 
Did you forget to indel alignment-quality to your bam-file?
