#!/bin/bash

# Here we are assigning our arguments to variables
# Note the order (1,2, etc.). Like with functions in R
# the order of our inputs matters so they are asigned to 
# the proper variable 
inbam=$1
prefix=$2
logfile=$3

# Here we are defining the reference genome
reference=/project2/bios10602/final-project/reference/GRCh38_full_analysis_set_plus_decoy_hla.fasta

# Here we are writing all of our current variables to a logfile
# using a logfile will provide us with a record of what we did
echo %%%  INPUT DATA OF CURRENT RUN > $logfile 2>&1
echo %%%  inbam: $inbam >> $logfile 2>&1
echo %%%  prefix: $prefix >> $logfile 2>&1
echo %%%  logFile: $logfile >> $logfile 2>&1

##########
# STEP 1 #
##########
echo %%%[`date`]%%% Indexing the input .bam file... >> $logfile 2>&1
/project2/bios10602/Lab8/software/samtools/samtools index $inbam
echo %%%[`date`]%%% Finished indexing the .bam file... >> $logfile 2>&1
echo %%%[`date`]%%% Creating the raw vcf file... >> $logfile 2>&1
# ${prefix}.raw.vcf will be one of the output files
/project2/bios10602/Lab8/software/samtools/samtools mpileup -r chr1 -t SP -uv -f $reference $inbam | /project2/bios10602/Lab8/software/samtools/bcftools call -mv > ${prefix}.raw.vcf
echo %%%[`date`]%%% Finished creating the raw vcf file... >> $logfile 2>&1

##########
# STEP 2 #
##########
echo %%%[`date`]%%% Filtering the raw vcf file... >> $logfile 2>&1
# ${prefix}.flt.vcf will be one of the output files
/project2/bios10602/Lab8/software/samtools/vcfutils.pl varFilter -D100 ${prefix}.raw.vcf > ${prefix}.flt.vcf
echo %%%[`date`]%%% Finished filtering the vcf... >> $logfile 2>&1
