#!/bin/bash

# Here we are assigning our arguments to variables
# These arguments are specified at the time we execute this script (with the SUBMIT.sh script), such as:
# bash Bwa_alignment_paired-end.sh reads_1.fastq reads_2.fastq my_alignment my_alignment.log

# Note the order (1,2, etc.). Like with functions in R, the order of our inputs matters so
# they are asigned to the proper variable. In this case, $1 represents the first input
# argument, $2 the second, etc.
fastq1=$1
fastq2=$2
prefix=$3
logfile=$4

# Define the reference genome to which our reads will be aligned
reference=../reference/GRCh38_full_analysis_set_plus_decoy_hla.fasta

# Write all of our current variables to a logfile
# using a logfile will provide us with a record of what we did
echo %%%  INPUT DATA OF CURRENT RUN >> $logfile 2>&1
echo %%%  fastq1: $fastq1 >> $logfile 2>&1
echo %%%  fastq2: $fastq2 >> $logfile 2>&1
echo %%%  prefix: $prefix >> $logfile 2>&1
echo %%%  logFile: $logfile >> $logfile 2>&1

#########################
# STEP 1: BWA Alignment #
#########################
echo %%%[`date`]%%% Running BWA alignment on fastq 1... >> $logfile 2>&1
# Alignment command for fastq1 below... will generate a .sai file for end 1
../software/bwa/bwa aln -q 5 -t 28 $reference $fastq1 > ${prefix}_1.sai #2>> $logfile 
echo %%%[`date`]%%% Finished BWA alignment on fastq 1! >> $logfile 2>&1

echo %%%[`date`]%%% Running BWA alignment on fastq 2... >> $logfile 2>&1
# Alignment command for fastq2 below... will generate a .sai file for end 2
../software/bwa/bwa aln -q 5 -t 28 $reference $fastq2 > ${prefix}_2.sai #2>> $logfile
echo %%%[`date`]%%% Finished BWA alignment on fastq 2! >> $logfile 2>&1

##########################
# STEP 2: BWA create bam #
##########################
echo %%%[`date`]%%% Running sampe... >> $logfile 2>&1
../software/bwa/bwa sampe -P $reference ${prefix}_1.sai ${prefix}_2.sai $fastq1 $fastq2 2>> $logfile | ../software/samtools/samtools view -b -S - > ${prefix}.bam  2>> $logfile
echo %%%[`date`]%%% Finished sampe! >> $logfile 2>&1

#################################
# STEP 3: Sort the bam file #
#################################
# Use samtools sort to sort and index our bam file
echo %%%[`date`]%%% Sorting and indexing ${prefix}.bam into ${prefix}.sorted.bam >> $logfile 2>&1
../software/samtools/samtools sort -@ 28 -m 1500M -T ${prefix} -o ${prefix}.sorted.bam ${prefix}.bam 
echo %%%[`date`]%%% Finished sorting! >> $logfile 2>&1
