#!/bin/sh
#SBATCH --job-name=alignmentGroup5
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --account=bios10007
#SBATCH --tasks-per-node=48

# define filepaths referenced by file variables
fastq1=/scratch/midway3/sophiaw09/final-project/data/SRR701471_1.fastq
fastq2=/scratch/midway3/sophiaw09/final-project/data/SRR701471_2.fastq
prefix=SRR701471
logfile=${prefix}.align.log

# point to the script with our analyses
/scratch/midway3/sophiaw09/final-project/scripts/Bwa_alignment_paired-end.sh $fastq1 $fastq2 $prefix $logfile
