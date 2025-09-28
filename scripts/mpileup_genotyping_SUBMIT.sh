#!/bin/sh
#SBATCH --job-name=genotypingGroup5
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --account=bios10007
#SBATCH --ntasks-per-node=48

#define filename variables for input and output (change cnet)
inbam=/scratch/midway3/sophiaw09/final-project/alignment/SRR701471.sorted.bam
prefix=SRR701471
logfile=${prefix}.genotype.log

#run mpileup genotyping script
/scratch/midway3/sophiaw09/final-project/scripts/mpileup_genotyping.sh $inbam $prefix $logfile
