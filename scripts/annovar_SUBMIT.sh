#!/bin/sh
#SBATCH --job-name=annovarGroup5
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --account=bios10007
#SBATCH --tasks-per-node=48

../scripts/annovar.sh ../genotyping/SRR701471.flt.vcf SRR701471
