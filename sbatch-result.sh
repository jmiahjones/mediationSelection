#!/bin/bash
#SBATCH --partition=standard -J med-result
#SBATCH -c 12 --mem=50G
#    #SBATCH -c 11 --mem=50480
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 5-00:00
#SBATCH --array=0-11

echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date

module load r/3.6.1/b1

totalarray=12
cores=12

Rscript --verbose ./R/create-results.R $((SLURM_ARRAY_TASK_ID + 1)) $totalarray $cores > ./logs/result/log-$SLURM_ARRAY_TASK_ID.out

