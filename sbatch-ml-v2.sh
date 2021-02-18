#!/bin/bash
#SBATCH --partition=standard -J med-sim
#SBATCH -c 24 --mem=50G
#    #SBATCH -c 11 --mem=50480
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 5-00:00
#SBATCH --array=0-1

echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date

module load r/3.6.1/b1

cores="default"
./run-sims-ml-sbatch.sh
