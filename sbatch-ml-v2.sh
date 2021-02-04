#!/bin/bash
#SBATCH --partition=standard -J med-sim
#    #SBATCH -c 20 --mem=50480
#SBATCH -c 11 --mem=50480
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 2-00:00
#SBATCH --array=1-3

# echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date

module load r/3.6.1/b1

./run-sims-ml-v2.sh
