#!/bin/bash
#SBATCH --partition=standard -J med-sim-par
#SBATCH -c 12 --mem=120G
#    #SBATCH -c 11 --mem=50480
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 5-00:00
#SBATCH --array=0-23

echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date

module load r/3.6.3/b1

cores=12 ./run-sims-ml-sbatch-innerpar.sh
