#!/bin/bash
#SBATCH --partition=standard -J med-allres
#SBATCH -c 1 --mem=20G
#    #SBATCH -c 11 --mem=50480
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH -t 5-00:00

hostname
date

module load r/3.6.1/b1

Rscript --verbose ./R/create-allres.R > ./logs/result/log-allres.out
exit 0
