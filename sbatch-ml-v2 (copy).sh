#!/bin/bash
#SBATCH --partition=debug -J test

#    echo Array index: $SLURM_ARRAY_TASK_ID
hostname
date
echo $PWD
module load r/3.6.1/b1


n=200
usesls=T
num_sims=2
suffix="foo"
weightgam=cvgam
cores=1


scenario=lll
coefsize=large

echo "scenario=$scenario"
echo "coefsize=$coefsize"


      Rscript --verbose ./test.R $n $num_sims $scenario $coefsize \
        $weightgam $usesl $suffix $cores

echo "Complete!"

# Rscript --verbose ./test.R

exit 0
