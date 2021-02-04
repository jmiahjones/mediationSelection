#!/bin/bash

n=500
num_sims=1000
suffix="boot-mgcv-oracle"

scenario=rnn
coefsize=large
weightgam=cvgam
usesl=T
# scenarios=(lll lnn nnn)
# scenario=${scenarios[$SLURM_ARRAY_TASK_ID-1]}


# for scenario in lll lln lnl nll nln nnl lnn nnn; do

for usesl in T F; do
  for n in 500; do
    for scenario in lnn lll nnn; do
      for coefsize in large small; do
        # for weightgam in 0.5 1 2; do
          echo "Beginning scenario $scenario."
          Rscript --verbose ./R/var-selection-sim-cv.R $n $num_sims $scenario $coefsize \
            $weightgam $usesl $suffix \
            > ./logs/$n-$num_sims-$scenario-$coefsize-$weightgam-$usesl-$suffix.log \
            2>&1
          echo "Completed scenario $scenario."
        # done
      done
    done
  done
done
