#!/bin/bash

ns=(500 4000)
usesls=(T F)
num_sims=1000
suffix="boot-mgcv-oracle"
# cores="default"

scenario=rnn
coefsize=large
weightgam=cvgam
usesl=T
# scenarios=(lll lnn nnn)
# scenario=${scenarios[$SLURM_ARRAY_TASK_ID-1]}

nidx=$(($SLURM_ARRAY_TASK_ID / 2))
useslidx=$(($SLURM_ARRAY_TASK_ID % 2))


n=${ns[$nidx]}
usesl=${usesls[$useslidx]}

echo "n=$n"
echo "usesl=$usesl"

# for scenario in lll lln lnl nll nln nnl lnn nnn; do

# for usesl in T F; do
  # for n in 500; do
    for scenario in lnn lll nnn; do
      for coefsize in large small; do
        # for weightgam in 0.5 1 2; do
          echo "Beginning scenario $scenario."
          Rscript --verbose ./R/var-selection-sim-cv.R $n $num_sims $scenario $coefsize \
            $weightgam $usesl $suffix $cores \
            > ./logs/$n-$num_sims-$scenario-$coefsize-$weightgam-$usesl-$suffix.log \
            2>&1
          echo "Completed scenario $scenario."
        # done
      done
    done
  # done
# done
