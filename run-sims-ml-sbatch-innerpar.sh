#!/bin/bash

n=4000
usesl=T
num_sims=1000
suffix="boot-mgcv-oracle"
weightgam=cvgam
# cores="default"

scenarios=(lll lnn nnn)
coefsizes=(large small)

scenidx=$(($SLURM_ARRAY_TASK_ID / 2))
coefidx=$(($SLURM_ARRAY_TASK_ID % 2))


scenario=${scenarios[$scenidx]}
coefsize=${coefsizes[$coefidx]}

echo "scenario=$scenario"
echo "coefsize=$coefsize"


# for scenario in lnn lll nnn; do
#   for coefsize in large small; do
#       echo "Beginning scenario $scenario."
      Rscript --verbose ./R/var-selection-sim-cv.R $n $num_sims $scenario $coefsize \
        $weightgam $usesl $suffix $cores \
        > ./logs/$n-$num_sims-$scenario-$coefsize-$weightgam-$usesl-$suffix.log \
        2>&1
#       echo "Completed scenario $scenario."
#   done
# done

