#!/bin/bash

ns=(500 1000 2000 4000)
usesl=T
num_sims=1000
suffix="ml0-9-1"
weightgam=cvgam
# cores="default"

scenarios=(lll lnn nnn)
coefsizes=(large small)

# go coefsize -> scenario -> n
coefidx=$(($SLURM_ARRAY_TASK_ID / 12))
scenidx=$((($SLURM_ARRAY_TASK_ID % 12) / 4))
nidx=$((($SLURM_ARRAY_TASK_ID % 12) % 4))

scenario=${scenarios[$scenidx]}
coefsize=${coefsizes[$coefidx]}
n=${ns[$nidx]}
echo "scenario=$scenario"
echo "coefsize=$coefsize"
echo "n=$n"
echo "cores=$cores"


# for scenario in lnn lll nnn; do
#   for coefsize in large small; do
#       echo "Beginning scenario $scenario."
      Rscript --verbose ./R/main.R $n $num_sims $scenario $coefsize \
        $weightgam $usesl $suffix $cores \
        > ./logs/$n-$num_sims-$scenario-$coefsize-$weightgam-$usesl-$suffix.log \
        2>&1
#       echo "Completed scenario $scenario."
#   done
# done

