#!/bin/bash

ns=(500 1000 2000 4000)
usesl=F
num_sims=1000
suffix="ml0-10-1"
weightgam=cvgam
# cores="default"

scenarios=(lll)
coefsizes=(smallalpha)

nidx=$SLURM_ARRAY_TASK_ID

scenario=${scenarios[0]}
coefsize=${coefsizes[0]}
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

