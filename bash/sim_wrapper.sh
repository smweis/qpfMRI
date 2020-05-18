#!/bin/sh

for RUN in {0001..0050}
do
  sbatch deploy_sim.sh 0 "yes_maxBOLD_fix_false_"$RUN $RANDOM 1.0 0.1 0.2
  sleep .1
  sbatch deploy_sim.sh 1 "yes_maxBOLD_fix_true_"$RUN $RANDOM 1.0 0.1 0.2
  sleep .1
done
