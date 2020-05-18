#!/bin/sh

for RUN in {0001..0050}
do
  sbatch deploy_sim.sh false $RUN
  sleep .1
  sbatch deploy_sim.sh true $RUN
  sleep .1
done
