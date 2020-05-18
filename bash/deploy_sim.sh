#!/bin/bash
#SBATCH --job-name=sim_test    # Job name
#SBATCH --qos=stevenweisberg-b
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=stevenweisberg@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=sim_test_%a.log   # Standard output and error log


# See details of /code/qpfMRI/main/compiledSimulate.m for how to call this.

sleep 1

pwd; hostname; date

ml matlab

export MCR_CACHE_ROOT=$SCRATCH

d=`date +%m-%d-%Y`

cd /ufrc/stevenweisberg/stevenweisberg/qpfMRIResults/compiled.$d

./run_compiledSimulate.sh /apps/matlab/mcr/2019b/v97 \
doeTemporalModel \
param1Lower .899 \
param1Interval .025 \
param1Upper 1.099 \
param2Lower .01 \
param2Interval .04 \
param2Upper .4 \
param3Lower .01 \
param3Interval .04 \
param3Upper .4 \
param4Lower .8 \
param4Interval .1 \
param4Upper 1.4 \
param5Lower .3 \
param5Interval .2 \
param5Upper 1.0 \
param1Simulated 1.004 \
param2Simulated .016 \
param3Simulated .118 \
param4Simulated 1.0 \
param5Simulated .1 \ 
qpPres $1 \ 
outNum $2