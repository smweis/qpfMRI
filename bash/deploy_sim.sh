#!/bin/bash
#SBATCH --job-name=sim_test    # Job name
#SBATCH --qos=stevenweisberg-b
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=stevenweisberg@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=sim_test_%A.log   # Standard output and error log


# See details of /code/qpfMRI/main/compiledSimulate.m for how to call this.

sleep 1

pwd; hostname; date

ml matlab

export MCR_CACHE_ROOT=$SCRATCH

d=`date +%m-%d-%Y`

cd /ufrc/stevenweisberg/stevenweisberg/qpfMRIResults/compiled.$d

./run_compiledSimulate.sh /apps/matlab/mcr/2019b/v97 \
logistic \
param1Lower .01 \
param1nDivisions 30 \
param1Upper 1.0 \
param2Lower .01 \
param2nDivisions 30 \
param2Upper 1.0 \
param3Lower .8 \
param3nDivisions 7 \
param3Upper 1.4 \
param4Lower .3 \
param4nDivisions 10 \
param4Upper 1.5 \
stimLower .01 \
stimUpper 1.0 \
stimnDivisions 30 \
noiseSDLower .1 \
noiseSDUpper .8 \
noiseSDnDivisions 7 \
qpPres $1 \ 
outNum $2