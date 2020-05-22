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

d='05-19-2020'

cd /ufrc/stevenweisberg/stevenweisberg/qpfMRIResults/compiled.$d

./run_compiledSimulate.sh /apps/matlab/mcr/2019b/v97 \
logistic \
param1Lower .01 \
param1nDivisions 20 \
param1Upper 1.0 \
param2Lower .01 \
param2nDivisions 20 \
param2Upper 1.0 \
param3Lower .8 \
param3nDivisions 7 \
param3Upper 1.4 \
param4Lower .3 \
param4nDivisions 10 \
param4Upper 1.5 \
param1Spacing lin \
param2Spacing lin \
param3Spacing lin \
param4Spacing lin \
stimDomainLower .01 \
stimDomainUpper 1.0 \
stimDomainnDivisions 30 \
qpPres $1 \ 
outNum $2