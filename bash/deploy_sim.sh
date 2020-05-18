#!/bin/bash
#SBATCH --job-name=sim_test    # Job name
#SBATCH --qos=stevenweisberg-b
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=stevenweisberg@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=sim_test_%a.log   # Standard output and error log


# Call this script with the following positional arguments:
# 1. Whether to use Q+ or not.
# 2. Name of the individual output file ('doe_[name].csv')

sleep 1

pwd; hostname; date

ml matlab

export MCR_CACHE_ROOT=$SCRATCH



cd /ufrc/stevenweisberg/stevenweisberg/compiledDoe/

./run_doeSimulate.sh /apps/matlab/mcr/2019b/v97 $4 $5 $6 1 .4 800 12 $1 $2 $3

date
