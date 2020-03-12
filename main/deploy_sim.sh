#!/bin/bash
#SBATCH --job-name=sim_test    # Job name
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=email@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --output=sim_test_%a.log   # Standard output and error log
#SBATCH --array=1-2
pwd; hostname; date

module load matlab

./run_doeSimulate.sh /apps/matlab/r2019b 1 .015 .75 1 0 800 12 true 1

date
