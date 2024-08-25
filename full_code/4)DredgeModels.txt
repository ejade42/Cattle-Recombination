#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=12:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=24000                                               # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=bigmem
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL

#SBATCH --output=slurm_output/slurm-%j-Dredge.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: fits a very large linear mixed model, then assesses all submodels (fixed effects only) and picks the best by AICc
##              then it manually checks AICc values for each combination of random effects.

echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a
echo \($(date +"%Y-%m-%d %T")\) fitting superset model, running dredge; echo
Rscript scripts/dredge_model.r
Rscript scripts/dredge_manual_followup.r
echo done