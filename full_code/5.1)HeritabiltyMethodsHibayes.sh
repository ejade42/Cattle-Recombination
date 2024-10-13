#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=100:00:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=12000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-8
#SBATCH --output=slurm_output/slurm-%A_%a-HeritabilityMethodsHibayes.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: uses the various hibayes methods to analyse the simulated data

## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a


## EXTRACT METHOD TO USE
i=$SLURM_ARRAY_TASK_ID
methods="BayesA BayesB BayesBpi BayesC BayesCpi BayesL BayesR BayesRR BSLMM"
methods=($methods)
method=${methods[$i]}
echo $method


## ANALYSIS ONLY - simulation is in 5.0 script
echo \($(date +"%Y-%m-%d %T")\) Currently looking at method $method; echo; echo
cd method_comparison
echo va,n,true_h2,est_h2,est_sd,est_lower,est_upper,h2_error > iterations_$method".csv"
Rscript ../scripts/heritability_methods_analysis_hibayes.r $method
echo; echo; echo \($(date +"%Y-%m-%d %T")\) Analysis complete for method $method
echo; echo; echo; echo; echo
echo \($(date +"%Y-%m-%d %T")\) done
