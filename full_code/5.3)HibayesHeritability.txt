#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=100:00:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=32000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-3
#SBATCH --output=slurm_output/slurm-%A_%a-HibayesHeritability.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: Estimates heritability for 4 groups/subsets of the sires, using hibayes


echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a


## IMPORTANT: BREED MUST BE ONE OF "all" "subsample" "holstein" "jersey" 
i=$SLURM_ARRAY_TASK_ID
breeds="all subsample holstein jersey"
breeds=($breeds)
breed=${breeds[$i]}
echo $breed


# Set up blank csv to store iteration data
cd heritability
echo sire_breeding_window,n,mean_h2,mean_sd,mean_lower,mean_upper,slope_h2,slope_sd,slope_lower,slope_upper > iterations_hibayes_$breed".csv"

for filter in {0..500}
do
    echo -----------------------------------------------------------------------------------------------
    echo FILTER CURRENTLY SET TO: $filter

    echo \($(date +"%Y-%m-%d %T")\) creating running hibayes analysis for mean, sd, and slope phenotypes
    Rscript ../scripts/hibayes_analysis.r $filter $breed
    
    echo \($(date +"%Y-%m-%d %T")\) finished for filter $filter
    echo -----------------------------------------------------------------------------------------------
    echo; echo; echo; echo; echo; echo
done

echo \($(date +"%Y-%m-%d %T")\) done