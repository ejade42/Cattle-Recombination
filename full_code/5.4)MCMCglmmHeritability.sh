#!/bin/bash
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=150:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=24000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-500
#SBATCH --output=slurm_output/slurm-%A_%a-MCMCglmm.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: Estimates heritability using MCMCglmm models


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a

cd heritability
mkdir -p mcmcglmm


chunk_size=1
baseline=$(($chunk_size*$SLURM_ARRAY_TASK_ID))
for filter in $(seq $baseline $(($baseline+$chunk_size-1)))
do
    echo "-------------------------------------------------------------------------"
    echo \($(date +"%Y-%m-%d %T")\) "Working on filter:" $filter
    for breed in all holstein jersey subsample
    do
        echo \($(date +"%Y-%m-%d %T")\) "Working on breed:" $breed
        Rscript ../scripts/heritability_mcmcglmm.r $filter $breed
        echo; echo;
        
    done
    echo; echo
done


echo \($(date +"%Y-%m-%d %T")\) done
