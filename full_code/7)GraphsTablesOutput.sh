#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=01:00:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_output/slurm-%j-GraphsTablesOutput.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: creates the final output figures etc


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a

mkdir -p ./output

echo --------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) WORKING ON: GLOBAL AGE EFFECT; echo
Rscript scripts/graphs_global_age_effect.r
echo --------------------------------------------------------------------------
echo; echo; echo; echo; echo


echo --------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) WORKING ON: HERITABILITY METHODS; echo
Rscript scripts/graphs_heritability_methods.r
echo --------------------------------------------------------------------------
echo; echo; echo; echo; echo


echo --------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) WORKING ON: HERITABILITY ESTIMATION; echo
Rscript scripts/graphs_heritability.r
Rscript scripts/graphs_heritability_mcmcglmm.r
echo --------------------------------------------------------------------------
echo; echo; echo; echo; echo


echo --------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) WORKING ON: GENOME PARTITIONING; echo
Rscript scripts/graphs_genome_partitioning.r
echo --------------------------------------------------------------------------
echo; echo; echo; echo; echo


echo --------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) WORKING ON: GWAS; echo
Rscript scripts/graphs_gwas.r
echo --------------------------------------------------------------------------
echo; echo;


echo \($(date +"%Y-%m-%d %T")\) done
