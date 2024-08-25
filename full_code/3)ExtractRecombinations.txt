#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=12:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=4000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=4                                         # Number of threads (level of paralellism)
#SBATCH --partition=large
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL

#SBATCH --output=slurm_output/slurm-%j-ExtractRecombinations.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: uses the linkphase output to create the csv dataset for modelling,
##              applies per-pair and per-chromosome filters, and
##              fits linear models to sires to estimate mean and slope phenotypes

module load R/4.2.1-gimkl-2022a
echo \($(date +"%Y-%m-%d %T")\) start


echo \($(date +"%Y-%m-%d %T")\) extracting max genetic distance per chromosome
cd scripts
chmod +x chromosome_total_distances.txt
./chromosome_total_distances.txt
cd ..
echo \($(date +"%Y-%m-%d %T")\) max genetic distance per chromosome extracted


echo \($(date +"%Y-%m-%d %T")\) extracting recombination, date, and breed information
Rscript scripts/extract_recombinations.r
echo \($(date +"%Y-%m-%d %T")\) recombination information extracted; echo


echo \($(date +"%Y-%m-%d %T")\) filtering recombinations
Rscript scripts/filter_recombinations.r
echo \($(date +"%Y-%m-%d %T")\) recombinations filtered; echo


echo \($(date +"%Y-%m-%d %T")\) fitting sire models
Rscript scripts/fit_sire_models.r
echo \($(date +"%Y-%m-%d %T")\) sire models fitted; echo


module unload R/4.2.1-gimkl-2022a
echo \($(date +"%Y-%m-%d %T")\) done