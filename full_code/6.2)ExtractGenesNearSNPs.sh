#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=01:00:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=4000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)

#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_output/slurm-%j-ExtractGenesNearSNPs.out

## AUTHOR:  Evelyn Jade and Anna Santure (R script)
## DESCRIPTION: Takes GWAS output and list of genes to find genes near low-p-value SNPs (p < 10^-4).



## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load R/4.2.1-gimkl-2022a

mkdir -p ./output


## Create list files for R
echo \($(date +"%Y-%m-%d %T")\) creating list files for R
cd files
echo GWAS_mean_filter_0_breed_all.mlma > gwas_file_all.list
echo GWAS_mean_filter_365_breed_all.mlma >> gwas_file_all.list
echo GWAS_slope_filter_0_breed_all.mlma >> gwas_file_all.list
echo GWAS_slope_filter_365_breed_all.mlma >> gwas_file_all.list

echo GWAS_mean_filter_0_breed_holstein.mlma > gwas_file_holstein.list
echo GWAS_mean_filter_365_breed_holstein.mlma >> gwas_file_holstein.list
echo GWAS_slope_filter_0_breed_holstein.mlma >> gwas_file_holstein.list
echo GWAS_slope_filter_365_breed_holstein.mlma >> gwas_file_holstein.list

echo GWAS_mean_filter_0_breed_jersey.mlma > gwas_file_jersey.list
echo GWAS_mean_filter_365_breed_jersey.mlma >> gwas_file_jersey.list
echo GWAS_slope_filter_0_breed_jersey.mlma >> gwas_file_jersey.list
echo GWAS_slope_filter_365_breed_jersey.mlma >> gwas_file_jersey.list
cd ..


## Run R script
echo; echo \($(date +"%Y-%m-%d %T")\) running R script
Rscript scripts/genes_near_SNPs.r


echo; echo \($(date +"%Y-%m-%d %T")\) done
