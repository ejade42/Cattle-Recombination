#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=0:30:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=2000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-2
#SBATCH --output=slurm_output/slurm-%A_%a-GWAS.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: performs a per-breed GWAS on the mean and slope phenotypes


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load PLINK/1.09b6.16
module load R/4.2.1-gimkl-2022a


## EXTRACT METHOD TO USE
i=$SLURM_ARRAY_TASK_ID
breeds="all holstein jersey"
breeds=($breeds)
breed=${breeds[$i]}
echo $breed


# Create input files in gwas directory
cd gwas
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sires_to_keep_$breed --make-bed --out sires_$breed --maf 0.01
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sires_to_keep_$breed --recode A --out sires_$breed"_recoded" --maf 0.01
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sires_to_keep_$breed --recode 12 --out sires_$breed --maf 0.01


for filter in 0 365
do
## MEAN
    ./gcta-1.94.1 --autosome-num 29 --bfile sires_$breed --keep recombinations_filter$filter"_mean.phen" --pheno recombinations_filter$filter"_mean.phen" --maf 0.01 --make-grm --out grm_mean_$breed"_"$filter 
    ./gcta-1.94.1 --autosome-num 29 --bfile sires_$breed --keep recombinations_filter$filter"_mean.phen" --grm grm_mean_$breed"_"$filter --pheno recombinations_filter$filter"_mean.phen" --maf 0.01 --mlma --out GWAS_mean_filter_$filter"_breed_"$breed --thread-num 8

    ## SLOPE
    ./gcta-1.94.1 --autosome-num 29 --bfile sires_$breed --keep recombinations_filter$filter"_slope.phen" --pheno recombinations_filter$filter"_slope.phen" --maf 0.01 --make-grm --out grm_slope_$breed"_"$filter
    ./gcta-1.94.1 --autosome-num 29 --bfile sires_$breed --keep recombinations_filter$filter"_slope.phen" --grm grm_slope_$breed"_"$filter --pheno recombinations_filter$filter"_slope.phen" --maf 0.01 --mlma --out GWAS_slope_filter_$filter"_breed_"$breed --thread-num 8
done

echo \($(date +"%Y-%m-%d %T")\) done
    