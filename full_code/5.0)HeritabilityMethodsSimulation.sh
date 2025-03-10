#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=00:15:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=4000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=4                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_output/slurm-%j-HeritabilityMethodsSimulation.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: simulates phenotypic data with known true h^2, in order to test how good the estimation tools are


## LOAD MODULES
echo \($(date +"%Y-%m-%d %T")\) start
module load PLINK/1.09b6.16
module load R/4.2.1-gimkl-2022a


# Create input files in heritability directory
mkdir -p ./heritability
mkdir -p ./method_comparison
cp input/software/gcta-1.94.1 heritability

cd heritability
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sire_models --make-bed --out sires --maf 0.01
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sire_models --recode A --out sires_recoded --maf 0.01
plink --cow --bfile ../qc/all_cows_10 --keep ../files/sire_models --recode 12 --out sires --maf 0.01
cd ..

## SIMULATE PHENOTYPES
filter=0
k=10
echo; echo \($(date +"%Y-%m-%d %T")\) simulating phenotypes
Rscript scripts/heritability_methods_simulation.r $filter $k
echo \($(date +"%Y-%m-%d %T")\) phenotypes simulated


## CREATE GRMS
./gcta-1.94.1 --autosome-num 29 --bfile sires --maf 0.01 --make-grm --out grm_all_sires
