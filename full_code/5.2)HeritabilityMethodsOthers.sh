#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=48:00:00                                            # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=16000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=8                                         # Number of threads (level of paralellism)
#SBATCH --partition=milan
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL

#SBATCH --output=slurm_output/slurm-%j-HeritabilityMethodsOthers.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: Uses the non-hibayes methods to analyse the simulated data


## LOAD MODULES
## ----------------------------------------------------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) start
module load PLINK/1.09b6.16
module load R/4.2.1-gimkl-2022a
## ----------------------------------------------------------------------------------------------------------------------



## CREATE INPUT FILES
## ----------------------------------------------------------------------------------------------------------------------
## Copy software and give exe permissons
cp input/software/gcta-1.94.1 method_comparison
cp input/software/gctb-2.05beta method_comparison
cp input/software/ldak5.2.linux method_comparison

cd method_comparison
chmod +x gcta-1.94.1
chmod +x gctb-2.05beta
chmod +x ldak5.2.linux


## Create blank results CSVs
echo; echo \($(date +"%Y-%m-%d %T")\) creating input files
echo va,n,true_h2,est_h2,p,h2_error > iterations_gcta.csv
echo va,n,true_h2,est_h2,h2_error > iterations_gctb.csv
echo va,n,true_h2,est_h2,sd,h2_error > iterations_ldak.csv


## Create PLINK/GRM input files (for GCTA/B and LDAK)
./gcta-1.94.1 --autosome-num 29 --bfile ../heritability/sires --make-grm --out sires_grm
echo \($(date +"%Y-%m-%d %T")\) input files created; echo; echo; echo; echo; echo
## ----------------------------------------------------------------------------------------------------------------------




## RUN LOOP TO ANALYSE EACH SIMULATED PHENOTYPE
## ----------------------------------------------------------------------------------------------------------------------
echo \($(date +"%Y-%m-%d %T")\) analysing simulated phenotype data \(GCTA, GCTB, LDAK\)
for i in {1..300}
do
    echo -----------------------------------------------------------------------------------------------
    echo Currently looking at phenotype $i"/300" \(10 phenotypes simulated for each of 30 h2 values\)
    
    # RUN GCTA, GCTB, LDAK ANALYSIS
    ./gcta-1.94.1 --autosome-num 29 --grm sires_grm --pheno simulated_phens --mpheno $i --reml --out sires_results
    ./gctb-2.05beta --bfile ../heritability/sires --pheno simulated_phens --mpheno $i --out sires_results
    ./ldak5.2.linux --grm sires_grm --pheno simulated_phens --mpheno $i --reml sires_results --kinship-details NO --constrain YES --max-threads 8
    
    # EXTRACT RESULTS TO CSV
    Rscript ../scripts/heritability_methods_iterate_gcta.r $i
    Rscript ../scripts/heritability_methods_iterate_gctb.r $i
    Rscript ../scripts/heritability_methods_iterate_ldak.r $i

    echo \($(date +"%Y-%m-%d %T")\) finished for phenotype $i
    echo -----------------------------------------------------------------------------------------------
    echo; echo; echo; echo; echo; echo
done
## ----------------------------------------------------------------------------------------------------------------------


echo; echo; echo; echo; echo; echo \($(date +"%Y-%m-%d %T")\) done
    