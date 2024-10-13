#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=01:00:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=1500                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=1                                         # Number of threads (level of paralellism)
#SBATCH --partition=large
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL

#SBATCH --output=slurm_output/slurm-%j-QC.out

## AUTHOR:  Evelyn Jade and Anna Santure
## DESCRIPTION: takes original cattle data, updates IDs and genetic distances, and performs quality control

module load R/4.2.1-gimkl-2022a
module load PLINK/1.09b6.16
mkdir -p ./qc
mkdir -p ./files
cd qc
echo \($(date +"%Y-%m-%d %T")\) start; echo

## Extract data to bed - comment out once done (takes ages)
echo \($(date +"%Y-%m-%d %T")\) extracting data 
#plink --cow --vcf ../input/allgens_cows_recoded.vcf.gz --make-bed --out all_cows_00
echo \($(date +"%Y-%m-%d %T")\) data extracted; echo

## Update genetic distances; remove variants not in Shen linkage map and variants with missing ID (i.e. "."):
echo \($(date +"%Y-%m-%d %T")\) recoding
plink --cow --bfile all_cows_00 --recode 12 --out all_cows_01_map
echo \($(date +"%Y-%m-%d %T")\) recoded; echo
echo \($(date +"%Y-%m-%d %T")\) updating genetic distances
Rscript ../scripts/update_genetic_distances.r
echo \($(date +"%Y-%m-%d %T")\) genetic distances updated; echo
echo \($(date +"%Y-%m-%d %T")\) recoding
plink --cow --bfile all_cows_00 --recode 12 --extract all_cows_01_snp_ids_to_keep --out all_cows_01
echo \($(date +"%Y-%m-%d %T")\) recoded
cp all_cows_01.map all_cows_01_intermediate.map
cp all_cows_01_edited.map all_cows_01.map
plink --cow --file all_cows_01 --make-bed --silent --out all_cows_01
echo \($(date +"%Y-%m-%d %T")\) genetic distances merged; echo; echo

## Updating IDs:
echo \($(date +"%Y-%m-%d %T")\) updating IDs
awk '$3 = 1' ../input/newID4plink > ../files/newID4plink_01
awk '$1 = 1' ../input/update_parents.txt > ../files/update_parents_01.txt
plink --cow --bfile all_cows_01 --update-ids ../files/newID4plink_01 --make-bed --out all_cows_02
plink --cow --bfile all_cows_02 --update-parents ../files/update_parents_01.txt --make-bed --out all_cows_02
plink --cow --bfile all_cows_02 --make-bed --out all_cows_03
echo \($(date +"%Y-%m-%d %T")\) parents and IDs updated; echo

## Update sex IDs:
echo \($(date +"%Y-%m-%d %T")\) updating sex IDs
Rscript ../scripts/update_sex_ids.r all_cows_03.fam
echo \($(date +"%Y-%m-%d %T")\) sex IDs updated; echo; echo


## QC steps:
echo \($(date +"%Y-%m-%d %T")\) performing quality control
plink --cow --bfile all_cows_03 --make-bed --mind 0.01 --out all_cows_04
plink --cow --bfile all_cows_04 --make-bed --geno 0.01 --out all_cows_05
plink --cow --bfile all_cows_05 --make-bed --maf 0.01  --out all_cows_06
##plink --cow --bfile all_cows_06 --make-bed --hwe 0.001 --out all_cows_07
plink --cow --bfile all_cows_06 --make-bed --out all_cows_07
echo \($(date +"%Y-%m-%d %T")\) QC steps completed; echo; echo


## Mendel errors:
echo \($(date +"%Y-%m-%d %T")\) calculating mendel errors
plink --cow --bfile all_cows_07 --make-bed --mendel-duos --mendel --out all_cows_mendelerrors_before
plink --cow --bfile all_cows_07 --make-bed --mendel-duos --me 1 0.01 --out all_cows_08
plink --cow --bfile all_cows_08 --make-bed --mendel-duos --me 0.01 1 --out all_cows_09_after_mendel
echo \($(date +"%Y-%m-%d %T")\) mendel errors calculated; echo
echo \($(date +"%Y-%m-%d %T")\) keeping individuals that passed Mendel check OR are one of our 3151 bulls
Rscript ../scripts/mendel_filter_plus_bulls.r
echo \($(date +"%Y-%m-%d %T")\) cattle to keep identified; echo
echo \($(date +"%Y-%m-%d %T")\) extracting final individuals
plink --cow --bfile all_cows_08 --make-bed --keep all_cows_09_cattle_ids_to_keep --out all_cows_09
plink --cow --bfile all_cows_09 --make-bed --mendel-duos --mendel --out all_cows_mendelerrors_after
plink --cow --bfile all_cows_09 --make-bed --mendel-duos --set-me-missing --out all_cows_10
echo \($(date +"%Y-%m-%d %T")\) mendel errors completed; echo; echo


## Extract to text file
echo \($(date +"%Y-%m-%d %T")\) extracting to text file
#plink --cow --bfile all_cows_10 --recode 12 --out all_cows_10
echo \($(date +"%Y-%m-%d %T")\) extracted


module unload PLINK/1.09b6.16
module unload R/4.2.1-gimkl-2022a
echo \($(date +"%Y-%m-%d %T")\) done