#!/bin/bash -e
#SBATCH -A uoa03818                                               # Project Account
#SBATCH --time=02:30:00                                           # Max runtime of the job (Hours:Minutes:Seconds)
#SBATCH --mem=2000                                                # Memory/CPU (in MB)
#SBATCH --cpus-per-task=2                                         # Number of threads (level of paralellism)
#SBATCH --partition=large
#SBATCH --mail-user=email@domain.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-29
#SBATCH --output=slurm_output/slurm-%A_%a-Linkphase.out

## AUTHOR:  Evelyn Jade
## DESCRIPTION: Uses LINKPHASE3 to phase the pedigree and count the recombinations (runs in parallel, with i=1 to i=29)

module load PLINK/1.09b6.16
cd qc
i=$SLURM_ARRAY_TASK_ID
echo \($(date +"%Y-%m-%d %T")\) start; echo

## Create directories called "chr_i" for each chromosome, containing data filtered to that chromosome
echo \($(date +"%Y-%m-%d %T")\) creating chromosome $i directory
folder_name="chr_"$i
mkdir -p ./$folder_name
echo \($(date +"%Y-%m-%d %T")\) chromosome $i directory created; echo
echo \($(date +"%Y-%m-%d %T")\) extracting data
plink --cow --bfile all_cows_10 --recode 12 --chr $i --out ./$folder_name/$folder_name"_01"
echo \($(date +"%Y-%m-%d %T")\) data extracted; echo

## Create the 3 Linkphase input files inside the folder for the each chromosome
echo \($(date +"%Y-%m-%d %T")\) creating linkphase inputs
cd $folder_name
cp ../../input/newIDpedF0_noheader $folder_name".pedigree"
awk '{for (i=1;i<=NF;i++) if ((i==2 || 7<=i)) printf("%s ", $i); print ""}' $folder_name"_01.ped" > $folder_name".genotype"
awk '{for (i=1;i<=NF;i++) if ((i!=4)) printf("%s ", $i); print ""}' $folder_name"_01.map" > $folder_name".markers"
echo \($(date +"%Y-%m-%d %T")\) inputs created; echo

## Copy "LINKPHASE3" program
echo \($(date +"%Y-%m-%d %T")\) copying linkphase program
cd ../..
cp ./input/software/LINKPHASE3 ./qc/$folder_name
echo \($(date +"%Y-%m-%d %T")\) linkphase copied; echo

## Create "linkin.txt" settings file
echo \($(date +"%Y-%m-%d %T")\) creating linkin.txt input file
cd qc/$folder_name
echo "#PEDIGREE_FILE" > linkin.txt
echo $folder_name".pedigree" >> linkin.txt
echo "#GENOTYPE_FILE" >> linkin.txt
echo $folder_name".genotype" >> linkin.txt
echo "#MARKER_FILE" >> linkin.txt
echo $folder_name".markers" >> linkin.txt
echo "#HALFSIB_PHASING" >> linkin.txt
echo "yes" >> linkin.txt
echo "#HMM_PHASING" >> linkin.txt
echo "yes" >> linkin.txt
echo "#N_TEMPLATES" >> linkin.txt
echo "50" >> linkin.txt
echo "#CHECK_PREPHASING" >> linkin.txt
echo "yes" >> linkin.txt
echo \($(date +"%Y-%m-%d %T")\) linkin.txt created; echo

module unload PLINK/1.09b6.16
echo \($(date +"%Y-%m-%d %T")\) done with creating linkphase input; echo; echo; echo


module load R/4.2.1-gimkl-2022a
cd ..
echo \($(date +"%Y-%m-%d %T")\) start executing linkphase; echo

## Run Linkphase
folder_name="chr_"$SLURM_ARRAY_TASK_ID
cd $folder_name
chmod +x LINKPHASE3
echo \($(date +"%Y-%m-%d %T")\) executing linkphase
./LINKPHASE3
echo \($(date +"%Y-%m-%d %T")\) linkphase executed; echo
    
## Take output to new directory
echo \($(date +"%Y-%m-%d %T")\) copying nrec file
mkdir -p ../../LINKPHASE_outputs
cp nrec_hmm.txt ../../LINKPHASE_outputs/$folder_name"_nrec_hmm.txt"
echo \($(date +"%Y-%m-%d %T")\) nrec file copied to LINKPHASE_outputs; echo

## Create input file for extract_recombination_distances.r script and run the script
## Then copy the output csv to the LINKPHASE_outputs directory
echo \($(date +"%Y-%m-%d %T")\) extracting recombination distances
Rscript ../../scripts/extract_recombination_distances.r $SLURM_ARRAY_TASK_ID
cp $folder_name"_recombination_distances.csv" ../../LINKPHASE_outputs/
echo \($(date +"%Y-%m-%d %T")\) recombination distances extracted; echo

module unload R/4.2.1-gimkl-2022a
echo \($(date +"%Y-%m-%d %T")\) done