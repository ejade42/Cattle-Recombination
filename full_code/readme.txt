Hi, here is all the code used in my analysis -Ev.
However, the input directory does not contain all of the required files, because the genotype information is commercially sensitive. The final output data, along with the scripts to fit models and generate graphs, is in the other directory within this github (i.e. in "outputs" rather than "full_code").


Assuming access to the input files, the analysis could be re-run as follows:

What you should have to start with:
- A bunch of .sh bash/slurm scripts in the main directory (plus the 0)settings file)
- An input directory containing the key initial data and software files
- A scripts directory containing a bunch of scripts, mostly R

Everything else will be generated as you run them.
In general, if a particular process produces just a few output files (e.g. the extracted dataset of recombinations or the per-sire linear model estimates), they will be saved to the files directory.
If there are many output files, they will be contained in their own directory (e.g. qc, LINKPHASE_outputs, method_comparison, heritability, genome_partitioning, gwas) - these should be generated as the scripts run.
Final figures and the GWAS results table will be available in the output directory, but some of the results are printed to the console as the final script runs, e.g. model coefficients and confints - this can be accessed from the slurm output file (or the R scripts can be run directly in terminal and the output will display).


The scripts labelled 1) up to 7) do all of the analysis here, and should be run in order to do the full analysis.
However, some can be run simultaneously.

These must be run in sequence (not simultaneously) before anything else can be done:
- 1)QC
- 2)Linkphase
- 3)ExtractRecombinations

Then, scripts labelled 4, 5, and 6 can all be run independently (i.e. in parallel).
However, the first one of each number must be run before the others can.
So, you can simultaneously run:
- 4)DredgeModels
- 5.0)HeritabilityMethodsSimulation
    - Once this is finished, you can simultaneously run 5.1, 5.2, and 5.3
    - Note that 5.3 will say it failed; that is because the "subsample" breed can't do filters above 365. As long as that got to 365 and the others finished, there is no error.
- 6.0)GenomePartitioning
    - Once this is finished, you can run 6.1
    - Once 6.1 is finished, you can run 6.2 (to use the GWAS results to find nearby genes). This is the first script to make a file in the output directory.

Finally, the output script which generates all other graphs and tables must be run last, after all others have finished.
Alternatively, the R scripts prefixed "graphs_" can be run independently to generate output (graphs/tables etc) for just that analysis.
Note that some output e.g. model estimates, P-values, and confidence intervals is visible only in the console/slurm output txt.
- 7)GraphsTablesOutput


IMPORTANT NOTE: 
Some of the scripts (2, 5.1, 5.3, and 6.1) are slurm "Array" scripts, which run in parallel.
The number of instances is specified by the #SBATCH --array= argument, which provides a range of numbers.
Each number in this range, inclusive, will generate its own instance of the script.
The instance number can be accessed from the $SLURM_ARRAY_TASK_ID variable.
To use these in a standard bash environment, run the whole script inside a loop, replacing the "i=$SLURM_ARRAY_TASK_ID" with "for i in {a..b}" iterating across the same range of values specified in the slurm array argument.


IF YOU WANT TO GENERATE THE PHENOTYPES FOR GWAS:
These are generated in line 58 of the 6.0)GenomePartitioning script.
This calls the R script "extract_hibayes_corrected_phenotypes.r", which takes a numerical argument in the command line representing filter.
This R script will filter to only animals with breeding window >= the specified argument, and generate a file of phenotypes for both mean and slope. It will save them as "recombinations_filter{filter}_{phenotype}.phen" in the directory in which the script runs.


Versions of packages are commented in each script. However, tidyverse 2.0.0 can load a range of different versions of its component packages. The core tidyverse package versions used are:
dplyr     1.1.2     
forcats   1.0.0     
ggplot2   3.4.2     
lubridate 1.9.2     
purrr     1.0.1
readr     2.1.4
stringr   1.5.0
tibble    3.2.1
tidyr     1.3.0


For the "input/gene_names_pos_forR.csv" file, the creator (JS) gave the following method:
"I downloaded the information about each gene from this webpage: 
https://www.ncbi.nlm.nih.gov/gene/?term=ARS-UCD1.2/bosTau9
Clicked on "Send to" and chose "File". 
This gave me a text file called gene_result.txt.
Then I just opened it in excel (as tab delimited) deleted the 7th-10th column  and converted to csv (saved  'gene_names_pos_forR.csv') then I could read it easily into R."
