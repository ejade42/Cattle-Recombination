Hi, we cannot publish all of the cattle genotype data as it is commercially sensitive.
The full code for the analysis is available in the other directory within this github (i.e. in "full_code" rather than "outputs"), so the analysis and filtering steps can be verified if necessary.

This folder contains a subset of processed, output data that can be made public, as well as all of the scripts that take this data and produce the graphs (as well as a few summary tables).

The create_outputs.sh script can be run on slurm or through the terminal, and will call the R scripts sequentially to produce the outputs (which are saved to the output directory). Note that while these scripts perform exactly the same analysis and generate the same graphs as the scripts in the main analysis/code-only folder, the file structure is slightly modified (e.g. for simplicity, all files used are condensed into the "input" directory in this folder, whereas they have a variety of locations in the main analysis folder). 

Therefore, the "input" directory here does not match the "input" directory in the other folder (though the "scripts" directory here is a subset of the "scripts" directory in the main analysis folder).


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
