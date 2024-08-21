## Runs in "genome partitioning", copied to "gwas", to extract the hibayes model's residuals in order to perform variance partitioning and GWAS on them

start_time <- Sys.time()

## SETTING UP
##-------------------------------------------------------------------------------------------
## Load packages
library("data.table")
library("hibayes")
library("tidyverse")
source("../scripts/functions.r")    ## load common functions
options(bigmemory.typecast.warning=FALSE)

print_blank()
fprint("Start time: {round(start_time)}")


## Define function to replace NA values with the column mean
impute_na_by_col <- function(dataframe) {
    n_cols <- ncol(dataframe)
    for (i in 1:n_cols) {
        na <- is.na(dataframe[,i])
        dataframe[na,i] <- mean(dataframe[,i], na.rm = T)
    }
    dataframe
}
##-------------------------------------------------------------------------------------------




## LOADING PHENOTYPES
##-------------------------------------------------------------------------------------------
## Load phenotype data and create pheno object
filter <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(filter) == 0) {stop("Please provide a filter setting when calling the script")}
fprint("Filter: {filter} days.")
sire_models <- read.csv("../files/sire_models.csv")
fprint("{nrow(sire_models)} sires present before filter.")
sire_models <- sire_models[sire_models$sire_breeding_window >= filter,]
fprint("{nrow(sire_models)} sires remain after filter.")


sire_models <- sire_models[,c("sires", "mean", "mean_sigma", "lm_slope", "sire_jersey_propn", "sire_birth_year", "sire_region", "halfsibs_total", "sire_breeding_window")]
colnames(sire_models)[2:4] <- c("mean", "sd", "slope")
##-------------------------------------------------------------------------------------------




## LOADING AND IMPUTING GENOTYPES
##-------------------------------------------------------------------------------------------
## Load genotype and marker data (from PLINK output files)
data <- read_plink(bfile="sires", threads=4, out=tempfile(), mode="A", impute=TRUE)
genotypes <- read.table("sires_recoded.raw", header = T)

fam  <- data$fam
map  <- data$map
geno <- genotypes[,7:ncol(genotypes)]


## NORMAL CODE FROM HERE
map[1:10,1:ncol(map)]
geno[1:10,1:5]
fprint("{sum(is.na(geno))} NAs in data")

## Impute NAs based on the column average                    
## Seems like I probably have to round otherwise I will get a non-integer genotype?
fprint("Imputing missing genotypes by filling in the column/SNP average genotype (rounded), start time: {round(Sys.time())}")
geno <- impute_na_by_col(geno)
warnings()
fprint("Genotypes imputed, end time: {round(Sys.time())}")
fprint("{sum(is.na(geno))} NAs in data")


## Remove any columns (SNPs) with NAs (missing genotypes)
fprint("Removing columns/SNPs with missing (NA) genotypes")
geno <- geno %>% select(where(~ !any(is.na(.))))
fprint("{sum(is.na(geno))} NAs in data")

n_vari <- ncol(geno)
n_indi <- nrow(geno)

fprint("{n_vari} SNPs and {n_indi} individuals remain after NA filter (independent of breeding window filter)")
print_blank()

##-------------------------------------------------------------------------------------------




## FITTING MODELS
##-------------------------------------------------------------------------------------------
phenotypes <- c("mean", "slope")
row <- c(filter, nrow(sire_models))
confidence <- 0.025


for (phenotype in phenotypes) {
    recombinations.phen <- read.csv("../files/sire_models.csv")[,c("family", "sires")]
    rownames(recombinations.phen) <- recombinations.phen$sires
    
    sire_models[,"phen"] <- sire_models[,phenotype]
    print_blank()
    fprint("Fitting model for phenotype {phenotype}")

    ## Fit model
    fit <- ibrm(phen ~ sire_jersey_propn + halfsibs_total + sire_breeding_window + (1 | sire_birth_year) + (1 | sire_region), 
                data = sire_models, M = geno, M.id = fam[,2], method = "BayesCpi",
                niter = 10000, nburn = 2000, thin = 1, printfreq = 500, threads = 8)
    quick <- 'fit <- ibrm(phen ~ sire_jersey_propn + halfsibs_total + sire_breeding_window + (1 | sire_birth_year) + (1 | sire_region), 
                data = sire_models, M = geno, M.id = fam[,2], method = "BayesCpi",
                niter = 100, nburn = 20, thin = 1, printfreq = 500, threads = 8)'

    alpha <- fit$alpha
    gen_effects <- as.matrix(geno) %*% alpha
    names(gen_effects) <- fam[,2]
    str(gen_effects)
    
    residuals <- fit$e$e
    names(residuals) <- remaining_individuals <- fit$e$id
    str(residuals)
    
    unexplained_variance <- gen_effects[remaining_individuals] + residuals[remaining_individuals]
    str(unexplained_variance)

    recombinations.phen[remaining_individuals,"phenotype"] <- unexplained_variance
    recombinations.phen <- drop_na(recombinations.phen) 
    
    write.csv(recombinations.phen, fstring("recombinations_filter{filter}_{phenotype}.csv"), row.names = F)
    write.table(recombinations.phen, fstring("recombinations_filter{filter}_{phenotype}.phen"), row.names = F, col.names = F)
}
##-------------------------------------------------------------------------------------------




## SAVING OUTPUT
##-------------------------------------------------------------------------------------------

print_end_message(start_time)
##-------------------------------------------------------------------------------------------
