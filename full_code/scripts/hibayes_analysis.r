## Runs in "heritability", as part of the script to iterate over different breeding window filters, for one breed only

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
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2) {stop("Please provide a filter and breed setting when calling the script")}
filter <- as.numeric(arguments[1])
breed <- as.character(arguments[2])

fprint("Filter: {filter} days.")
fprint("Breed: {breed} only.")


sire_models <- read.csv("../files/sire_models.csv")
fprint("{nrow(sire_models)} sires present before filter.")
sire_models <- sire_models[sire_models$sire_breeding_window >= filter,]
fprint("{nrow(sire_models)} sires remain after filter.")

if (tolower(breed) == "all") {
    
} else if (tolower(breed) == "subsample") {
    sample_size <- nrow(filter(sire_models, sire_breeding_window >= 365))
    cows_to_keep <- sample(1:nrow(sire_models), sample_size, replace = F)
    sire_models <- sire_models[cows_to_keep,]
    fprint("Sires randomly sampled (without replacement) down to {nrow(sire_models)} sires, the number left from all sires at filter = 365.")
    if (filter > 365) {
        stop("Subsample only runs up to filter = 365. Ending process.")
    }
} else if (tolower(breed) == "holstein") {
    sire_models <- sire_models[sire_models$sire_jersey_propn <= 1/16,]
} else if (tolower(breed) == "jersey") {
    sire_models <- sire_models[sire_models$sire_jersey_propn >= 15/16,]

} else {stop('Please provide "all", "subsample", "jersey", or "holstein" as the breed')}
fprint('{nrow(sire_models)} sires remain of the breed "{breed}".')


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
    sire_models[,"phen"] <- sire_models[,phenotype]
    print_blank()
    fprint("Fitting model for phenotype {phenotype}")

    ## Fit model
    fit <- ibrm(phen ~ sire_jersey_propn + halfsibs_total + sire_breeding_window + (1 | sire_birth_year) + (1 | sire_region), 
                data = sire_models, M = geno, M.id = fam[,2], method = "BayesCpi",
                niter = 10000, nburn = 2000, thin = 1, printfreq = 500, threads = 8)

    h2_est   <- summary(fit)$VGR[2,1]
    h2_sd    <- summary(fit)$VGR[2,2]
    h2_ci    <- quantile(as.vector(fit$MCMCsamples$h2), c(confidence, 1-confidence))
    h2_lower <- h2_ci[1]
    h2_upper <- h2_ci[2]
    
    row <- c(row, h2_est, h2_sd, h2_lower, h2_upper)
}
##-------------------------------------------------------------------------------------------




## SAVING OUTPUT
##-------------------------------------------------------------------------------------------
fprint("Updating iterations_hibayes file")
iterations_hibayes <- read.csv(fstring("iterations_hibayes_{breed}.csv"))
iterations_hibayes[nrow(iterations_hibayes)+1,] <- row
write.csv(iterations_hibayes, fstring("iterations_hibayes_{breed}.csv"), row.names = F)
print_end_message(start_time)
##-------------------------------------------------------------------------------------------
