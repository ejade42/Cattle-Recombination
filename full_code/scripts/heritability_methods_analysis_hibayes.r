## Runs in method_comparison

start_time <- Sys.time()


## SETTING UP
##-------------------------------------------------------------------------------------------
## Load packages
library("data.table")            ## v1.14.10
library("hibayes")               ## v3.0.1
library("tidyverse")             ## v2.0.0
source("../scripts/functions.r") ## load common functions
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




## LOADING SIMULATED PHENOTYPES
##-------------------------------------------------------------------------------------------
## Load simulated phenotype data
simulated_phens <- read.table("simulated_phens")
n_sims <- 300
colnames(simulated_phens) <- c("family", "sires", paste0("phen_", 1:n_sims))

# Load real sire data
sire_models <- read.csv("../files/sire_models.csv")

# Merge together
phenotypes <- merge(sire_models, simulated_phens)
phenotypes <- phenotypes[,2:ncol(phenotypes)]
#write.csv(phenotypes, "combined_phenotypes.csv", row.names = F)
##-------------------------------------------------------------------------------------------




## LOADING AND IMPUTING GENOTYPES
##-------------------------------------------------------------------------------------------
## Load genotype and marker data (from PLINK output files)
data <- read_plink(bfile="../heritability/sires", threads=4, out=tempfile(), mode="A", impute=TRUE)
genotypes <- read.table("../heritability/sires_recoded.raw", header = T)

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

fprint("{n_vari} SNPs and {n_indi} individuals remain")
print_blank()
##-------------------------------------------------------------------------------------------




## FITTING MODELS TO SIMULATED DATA
##-------------------------------------------------------------------------------------------
## Simulation parameters (MUST MATCH power_analysis_simulations.r)
va_values <- c(seq(4, 2, -0.25), seq(1.9, 0.1, -0.1), 0.05, 0.005)
N <- 10

va_vec <- rep(va_values, each = 10)
n_vec  <- rep(1:N, times = 30)


## Read in Bayes alphabet method to use from filter
method <- commandArgs(trailingOnly = TRUE)
if (length(method) == 0) {stop("Please provide a method setting when calling the script")}
filename <- paste0("iterations_", method, ".csv")
fprint("Currently using method: {method}")
sim_hibayes_iterations <- read.csv(filename)


# Fit model to each simulated phenotype
loop_start_time <- Sys.time()
fprint("Fitting simulated phenotype models, start time: {round(loop_start_time)}")

confidence <- 0.025
for (i in 1:n_sims) {
    print_blank()
    phenotypes[,"phen"] <- phenotypes[,paste0("phen_", i)]

    ## Fit model
    fit <- 'ibrm(phen ~ sire_jersey_propn + halfsibs_total + sire_breeding_window + (1 | sire_birth_year) + (1 | sire_region), 
                data = phenotypes, M = geno, M.id = fam[,2], method = method,
                niter = 10000, nburn = 2000, thin = 1, printfreq = 500, threads = 8)'
    fit <- ibrm(phen ~ 1,                                                                      ## It doesn't make sense to fit covariates as the phenotypes are randomly simulated with no covariate effects
                data = phenotypes, M = geno, M.id = fam[,2], method = method,
                niter = 10000, nburn = 2000, thin = 1, printfreq = 500, threads = 8)

    va <- va_vec[i]
    n  <- n_vec[i]
    true_h2   <- round(va / 5, 3)
    est_h2    <- summary(fit)$VGR[2,1]
    est_sd    <- summary(fit)$VGR[2,2]
    est_ci    <- quantile(as.vector(fit$MCMCsamples$h2), c(confidence, 1-confidence))
    est_lower <- est_ci[1]
    est_upper <- est_ci[2]
    h2_error  <- est_h2 - true_h2

    
    row <- c(va, n, true_h2, est_h2, est_sd, est_lower, est_upper, h2_error)
    
    sim_hibayes_iterations[nrow(sim_hibayes_iterations)+1,] <- row
    write.csv(sim_hibayes_iterations, filename, row.names = F)


    fprint("Fitted model for simulated phenotype {i}/{n_sims} ({round(i/n_sims * 100, 2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")
}

print_end_message(start_time)
##-------------------------------------------------------------------------------------------
