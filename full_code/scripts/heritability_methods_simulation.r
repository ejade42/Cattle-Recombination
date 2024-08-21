## Runs in main directory

start_time <- Sys.time()


## SETTING UP
##-------------------------------------------------------------------------------------------
## Load packages
library("data.table")            ## v1.14.10
library("hibayes")               ## v3.0.1
library("tidyverse")             ## v2.0.0
source("scripts/functions.r")    ## load common functions
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




## LOADING AND IMPUTING GENOTYPES
##-------------------------------------------------------------------------------------------
## Read arguments for filter and k settings
arguments <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(arguments) < 2) {stop("Please provide a filter and k setting when calling the script")}

filter <- arguments[1]
k <- arguments[2]
fprint("Filter: {filter} days")
fprint("Causative SNPs: k = {k}")


## Load sire models data
sire_models   <- read.csv("files/sire_models.csv")
original_size <- nrow(sire_models)
sire_models   <- filter(sire_models, sire_breeding_window >= filter)
fprint("With filter set to {filter} days, {nrow(sire_models)} sires remain out of {original_size} originally present.")


## Load genotype and marker data (from PLINK output files)
data <- read_plink(bfile="heritability/sires", threads=4, out=tempfile(), mode="A", impute=TRUE)
genotypes <- read.table("heritability/sires_recoded.raw", header = T)

fam  <- data$fam
fam  <- fam[fam[,2] %in% sire_models$sires,]
map  <- data$map
geno <- genotypes[genotypes[,2] %in% sire_models$sires, 7:ncol(genotypes)]

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




## SIMULATING PHENOTYPES
##-------------------------------------------------------------------------------------------
fprint("Simulating phenotypes")

## Define simulation function
do_hibayes_sim <- function(h2, k) {  ## k = Causative SNP number (try different powers of 10)
    
    X = geno
    n = n_indi
    m = dim(geno)[2]
    qtl = sort(sample(1:m, k))

    betaTrue = array(0,m)
    betaTrue[qtl] = rnorm(k)

    X = as.matrix(X)

    g = X%*%betaTrue


    vg = var(g)
    ve = (1-h2)/h2 * vg

    y = g + rnorm(n,0,sqrt(ve))
    y
}


## Create matrix to store simulated data
all_phens <- data.frame(family = rep(1, nrow(fam)), sire = fam[,2])

## MUST MATCH OTHER SCRIPTS
va_values <- c(seq(4, 2, -0.25), seq(1.9, 0.1, -0.1), 0.05, 0.005)    # Heritability values to test - "true h^2" is va/5
N <- 10                                                               # Number of simulations per va value (= K in MCMCglmm simulation)


## Run simulation loop
loop_start_time <- Sys.time()
fprint("Simulating phenotypes, start time: {round(loop_start_time)}")

for (i in 1:length(va_values)) {
    va <- va_values[i]
    h2 <- va/5
    phens <- matrix(0, nrow = nrow(fam), ncol = N)
    
    for (n in 1:N) {
        y <- do_hibayes_sim(h2, k)
        phens[,n] <- y
        
        if (n %% 1 == 0) {
            completed <- N * (i - 1) + n
            max       <- N * length(va_values)
            fprint("Completed simulation {n} for va value {i} (h^2 = {h2}). {completed} out of {max} total phenotypes simulated ({round(completed/max * 100, 2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")
        }
    }
    
    colnames(phens) <- paste0(as.character(h2), "_", as.character(1:N))
    all_phens <- cbind(all_phens, phens)
}
##-------------------------------------------------------------------------------------------




## Save simulated phenotypes
write.csv(all_phens, "method_comparison/simulated_phens.csv", row.names = F)
write.table(all_phens, "method_comparison/simulated_phens", row.names = F, col.names = F)
print_end_message(start_time)
