library(MCMCglmm)    ## v2.36
library(genio)       ## v1.1.2
library(MASS)        ## v7.3.57
library(corpcor)     ## v1.6.10
library(tidyverse)   ## v2.0.0
source("../scripts/functions.r")


## Read arguments
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2) {stop("Please provide a filter and breed setting when calling the script")}
filter <- as.numeric(arguments[1])
breed <- as.character(arguments[2])

fprint("Filter: {filter} days.")
fprint("Breed: {breed} only.")

## Set MCMCglmm parameters
iteration_number <- 10000
thinning_number  <- 1
burnin_number    <- 2000

## Read GRM
grm <- read_grm("grm_all_sires")$kinship
grm_inv <- grm %>% 
    solve() %>%
    make.positive.definite() %>%
    as(., "dgCMatrix")


## Function to extract quantiles, mean, and sd from model
quantile_choices <- c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1)
extract_heritability_samples <- function(model_object) {
    vcv_genetic  <- model_object$VCV[, "parent_id"]
    vcv_residual <- model_object$VCV[, "units"]

    heritability_samples <- vcv_genetic / (vcv_genetic + vcv_residual)
    output <- quantile(heritability_samples, probs = quantile_choices)
    output <- c(mean(heritability_samples), sd(heritability_samples), output, length(heritability_samples))
    
    names(output) <- c("mean", "sd", paste0("quant_", quantile_choices), "sample_count")
    return(output)
}




## REPEATED MEASURES MODEL
## -------------------------------------------------------------------------------------------------------------------------
## Read data and subset to breed
fprint("Processing repeated measures per sire information...")
recombination_data <- read.csv("../files/recombinations_filtered.csv")
recombination_data$parent_id_dup <- recombination_data$parent_id

fprint("{length(unique(recombination_data$parent_id))} sires and {nrow(recombination_data)} offspring present before filter.")
recombination_data <- recombination_data[recombination_data$sire_breeding_window >= filter,]
fprint("{length(unique(recombination_data$parent_id))} sires and {nrow(recombination_data)} offspring present after filter.")


if (tolower(breed) == "all") {
    invisible()
} else if (tolower(breed) == "subsample") {
    num_sires <- length(unique(recombination_data[recombination_data$sire_breeding_window >= 365, "parent_id"]))

    sires_to_keep <- sample(unique(recombination_data$parent_id), num_sires, replace = F)
    recombination_data <- recombination_data[recombination_data$parent_id %in% sires_to_keep,]
    fprint("Sires randomly sampled (without replacement) down to {num_sires} sires, the number left from all sires at filter = 365.")
    if (filter > 365) {
        stop("Subsample only runs up to filter = 365. Ending process.")
    }
} else if (tolower(breed) == "holstein") {
    recombination_data <- recombination_data[recombination_data$sire_jersey_propn <= 1/16,]
} else if (tolower(breed) == "jersey") {
    recombination_data <- recombination_data[recombination_data$sire_jersey_propn >= 15/16,]

} else {stop('Please provide "all", "subsample", "jersey", or "holstein" as the breed')}

num_sires     <- length(unique(recombination_data$parent_id))
num_offspring <- nrow(recombination_data)
fprint('{num_sires} sires  of the breed "{breed}" remain, with {num_offspring} offspring in this data.')



## Fit repeated measures model
repeated_model <- MCMCglmm(fixed = recombination_rate ~ sire_jersey_propn + sire_breeding_window,
                           random = ~ parent_id + parent_id_dup + sire_birth_year + sire_region,
                           data = recombination_data, ginverse = list(parent_id = grm_inv),
                           nitt = iteration_number, thin = thinning_number, burnin = burnin_number,
                           pr = TRUE, pl = TRUE)

repeated_output <- extract_heritability_samples(repeated_model)
information <- c(breed = breed, breeding_window = filter, num_sires = num_sires, num_offspring = num_offspring)


## Write output to file
write.csv(data.frame(t(c(information, repeated_output))), fstring("mcmcglmm/mcmcglmm_repeated_{breed}_{sprintf('%03d', filter)}.csv"), row.names = F)
## -------------------------------------------------------------------------------------------------------------------------





## -------------------------------------------------------------------------------------------------------------------------
## Read data and subset to breed
fprint("Processing single measures per sire information...")
sire_models <- read.csv("../files/sire_models.csv")
fprint("{nrow(sire_models)} sires present before filter.")
sire_models <- sire_models[sire_models$sire_breeding_window >= filter,]
fprint("{nrow(sire_models)} sires remain after filter.")

if (tolower(breed) == "all") {
    invisible()
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

num_sires <- nrow(sire_models)
fprint('{num_sires} sires remain of the breed "{breed}".')


sire_models <- sire_models[,c("sires", "mean", "mean_sigma", "lm_slope", "sire_jersey_propn", "sire_birth_year", "sire_region", "halfsibs_total", "sire_breeding_window")]
colnames(sire_models)[2:4] <- c("mean", "sd", "slope")
sire_models$parent_id <- sire_models$parent_id_dup <- sire_models$sires



## Fit mean model
mean_model <- MCMCglmm(fixed = mean ~ sire_jersey_propn + halfsibs_total + sire_breeding_window,
                       random = ~ parent_id + sire_birth_year + sire_region,
                       data = sire_models, ginverse = list(parent_id = grm_inv),
                       nitt = iteration_number, thin = thinning_number, burnin = burnin_number,
                       pr = TRUE, pl = TRUE)

mean_output <- extract_heritability_samples(mean_model)
information <- c(breed = breed, breeding_window = filter, num_sires = num_sires)

## Write output to file
write.csv(data.frame(t(c(information, mean_output))), fstring("mcmcglmm/mcmcglmm_mean_{breed}_{sprintf('%03d', filter)}.csv"), row.names = F)


## Fit slope model
sire_models <- drop_na(sire_models, slope)
num_sires   <- nrow(sire_models)
fprint('{num_sires} sires remain with non-missing slope information.')

slope_model <- MCMCglmm(fixed = slope ~ sire_jersey_propn + halfsibs_total + sire_breeding_window,
                       random = ~ parent_id + sire_birth_year + sire_region,
                       data = sire_models, ginverse = list(parent_id = grm_inv),
                       nitt = iteration_number, thin = thinning_number, burnin = burnin_number,
                       pr = TRUE, pl = TRUE)

slope_output <- extract_heritability_samples(slope_model)
information <- c(breed = breed, breeding_window = filter, num_sires = num_sires)

## Write output to file
write.csv(data.frame(t(c(information, slope_output))), fstring("mcmcglmm/mcmcglmm_slope_{breed}_{sprintf('%03d', filter)}.csv"), row.names = F)
## -------------------------------------------------------------------------------------------------------------------------
