# Runs in main directory
# Fits per-sire models estimating key traits, notably his average recombination rate and his rate of change with age

## load packages
start_time <- Sys.time()
library(tidyverse)
source("scripts/functions.r")
print_blank()
fprint("Start time: {round(start_time)}")


## read in data
recombinations.df <- read.csv("files/recombinations_filtered.csv")
sires <- unique(recombinations.df$parent_id)


## create blank variables
mean                 <- NULL
mean_error           <- NULL
mean_sigma           <- NULL
lm_intercept         <- NULL
lm_intercept_error   <- NULL
lm_slope             <- NULL
lm_slope_error       <- NULL
lm_sigma             <- NULL
sire_jersey_propn    <- NULL
sire_breed           <- NULL
halfsibs_in_data     <- NULL
halfsibs_total       <- NULL
sire_breeding_window <- NULL
rounded_ages         <- NULL
children_per_age     <- NULL
sire_birth_date      <- NULL
sire_birth_year      <- NULL
sire_region          <- NULL
sire_avg_age_at_col  <- NULL


## fit null and linear models for each sire, save key info
for (i in 1:length(sires)) {
    his_data <- recombinations.df[recombinations.df$parent_id == sires[i],]
    
    null_model <- lm(recombination_rate ~ 1, data = his_data)
    age_model  <- lm(recombination_rate ~ sire_age_at_collection, data = his_data)
    null_model_coef <- summary(null_model)$coefficients
    age_model_coef  <- summary(age_model)$coefficient
    if (nrow(age_model_coef) == 1) { 
        age_model_coef <- rbind(age_model_coef, sire_age_at_collection = rep(NA, 4))
    }
    
    mean                 <- c(mean,               null_model_coef[1,1])
    mean_error           <- c(mean_error,         null_model_coef[1,2])
    mean_sigma           <- c(mean_sigma,         sigma(null_model))
    lm_intercept         <- c(lm_intercept,       age_model_coef[1,1])
    lm_intercept_error   <- c(lm_intercept_error, age_model_coef[1,2])
    lm_slope             <- c(lm_slope,           age_model_coef[2,1])
    lm_slope_error       <- c(lm_slope_error,     age_model_coef[2,2])
    lm_sigma             <- c(lm_sigma,           sigma(age_model))
    halfsibs_in_data     <- c(halfsibs_in_data,     nrow(his_data))
    halfsibs_total       <- c(halfsibs_total,       his_data[1, "num_halfsibs"])
    sire_jersey_propn    <- c(sire_jersey_propn,    his_data[1, "sire_jersey_propn"])
    sire_breed           <- c(sire_breed,           his_data[1, "sire_breed"])
    sire_breeding_window <- c(sire_breeding_window, his_data[1, "sire_breeding_window"])
    rounded_ages         <- c(rounded_ages,         his_data[1, "sire_rounded_ages"])
    children_per_age     <- c(children_per_age,     his_data[1, "sire_children_per_age"])
    sire_birth_date      <- c(sire_birth_date,      his_data[1, "sire_birth_date"])
    sire_birth_year      <- c(sire_birth_year,      his_data[1, "sire_birth_year"])
    sire_region          <- c(sire_region,          his_data[1, "sire_region"])
    sire_avg_age_at_col  <- c(sire_avg_age_at_col,  mean(his_data$sire_age_at_collection))
    
    if (i %% 250 == 0) {fprint("Completed up to sire {i} of {length(sires)} ({round(i / length(sires) * 100, 2)}% complete). Total time elapsed: {get_time_diff_min(start_time, Sys.time())} min")}
}
family <- rep(1, length(sires))


## combine into dataframe, write output to csv
sire_model_data <- data.frame(family, sires, mean, mean_error, mean_sigma, lm_intercept, lm_intercept_error, lm_slope, lm_slope_error, lm_sigma, sire_jersey_propn, sire_breed, sire_birth_date, sire_birth_year, sire_region, halfsibs_in_data, halfsibs_total, sire_breeding_window, sire_avg_age_at_col, rounded_ages, children_per_age)
write.csv(sire_model_data, "files/sire_models.csv", row.names = F)
write.table(sire_model_data, "files/sire_models", row.names = F, col.names = F)

print_end_message(start_time)