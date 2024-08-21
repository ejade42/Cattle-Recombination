# Runs in main directory

# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
library(lme4)       ## v1.1-34 
library(car)        ## v3.1-2 
library(effects)    ## v4.2-2 
library(MuMIn)      ## v1.47.5 
source("scripts/functions.r")

print_blank()
fprint("Start time: {round(start_time)}")


# Load in data
fprint("Loading recombination data ({round(Sys.time())})")
recombinations.df <- read.csv("files/recombinations_filtered.csv")


# Fit big model
model_fit_start_time <- Sys.time()
fprint("Fitting model ({round(model_fit_start_time)})")
mixed_model_full <- lmer(recombination_rate ~ 1 + (sire_age_at_collection + sire_jersey_propn + kept_markers_adjusted + as.Date(sire_birth_date) + sire_breeding_window)^2 + I(sire_age_at_collection^2) + (1 | parent_id) + (1 | sire_region) + (1 | sire_birth_year), data = recombinations.df)
fprint("Model fitted in {get_time_diff_min(model_fit_start_time, Sys.time())} min")


# Print model results
print_blank()
fprint("Model results:")
summary(mixed_model_full)
Anova(mixed_model_full)


# Dredge models
dredge_start_time <- Sys.time()
print_blank()
fprint("Dredging models ({round(dredge_start_time)})")
suppressWarnings({
options(na.action = "na.fail")
all.fits <- dredge(mixed_model_full, trace = 2)
})
print_blank()
fprint("Dredge completed in {get_time_diff_min(dredge_start_time, Sys.time())} min")


# Write output to csv
print_blank()
fprint("Writing output ({round(Sys.time())})")
write.csv(all.fits, "files/dredged_global_models.csv")
print_end_message(start_time)