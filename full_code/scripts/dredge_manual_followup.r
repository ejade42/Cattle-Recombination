## Runs in main directory
## The MuMIn dredge has picked the best combination of fixed effects (age & breed)
## Now we need to pick which random effects to include
## We will do this by fitting the eight possible models, and picking the best via AICc


# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
library(lme4)       ## v1.1-34 
library(MuMIn)      ## v1.47.5 
source("scripts/functions.r")
print_blank()
fprint("Start time: {round(start_time)}")


# Load in data
fprint("Loading recombination data ({round(Sys.time())})")
recombinations.df <- read.csv("files/recombinations_filtered.csv")


# Fit models
fprint("Fitting models ({round(Sys.time())})")
three_rand <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_region) + (1 | sire_birth_year), data = recombinations.df)
two_rand_1 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_region), data = recombinations.df)
two_rand_2 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = recombinations.df)
two_rand_3 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | sire_region) + (1 | sire_birth_year), data = recombinations.df)
one_rand_1 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | parent_id), data = recombinations.df)
one_rand_2 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | sire_region), data = recombinations.df)
one_rand_3 <- lmer(recombination_rate ~ sire_age_at_collection + sire_jersey_propn + (1 | sire_birth_year), data = recombinations.df)
no_rand    <- lm(recombination_rate ~ sire_age_at_collection + sire_jersey_propn, data = recombinations.df)


# Compare models
fprint("Comparing models ({round(Sys.time())})")
comparison <- AICc(three_rand, two_rand_1, two_rand_2, two_rand_3, one_rand_1, one_rand_2, one_rand_3, no_rand)
comparison
comparison[order(comparison$AICc),]

write.csv(comparison, "files/dredge_random_comparison.csv", row.names = F)