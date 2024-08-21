# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
library(lme4)       ## v1.1-34 
library(car)        ## v3.1-2 
library(effects)    ## v4.2-2 
library(gglm)       ## v1.0.2
source("scripts/functions.r")



## LOAD DATA AND PRINT SUMMARY STATISTICS
recombinations.df <- read.csv("input/recombinations_filtered.csv")
year_scale <- 365
total_morgans_genome <- max(recombinations.df$kept_distance_M)
recombinations.df$recombinations <- recombinations.df$recombination_rate * total_morgans_genome

fprint("There are {nrow(recombinations.df)} children remaining")
fprint("There are {length(unique(recombinations.df$parent_id))} sires remaining")
fprint("The average recombination rate is {mean(recombinations.df$recombination_rate)}")
fprint("The total genetic distance is {round(total_morgans_genome, 3)} Morgans")



## PLOT GRAPH (FIG. 2)
graph_model <- lmer(recombinations ~ 1 + sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = recombinations.df)
effects_recombination <- effect(term = "sire_age_at_collection", mod = graph_model, xlevels=5)
x_recom <- as.data.frame(effects_recombination)
line_colour = "#FFB500"
point_colour = "#204085"

recombination_plot <- ggplot() +
    geom_line(data = x_recom[1:4,], aes(x=sire_age_at_collection/year_scale, y=fit), colour = point_colour, linewidth = 2.5, alpha = 0.5) +
    geom_point(data = recombinations.df, aes(sire_age_at_collection/year_scale, recombinations), colour = point_colour, alpha = 0.1, size = 0.75) +
    geom_line(data = x_recom, aes(x=sire_age_at_collection/year_scale, y=fit), colour = line_colour, linewidth = 1.25) +
    xlab("Sire age at collection (years)") + ylab("Expected autosomal recombinations across genome\n(inferred from recombination rate)") +
    scale_x_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,11)) +
    scale_y_continuous(limits = c(0,50)) + theme_bw() + theme(panel.grid = element_line(size=0.15))
recombination_plot
ggsave("output/global_age_effect.png", dpi = 600, width = 7, height = 5)



## CREATE DATAFRAME FOR TABLE OUTPUT
output_table <- data.frame(breed = rep(c("All", "Hol", "Jer"), each = 3), n = rep(0, times = 9), variable = rep(c("Intercept", "Age", "Breed"), times = 3), estimate = rep(0, times = 9), lower = rep(0, times = 9), upper = rep(0, times = 9), p = rep(c(NA, 0, 0), times = 3))

## PRINT MODEL OUTPUT, ALL CATTLE
print_blank(); print_blank(); print_blank(); print_blank()
fprint("FINAL MODEL OUTPUT: ALL CATTLE")
fprint("There are {nrow(recombinations.df)} children remaining")
fprint("There are {length(unique(recombinations.df$parent_id))} sires remaining")
final_model <- lmer(recombination_rate ~ 1 + sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = recombinations.df)
(coef  <- summary(final_model)$coef)
(anova <- Anova(final_model))

confint <- confint(final_model)
print_blank(); fprint("Confidence intervals: normal (recombination rate per day):")
confint[4:6,]
print_blank(); fprint("Age confidence interval: scaled by year (recombination rate per year):")
confint[5,] * 365
print_blank(); fprint("Age confidence interval: scaled by year and morgans (recombinations per year):")
confint[5,] * 365 * total_morgans_genome

output_table[1:3, "n"]                <- nrow(recombinations.df)
output_table[1:3, "estimate"]         <- coef[1:3, 1]
output_table[1:3, c("lower","upper")] <- confint[4:6, 1:2]
output_table[2:3, "p"]                <- anova[1:2, 3]


## PRINT MODEL OUTPUT, HOLSTEIN
print_blank(); print_blank(); print_blank(); print_blank()
fprint("FINAL MODEL OUTPUT: HOLSTEIN ONLY")
hol <- filter(recombinations.df, sire_jersey_propn <= 1/16)
fprint("There are {nrow(hol)} children remaining")
fprint("There are {length(unique(hol$parent_id))} sires remaining")
final_model <- lmer(recombination_rate ~ 1 + sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = hol)
(coef  <- summary(final_model)$coef)
(anova <- Anova(final_model))

confint <- confint(final_model)
print_blank(); fprint("Confidence intervals: normal (recombination rate per day):")
confint[4:6,]
print_blank(); fprint("Age confidence interval: scaled by year (recombination rate per year):")
confint[5,] * 365
print_blank(); fprint("Age confidence interval: scaled by year and morgans (recombinations per year):")
confint[5,] * 365 * total_morgans_genome

output_table[4:6, "n"]                <- nrow(hol)
output_table[4:6, "estimate"]         <- coef[1:3, 1]
output_table[4:6, c("lower","upper")] <- confint[4:6, 1:2]
output_table[5:6, "p"]                <- anova[1:2, 3]


## PRINT MODEL OUTPUT, JERSEY
print_blank(); print_blank(); print_blank(); print_blank()
fprint("FINAL MODEL OUTPUT: JERSEY ONLY")
jer <- filter(recombinations.df, sire_jersey_propn >= 15/16)
fprint("There are {nrow(jer)} children remaining")
fprint("There are {length(unique(jer$parent_id))} sires remaining")
final_model <- lmer(recombination_rate ~ 1 + sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = jer)
(coef  <- summary(final_model)$coef)
(anova <- Anova(final_model))

confint <- confint(final_model)
print_blank(); fprint("Confidence intervals: normal (recombination rate per day):")
confint[4:6,]
print_blank(); fprint("Age confidence interval: scaled by year (recombination rate per year):")
confint[5,] * 365
print_blank(); fprint("Age confidence interval: scaled by year and morgans (recombinations per year):")
confint[5,] * 365 * total_morgans_genome

output_table[7:9, "n"]                <- nrow(jer)
output_table[7:9, "estimate"]         <- coef[1:3, 1]
output_table[7:9, c("lower","upper")] <- confint[4:6, 1:2]
output_table[8:9, "p"]                <- anova[1:2, 3]


## Save model information to csv
output_table[, 4:7] <- signif(output_table[, 4:7], 4)
write.csv(output_table, "output/global_age_effect.csv", row.names = F)


## Save gglm diagnostic plot of the all-cattle model
print_blank(); print_blank(); print_blank(); print_blank()
fprint("Saving all-cattle model diagnostic plot via gglm...")
final_model <- lmer(recombination_rate ~ 1 + sire_age_at_collection + sire_jersey_propn + (1 | parent_id) + (1 | sire_birth_year), data = recombinations.df)
gglm::gglm(final_model, theme = theme_bw())
ggsave("output/global_age_effect_gglm.png", dpi = 600, width = 8, height = 8)    ## save diagnostic plot, Fig. S1


# End message
print_blank()
print_end_message(start_time)