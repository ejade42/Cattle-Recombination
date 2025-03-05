# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
source("scripts/functions.r")



## Combine individual hibayes heritability results into one dataframe
hibayes_all <- read.csv("heritability/iterations_hibayes_all.csv")
hibayes_sub <- read.csv("heritability/iterations_hibayes_subsample.csv")
hibayes_hol <- read.csv("heritability/iterations_hibayes_holstein.csv")
hibayes_jer <- read.csv("heritability/iterations_hibayes_jersey.csv")

hibayes_all$breed <- "(a) All"
hibayes_sub$breed <- "(b) All, randomly subsampled"
hibayes_hol$breed <- "(c) Holstein"
hibayes_jer$breed <- "(d) Jersey"


hibayes <- rbind(hibayes_all, hibayes_sub, hibayes_hol, hibayes_jer)
write.csv(hibayes, "files/heritability_analysis_combined.csv", row.names = F)



## Make heritability plot
hibayes <- read.csv("files/heritability_analysis_combined.csv")
nmax <- 2000
override.linetype <- c(1, 1, 2)
ggplot(hibayes, aes(x = sire_breeding_window)) +
    scale_colour_manual(values = c("mean" = "#00BBFF", "slope" = "orange", "sample size" = "darkgrey"),
                        breaks = c("mean", "slope", "sample size")) + 
    labs(colour = "", linetype = "") +
    scale_linetype_manual(values = c("sample size" = 2), guide = F) + 
    geom_line(aes(y = n / nmax, linetype = "sample size", col = "sample size")) +
    geom_line(aes(y = mean_h2, col = "mean")) +
    geom_smooth(aes(y = mean_h2), col = "#005599", se = F, span = 0.5) +
    geom_line(aes(y = slope_h2, col = "slope")) + 
    geom_smooth(aes(y = slope_h2), col = "#AF7C00", se = F, span = 0.5) +    # "#9B0000"
    guides(colour = guide_legend(override.aes = list(linetype = override.linetype))) +
    ylab(bquote("Estimated"~h^2)) + xlab("Sire breeding window filter (days)") + 
    scale_y_continuous(sec.axis = sec_axis(~.*nmax, name="Sires remaining in data"), limits = c(0,1)) +
    scale_x_continuous(limits = c(0, 500)) +
    facet_wrap("breed") +
    theme_bw()
ggsave("output/heritability.png", dpi = 600, width = 10, height = 8)