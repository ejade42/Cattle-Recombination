## Load in packages/functions
library(tidyverse)  ## v2.0.0


## Read in data
heritability <- read.csv("input/heritability_analysis_mcmcglmm.csv")


## Relabel groups
heritability[heritability$breed == "all", "graph_breed"]       <- "(a) All"
heritability[heritability$breed == "subsample", "graph_breed"] <- "(b) All, randomly subsampled"
heritability[heritability$breed == "holstein", "graph_breed"]  <- "(c) Holstein"
heritability[heritability$breed == "jersey", "graph_breed"]    <- "(d) Jersey"


## Make MCMCglmm heritability plot
nmax <- 2000
override.linetype <- c(1, 1, 1, 2)
ggplot(heritability, aes(x = sire_breeding_window)) +
    scale_colour_manual(values = c("mean" = "#00BBFF", "slope" = "orange", "repeated" = "#00A240", "sample size" = "darkgrey"),
                        breaks = c("mean", "slope", "repeated", "sample size")) + 
    labs(colour = "", linetype = "") +
    scale_linetype_manual(values = c("sample size" = 2), guide = "none") + 
    geom_line(aes(y = n_sires / nmax, linetype = "sample size", col = "sample size")) +
    geom_line(aes(y = mean_h2, col = "mean")) +
    geom_smooth(aes(y = mean_h2), col = "#005599", se = F, span = 0.5) +
    geom_line(aes(y = slope_h2, col = "slope")) + 
    geom_smooth(aes(y = slope_h2), col = "#AF7C00", se = F, span = 0.5) +    # "#9B0000"
    geom_line(aes(y = repeated_h2, col = "repeated")) + 
    geom_smooth(aes(y = repeated_h2), col = "#005420", se = F, span = 0.5) +    # "#9B0000"
    guides(colour = guide_legend(override.aes = list(linetype = override.linetype))) +
    ylab(bquote("Estimated"~h^2)) + xlab("Sire breeding window filter (days)") + 
    scale_y_continuous(sec.axis = sec_axis(~.*nmax, name="Sires remaining in data"), limits = c(0,1)) +
    scale_x_continuous(limits = c(0, 500)) +
    facet_wrap("graph_breed") +
    theme_bw()
ggsave("output/heritability_mcmcglmm.png", dpi = 600, width = 10, height = 8)
