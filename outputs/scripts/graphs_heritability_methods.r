# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
source("scripts/functions.r")


## Read in data
method_comparison.df <- read.csv("input/heritability_methods_comparison.csv")


## Plot the three panels (different visualisations of the comparison) - Supplementary Figure S2a,b,c
ggplot(method_comparison.df, aes(x = method, y = h2_err)) +
    geom_jitter(aes(colour = method), alpha = 0.1, size = 1.25) +
    geom_boxplot(aes(colour = method), outlier.shape = NA) +
    geom_hline(aes(yintercept=0), col = "black", alpha = 0.9, linewidth = 0.25) +
    xlab("Heritability estimation method") + ylab("Error (estimated h^2 - true h^2)") +
    theme_bw()
ggsave("output/heritability_methods_1.png", dpi = 600, width = 9, height = 7)


ggplot(method_comparison.df, aes(x = method, y = h2_err)) +
    geom_boxplot(aes(colour = method), outlier.shape = NA) +
    xlab("Heritability estimation method") + ylab("Error (estimated h^2 - true h^2)") +
    geom_hline(aes(yintercept=0), col = "black", alpha = 0.9, linewidth = 0.25) +
    theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_wrap("true_h2")
ggsave("output/heritability_methods_2.png", dpi = 600, width = 12, height = 10)


ggplot(method_comparison.df, aes(x = true_h2, y = h2_err)) +
    geom_point(aes(colour = method), alpha = 0.75, size = 0.7, show.legend = F) +
    geom_smooth(col = "darkblue", se = F) +
    geom_hline(aes(yintercept=0), col = "black", alpha = 0.9, linewidth = 0.25) +
    facet_wrap("method", nrow = 4, ncol = 3) +
    xlab("True h^2 (simulation parameter)") + ylab("Error (estimated h^2 - true h^2)") +
    theme_bw()
ggsave("output/heritability_methods_3.png", dpi = 600, width = 8, height = 10)

