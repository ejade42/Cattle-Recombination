# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
source("scripts/functions.r")



## Set up background variables
methods <- c("BayesA", "BayesB", "BayesBpi", "BayesC", "BayesCpi", "BayesL", "BayesR", "BayesRR", "BSLMM", "gcta", "gctb", "ldak")
directory = "method_comparison/"

va_values <- c(seq(4, 2, -0.25), seq(1.9, 0.1, -0.1), 0.05, 0.005)
N <- 10
va_vec <- rep(va_values, each = 10)
n_vec  <- rep(1:N, times = 30)
true_h2_vec <- round(va_vec / 5, 3)
method_comparison.df <- NULL



## Read data and combine it all into one dataframe
for (method in methods) {
    filename <- fstring("iterations_{method}.csv")
    data   <- read.csv(paste0(directory, filename))
    h2_est <- data[,"est_h2"]
    h2_err <- data[,"h2_error"]
    
    na_vec <- rep(NA, 300 - nrow(data))
    h2_est <- c(h2_est, na_vec)
    h2_err <- c(h2_err, na_vec)
        
    template.df <- data.frame(va = va_vec, n = n_vec, true_h2 = true_h2_vec)
    template.df[,"h2_est"] <- h2_est
    template.df[,"h2_err"] <- h2_err
    template.df[,"method"] <- method
    
    method_comparison.df <- rbind(method_comparison.df, template.df)
}

write.csv(method_comparison.df, "files/heritability_methods_comparison.csv", row.names = F)



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

