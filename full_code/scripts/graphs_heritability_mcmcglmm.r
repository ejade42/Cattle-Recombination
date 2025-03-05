## Runs in main cows directory


path_to_files <- "heritability/mcmcglmm/"
all_filenames <- list.files(path = path_to_files)


construct_filename <- function(response, breed, filter) {
    return(paste0("mcmcglmm_", response, "_", breed, "_", sprintf('%03d', filter), ".csv"))
}

measures_per_response <- c("_h2", "_sd", "_lower", "_upper")
columns <- c("sire_breeding_window", "breed", "n_sires", "n_offspring_for_repeated", paste0("repeated", measures_per_response), paste0("mean", measures_per_response), paste0("slope", measures_per_response))
heritability_analysis_mcmcglmm <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(heritability_analysis_mcmcglmm) <- columns


for (breed in c("all", "holstein", "jersey", "subsample")) {
    for (filter in 0:500) {
        i <- nrow(heritability_analysis_mcmcglmm) + 1
        heritability_analysis_mcmcglmm[i, "sire_breeding_window"] <- filter
        heritability_analysis_mcmcglmm[i, "breed"]                <- breed
        
        for (response in c("repeated", "mean", "slope")) {
            filename <- construct_filename(response, breed, filter)
            if (file.exists(paste0(path_to_files, filename))) {
                new_data <- read.csv(paste0(path_to_files, filename))
                
                if (!is.na(new_data$num_sires)) {
                    heritability_analysis_mcmcglmm[i, "n_sires"] <- new_data$num_sires
                }
                if (response == "repeated") {
                    heritability_analysis_mcmcglmm[i, "n_offspring_for_repeated"] <- new_data$num_offspring
                }
                heritability_analysis_mcmcglmm[i, paste0(response, "_h2")]      <- new_data$mean       
                heritability_analysis_mcmcglmm[i, paste0(response, "_sd")]      <- new_data$sd
                heritability_analysis_mcmcglmm[i, paste0(response, "_lower")]   <- new_data$quant_0.025
                heritability_analysis_mcmcglmm[i, paste0(response, "_upper")]   <- new_data$quant_0.975
            }
        
        }
    }     
}


write.csv(heritability_analysis_mcmcglmm, "files/heritability_analysis_mcmcglmm.csv", row.names = F)



## Make heritability plot
library(tidyverse)  ## v2.0.0
heritability <- read.csv("files/heritability_analysis_mcmcglmm.csv")
heritability[heritability$breed == "all", "graph_breed"]       <- "(a) All"
heritability[heritability$breed == "subsample", "graph_breed"] <- "(b) All, randomly subsampled"
heritability[heritability$breed == "holstein", "graph_breed"]  <- "(c) Holstein"
heritability[heritability$breed == "jersey", "graph_breed"]    <- "(d) Jersey"


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
