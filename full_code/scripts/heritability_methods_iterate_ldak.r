## Runs in method_comparison

# function takes output of readLines(), and extracts the h2 and sd value from the .reml file
extract_values <- function(results) {
    line <- strsplit(results[20], " ", fixed = T)
    h2 <- as.numeric(line[[1]][2])
    sd <- as.numeric(line[[1]][3])

    c(h2, sd)
}


## READ DATA AND ARGUMENT
iterations_ldak <- read.csv("iterations_ldak.csv")
results_file    <- readLines("sires_results.reml")
results         <- extract_values(results_file)
i <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(i) == 0) {stop("Please provide an iteration number when calling the script")}


## SIMULATION PARAMETERS (MUST MATCH power_analysis_simulations.r)
va_values <- c(seq(4, 2, -0.25), seq(1.9, 0.1, -0.1), 0.05, 0.005)
N <- 10


## EXTRACT VALUES TO STORE
va_vec <- rep(va_values, each = 10)
n_vec  <- rep(1:N, times = 30)

va <- va_vec[i]
n  <- n_vec[i]
true_h2  <- round(va / 5, 3)
est_h2   <- results[1]
sd       <- results[2]
h2_error <- est_h2 - true_h2



## SAVE DATA TO CSV
row <- c(va, n, true_h2, est_h2, sd, h2_error)
iterations_ldak[nrow(iterations_ldak)+1,] <- row
write.csv(iterations_ldak, "iterations_ldak.csv", row.names = F)
