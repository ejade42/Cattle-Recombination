## Runs in method_comparison

# function takes output of readLines(), and extracts the h2, se, and p value from the .hsq file
extract_values <- function(results) {
    line <- strsplit(results[5], "\t", fixed = T)
    h2 <- as.numeric(line[[1]][2])
    se <- as.numeric(line[[1]][3])

    line <- strsplit(results[10], "\t", fixed = T)
    p <- as.numeric(line[[1]][2])
    c(h2, se, p)
}


## READ DATA AND ARGUMENT
gcta_iterations <- read.csv("iterations_gcta.csv")
results_file    <- readLines("sires_results.hsq")
results <- extract_values(results_file)
i <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(i) == 0) {stop("Please provide an iteration number when calling the script")}


## SIMULATION PARAMETERS (MUST MATCH heritability_methods_simulation.r)
va_values <- c(seq(4, 2, -0.25), seq(1.9, 0.1, -0.1), 0.05, 0.005)
N <- 10


## EXTRACT VALUES TO STORE
va_vec <- rep(va_values, each = 10)
n_vec  <- rep(1:N, times = 30)

va <- va_vec[i]
n  <- n_vec[i]
true_h2  <- round(va / 5, 3)
est_h2   <- results[1]
p        <- results[3]
h2_error <- est_h2 - true_h2


## SAVE DATA TO CSV
row <- c(va, n, true_h2, est_h2, p, h2_error)
gcta_iterations[nrow(gcta_iterations)+1,] <- row
write.csv(gcta_iterations, "iterations_gcta.csv", row.names = F)
