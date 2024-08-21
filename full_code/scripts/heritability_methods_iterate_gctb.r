## Runs in method_comparison


## READ DATA AND ARGUMENT
gctb_iterations <- read.csv("iterations_gctb.csv")
results         <- read.table("sires_results.parRes")
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
est_h2   <- results[6,1]
h2_error <- est_h2 - true_h2



## SAVE DATA TO CSV
row <- c(va, n, true_h2, est_h2, h2_error)
gctb_iterations[nrow(gctb_iterations)+1,] <- row
write.csv(gctb_iterations, "iterations_gctb.csv", row.names = F)
