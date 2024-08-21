# runs inside the chr_i directories inside qc

## Set up
start_time <- Sys.time()
source("../../scripts/functions.r")   ## load common functions
fprint("Start time: {round(start_time)}")

# read chromosome number from provided argument, create filename for the corresponding map file
chr_number <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(chr_number) != 1) {stop("Please provide a chromosome number when calling the script")}
fprint("Working on chromosome {chr_number}")
map_filename <- fstring("chr_{chr_number}.markers")

# load recombinations from Linkphase output and rename columns to something human-readable
recombinations.df <- read.table("recombinations_hmm")
colnames(recombinations.df) <- c("individual_id", "parent_id", "marker_1", "marker_2")
map.df <- read.table(map_filename)

# create variables in recombinations.df
recombinations.df$marker_1_position <- map.df[recombinations.df$marker_1, 3]
recombinations.df$marker_2_position <- map.df[recombinations.df$marker_2, 3]
recombinations.df$crossover_interval_centre <- (recombinations.df$marker_1_position + recombinations.df$marker_2_position)/2
recombinations.df$crossover_interval_width <- recombinations.df$marker_2_position - recombinations.df$marker_1_position
recombinations.df$distance_next_centre <- rep(NA, nrow(recombinations.df))
recombinations.df$distance_lowerlimit_nextupperlimit <- rep(NA, nrow(recombinations.df))
recombinations.df$distance_upperlimit_nextlowerlimit <- rep(NA, nrow(recombinations.df))

# calculate crossover distances for each individual (this script runs per-chromosomes), using three different distance measures
rows <- nrow(recombinations.df)
for (i in 1:(rows-1)) {
    if (recombinations.df[i,"individual_id"] == recombinations.df[i+1,"individual_id"]) {
        recombinations.df[i,"distance_next_centre"] <- recombinations.df[i+1,"crossover_interval_centre"] - recombinations.df[i,"crossover_interval_centre"]
        recombinations.df[i,"distance_lowerlimit_nextupperlimit"] <- recombinations.df[i+1,"marker_2_position"] - recombinations.df[i,"marker_1_position"]
        recombinations.df[i,"distance_upperlimit_nextlowerlimit"] <- recombinations.df[i+1,"marker_1_position"] - recombinations.df[i,"marker_2_position"]
    }
    if (i %% 5000 == 0) {
        fprint("calculated crossover distances for crossover row {i} of {rows} ({round(i/rows * 100,2)}% complete)")
    }
}

# write output to chr_i_recombination_distances.csv file
output_filename <- fstring("chr_{chr_number}_recombination_distances.csv")
write.csv(recombinations.df, output_filename, row.names = F)
print_end_message(start_time)
