# runs in qc directory
# takes .fam file, and sets sex of each individual to 2 (female) if that individual is present in the "dam" column and 1 (male) otherwise.

filename <- commandArgs(trailingOnly = TRUE)
if (length(filename) == 0) {stop("Please provide a filename when calling the script")}


start_time <- Sys.time()
source("../scripts/functions.r")   ## load common functions
fprint("Start time: {round(start_time)}")
fprint("Updating sex IDs for: {filename}")


original_fam <- read.table(filename)
colnames(original_fam) <- c("family_id", "individual_id", "sire", "dam", "sex", "phenotype")

original_fam[,"sex"] <- rep(1, nrow(original_fam))

for (i in 1:nrow(original_fam)) {
    position <- which(original_fam[i,"dam"] == original_fam[,"individual_id"])
    if (length(position) != 0) {
        original_fam[position[[1]],"sex"] <- 2
    }
    if (i %% 10000 == 0) {
        fprint("Sex ID updated for individual {i} of {nrow(original_fam)} ({round(i/nrow(original_fam) * 100,2)}% complete). Time elapsed: {get_time_diff_min(start_time, Sys.time())} min")
    }
}



write.table(original_fam, filename, row.names = F, col.names = F)
print_end_message(start_time)