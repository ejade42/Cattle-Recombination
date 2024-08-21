# runs in main project directory, reads and summarises the LINKPHASE outputs

## SETUP
##----------------------------------------------------------------------------------------------------------------------------------
start_time <- Sys.time()
source("scripts/functions.r")    ## load common functions
print_blank()
fprint("Start time: {round(start_time)}")

# load settings
settings <- readLines("0)settings.txt")

# read in data
recombinations.df <- read.table(fstring("LINKPHASE_outputs/chr_1_nrec_hmm.txt"))[,1:6]
colnames(recombinations.df) <- c("individual_id", "parent_id", "parent_sex", "num_halfsibs", "parent_phasing", "mate_genotype")
num_chromosomes <- 29
##----------------------------------------------------------------------------------------------------------------------------------




## READ IN RECOMBINATIONS AND INFORMATIVE MARKERS
##----------------------------------------------------------------------------------------------------------------------------------
# Extract recombinations
for (i in 1:num_chromosomes) {
    filename <- fstring("chr_{i}_nrec_hmm.txt")
    path     <- fstring("LINKPHASE_outputs/{filename}")
    col_name <- fstring("recombinations_chr{i}")
    recombinations <- read.table(path)[,10]
    recombinations.df[col_name] <- recombinations
}
print("read in recombinations", quote=F)

# Extract informative markers
for (i in 1:num_chromosomes) {
    filename <- fstring("chr_{i}_nrec_hmm.txt")
    path     <- fstring("LINKPHASE_outputs/{filename}")
    col_name <- fstring("informative_markers_chr{i}")
    recombinations <- read.table(path)[,9]
    recombinations.df[col_name] <- recombinations
}
print("read in informative markers", quote=F)

col_first_recombinations_chri <- which(colnames(recombinations.df) == "recombinations_chr1")
col_final_recombinations_chri <- which(colnames(recombinations.df) == fstring("recombinations_chr{num_chromosomes}"))
recombinations.df["total_recombinations"] <- rowSums(recombinations.df[,col_first_recombinations_chri:col_final_recombinations_chri])
col_first_informative_markers_chri <- which(colnames(recombinations.df) == "informative_markers_chr1")
col_final_informative_markers_chri <- which(colnames(recombinations.df) == fstring("informative_markers_chr{num_chromosomes}"))
recombinations.df["total_informative_markers"] <- rowSums(recombinations.df[,col_first_informative_markers_chri:col_final_informative_markers_chri])
print("merged", quote=F)
##----------------------------------------------------------------------------------------------------------------------------------




## CALCULATE MIN CROSSOVER DISTANCES
##----------------------------------------------------------------------------------------------------------------------------------
iterations       <- num_chromosomes
threshold        <- as.numeric(settings[4])
distance_measure <- as.character(settings[2])
fprint('OPTION: calculating distances using "{distance_measure}", with double-crossover threshold {threshold} cM.')
loop_start_time <- Sys.time()
print_blank()
fprint("Processing distances between crossovers, start time: {round(loop_start_time)}")

recombinations.df["min_crossover_distance"]  <- rep(Inf, nrow(recombinations.df))
recombinations.df["min_distance_chromosome"] <- rep(NA, nrow(recombinations.df))
recombinations.df["total_below_threshold"]   <- rep(0, nrow(recombinations.df))

for (i in 1:iterations) {
    col_name <- fstring("min_chr{i}")
    recombinations.df[col_name] <- rep(Inf, nrow(recombinations.df))
}
total_iterations <- iterations * nrow(recombinations.df)
for (i in 1:iterations) {
    filename <- fstring("chr_{i}_recombination_distances.csv")
    path     <- fstring("LINKPHASE_outputs/{filename}")
    distances.df <- read.csv(path)
    col_name     <- fstring("min_chr{i}")
    col_name_2   <- fstring("num_below_threshold_chr{i}")
    recombinations.df[col_name_2] <- rep(0, nrow(recombinations.df))
    for (j in 1:nrow(recombinations.df)) {
        individual_id <- recombinations.df[j,"individual_id"]
        parent_id <- recombinations.df[j,"parent_id"]
        min_distance_stored <- recombinations.df[j,"min_crossover_distance"]
        min_distance_chromosome <- recombinations.df[j,"min_distance_chromosome"]
        min_per_chr <- Inf
        if (individual_id %in% distances.df$individual_id & parent_id %in% distances.df$parent_id) {
            distances_filtered.df <- distances.df[distances.df$individual_id == individual_id & distances.df$parent_id == parent_id,]
            min_per_chr  <- min(distances_filtered.df[,distance_measure], na.rm = T)
            if (min_per_chr < min_distance_stored) {
                recombinations.df[j,"min_crossover_distance"]  <- min_per_chr
                recombinations.df[j,"min_distance_chromosome"] <- i
            }
            recombinations.df[j,col_name]   <- min_per_chr
            recombinations.df[j,col_name_2] <- sum(distances_filtered.df[,distance_measure] < threshold, na.rm = T)
        }
        if (j %% 20000 == 0) {
            this_iteration   <- (i-1)*nrow(recombinations.df) + j
            fprint("chr {i} parent-child pair {j} (iteration {this_iteration} of {total_iterations} = {round(this_iteration/total_iterations * 100,2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")
            fprint("min_per_chr: {min_per_chr}, global_min: {recombinations.df[j,'min_crossover_distance']} chr {recombinations.df[j,'min_distance_chromosome']}")
        }
    }
}
print_blank()

# Sum to find total number of crossovers closer together than the specified threshold
col_first_belowthreshold_chri <- which(colnames(recombinations.df) == "num_below_threshold_chr1")
col_final_belowthreshold_chri <- which(colnames(recombinations.df) == fstring("num_below_threshold_chr{iterations}"))
recombinations.df["total_below_threshold"] <- rowSums(recombinations.df[,col_first_belowthreshold_chri:col_final_belowthreshold_chri])
##----------------------------------------------------------------------------------------------------------------------------------




## INTEGRATE SIRE DATE AND BREED INFORMATION
##----------------------------------------------------------------------------------------------------------------------------------
# Load and subset the offspring batch collection data
print("reading in birth date and batch collection date data", quote=F)
offspring.df <- read.csv(fstring("input/matched_mating_data_final_117621_anmls.csv"))
offspring.df <- offspring.df[,c("progeny_lic_animal_key", "batch_collection_date", "mating.date")]
colnames(offspring.df) <- c("long_id", "batch_collection_date", "mating_date")

# Load and subset the sire birth date data
sires.df <- read.csv(fstring("input/sire_(3151_animals_of_interest)_locations.csv"))
sires.df <- sires.df[,c("lic_animal_key", "birth_date", "statistical_region_descr")]
colnames(sires.df) <- c("long_id", "birth_date", "region")

# Load and subset the sire breed data
breeds.df <- read.csv(fstring("input/bull_breed_mix.csv"))
breeds.df <- breeds.df[,c("lic_animal_key", "prop_Jersey")]
colnames(breeds.df) <- c("long_id", "propn_jersey")

# Load and subset the big and small ID numbers we used for PLINK
plink_ids.df <- read.table(fstring("files/newID4plink_01"))
plink_ids.df <- plink_ids.df[,c(4,1)]
colnames(plink_ids.df) <- c( "individual_id", "long_id")

# Merge PLINK IDs with birth/batch collection dates and breed data
print("merging data", quote=F)
offspring_merged <- merge(plink_ids.df, offspring.df)[,c("individual_id", "batch_collection_date", "mating_date")]
sires_merged     <- merge(plink_ids.df, sires.df)[,c("individual_id", "birth_date", "region")]
breeds_merged    <- merge(plink_ids.df, breeds.df)[,c("individual_id", "propn_jersey")]
colnames(sires_merged)  <- c("parent_id", "sire_birth_date", "region")
colnames(breeds_merged) <- c("parent_id", "propn_jersey")

# Integrate dates and breeds into recombinations dataset
print_blank()
loop_start_time <- Sys.time()
fprint("Processing date and breeding information, start time: {round(loop_start_time)}")
for (i in 1:nrow(recombinations.df)) {
    parent_id <- recombinations.df[i,"parent_id"]
    if (parent_id %in% sires_merged$parent_id) {
        sire_birth_date <- sires_merged[sires_merged$parent_id == parent_id,2][1]
        sire_region     <- sires_merged[sires_merged$parent_id == parent_id,3][1]
        recombinations.df[i,"sire_birth_date"] <- dob <- as.Date(sire_birth_date, format = "%d%b%Y")
        recombinations.df[i,"sire_region"]     <- sire_region
        recombinations.df[i,"sire_birth_year"] <- format(dob, "%Y")
    }
    individual_id <- recombinations.df[i,"individual_id"]
    if (individual_id %in% offspring_merged$individual_id) {
        batch_collection_date <- offspring_merged[offspring_merged$individual_id == individual_id,2][1]
        mating_date           <- offspring_merged[offspring_merged$individual_id == individual_id,3][1]
        recombinations.df[i,"batch_collection_date"] <- as.Date(batch_collection_date, format = "%Y-%m-%d")
        recombinations.df[i,"mating_date"]     <- as.Date(mating_date, format = "%d/%m/%Y")
        recombinations.df[i,"date_difference"] <- as.Date(mating_date, format = "%d/%m/%Y") - as.Date(batch_collection_date, format = "%Y-%m-%d")
    }
    if (parent_id %in% breeds_merged$parent_id) {
        sire_jersey_proportion <- breeds_merged[breeds_merged$parent_id == parent_id,2][1]
        recombinations.df[i,"sire_jersey_propn"] <- sire_jersey_proportion
    }
    if (i %% 5000 == 0) {
        fprint("dates and breed obtained for parent-child pair {i} of {nrow(recombinations.df)} ({round(i/nrow(recombinations.df) * 100,2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")
    }
}
print_blank()

# Make new variables based on date information
recombinations.df["sire_age_at_collection"] <- as.Date(recombinations.df$batch_collection_date) - as.Date(recombinations.df$sire_birth_date)
recombinations.df["sire_age_scaled"] <- recombinations.df$sire_age_at_collection / 365


# Make new factor variable for sire breed (Holstein/Mix/Jersey)
if (length(settings) >= 26) {breed_threshold <- as.numeric(settings[26])} else {breed_threshold <- 0.25}

recombinations.df$sire_breed <- ifelse(recombinations.df$sire_jersey_propn >= 1-breed_threshold,
                                       "Jersey",
                                       ifelse(recombinations.df$sire_jersey_propn <= breed_threshold, "Holstein", "Mix"))
recombinations.df$sire_breed <- factor(recombinations.df$sire_breed, levels = c("Holstein", "Mix", "Jersey"))


# Calculate breeding window pre-filtering
for (parent in unique(recombinations.df$parent_id)) {
    their_offspring <- recombinations.df[recombinations.df$parent_id == parent,]
    their_ages <- unique(their_offspring$sire_age_at_collection)
    breeding_window <- max(their_ages) - min(their_ages)
    recombinations.df[recombinations.df$parent_id == parent,"breeding_window_pre_filter"] <- breeding_window
}
##----------------------------------------------------------------------------------------------------------------------------------




fprint("All recombination information extracted.")
fprint("{nrow(recombinations.df)} parent-child pairs and {length(unique(recombinations.df$parent_id))} parents present before filtering.")
write.csv(recombinations.df, "files/recombinations.csv", row.names = F)
print_end_message(start_time, F)