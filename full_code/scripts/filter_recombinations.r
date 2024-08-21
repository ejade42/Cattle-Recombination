# runs in Cows main directory - can run elsewhere, but will read and write recombinations csvs to that location


## SETUP
##----------------------------------------------------------------------------------------------------------------------------------
start_time <- Sys.time()
source("scripts/functions.r")   ## load common functions
print_blank()
fprint("Start time: {round(start_time)}")

# function to print info on how many sires/offspring remain, 
print_info <- function(nrow_previous, sires_previous, chr_previous, step_name, fathers = T) {
    nrow_current <- nrow(recombinations.df)
    difference <- nrow_previous - nrow_current
    message <- "father-child pairs remain."
    if (fathers == F) {message <- "parent-child pairs remain."}
    print(paste(step_name, nrow_current, message, difference, "removed this step."), quote = F)
    
    sires_current <- length(unique(recombinations.df$parent_id))
    difference <- sires_previous - sires_current
    spaces <- paste(rep(" ", nchar(step_name)), collapse = "")
    message <- "sires remain."
    if (fathers == F) {message <- "parents remain."}
    print(paste(spaces, sires_current, message, difference, "removed this step."), quote = F)
    
    chr_current <- sum(recombinations.df$num_kept_chromosomes)
    difference <- chr_previous - chr_current
    spaces <- paste(rep(" ", nchar(step_name)), collapse = "")
    message <- "total chromosomes remain."
    print(paste(spaces, chr_current, message, difference, "removed this step."), quote = F)   
}

# load settings
settings <- readLines("0)settings.txt")

# load recombination data
recombinations.df <- read.csv("files/recombinations.csv")  ## this will be read from current directory, not "Cows"
n_chromosomes <- 29
recombinations.df["num_kept_chromosomes"]     <- rep(n_chromosomes, nrow(recombinations.df))

nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)
print_info(nrow_last_step, sires_last_step, chr_last_step, "BEFORE FILTERING:", F)

##----------------------------------------------------------------------------------------------------------------------------------




## FILTERS FOR FATHERS/PATERNAL EVENTS ONLY, AND CHECK ALL NECESSARY DATE/BREED INFO PRESENT
##----------------------------------------------------------------------------------------------------------------------------------
# remove all maternal events/keep only paternal events
recombinations.df <- recombinations.df[recombinations.df$parent_sex == 1,]
print_info(nrow_last_step, sires_last_step, chr_last_step, "PATERNAL EVENTS ONLY:")
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# filter to only our 3151 sires of interest
sires_of_interest <- read.table("files/sires_of_interest")[,1]
recombinations.df <- recombinations.df[recombinations.df$parent_id %in% sires_of_interest,]
print_info(nrow_last_step, sires_last_step, chr_last_step, "3151 SIRES OF INTEREST ONLY:")
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# Filter to only sires where we have breed information
# This filters to only the sires of interest anyway as we only have breed info for those 3151 sires, just nice to make it more explicit
recombinations.df <- recombinations.df[is.na(recombinations.df$sire_jersey_propn) == F,]
print_info(nrow_last_step, sires_last_step, chr_last_step, "SIRE BREED INFORMATION PRESENT:")
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# Filter to only parent-child pairs with both dates present
recombinations.df <- recombinations.df[is.na(recombinations.df$sire_birth_date) == F & is.na(recombinations.df$batch_collection_date) == F,]
print_info(nrow_last_step, sires_last_step, chr_last_step, "BOTH DATES PRESENT:")
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# Filter to only non-negative ages at time of collection
recombinations.df <- recombinations.df[recombinations.df$sire_age_at_collection >= 0,]
print_info(nrow_last_step, sires_last_step, chr_last_step, "NON-NEGATIVE AGE:")
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)
##----------------------------------------------------------------------------------------------------------------------------------




## FILTERS FOR CHROMOSOMES (e.g. informative markers)
##----------------------------------------------------------------------------------------------------------------------------------
# decide whether informative marker filtering will be done per-chromosome
if (length(settings) >= 36) {do_chromosome_filtering_this_time <- as.logical(settings[36])} else {do_chromosome_filtering_this_time <- FALSE}
if (do_chromosome_filtering_this_time == TRUE) {
    print_blank()
    fprint("OPTION: doing chromosome filtering (this will take a while...)")
    do_perchr_informative_markers <- as.logical(settings[8])
    fprint("OPTION: filter informative markers per-chromosome = {do_perchr_informative_markers}")

    if (do_perchr_informative_markers) {
        informative_marker_proportion <- as.numeric(settings[10])
        fprint("OPTION: filtering out lowest {informative_marker_proportion*100}% of informative markers per-chromosome.")
        if (informative_marker_proportion != 0) {
            marker_thresholds <- NULL
            for (i in 1:n_chromosomes) {
                col_name <- fstring("informative_markers_chr{i}")
                markers_this_chromosome <- sort(recombinations.df[,col_name])
                five_pc_point <- round(informative_marker_proportion * nrow(recombinations.df))
                five_pc_marker_num <- markers_this_chromosome[five_pc_point]
                marker_thresholds <- c(marker_thresholds, five_pc_marker_num)
            }
        } else {marker_thresholds <- rep(0,n_chromosomes)}
        fprint("thresholds calculated for lower {informative_marker_proportion*100}% of informative markers per chromosome")
        fprint("thresholds: {paste(marker_thresholds, collapse = ' ')}")
    }

    # load total genetic distance for each chromosome
    chromosome_distances <- read.table(fstring("files/chromosome_total_distances.txt"))[3]
    fprint("genetic distances per chromosome loaded")

    # remove any chromosomes with double crossovers below the minimum distance
    # remove any chromosomes in the lower 5% of informative markers for that chromosome number
    recombinations.df["kept_chromosomes"]         <- rep(NULL, nrow(recombinations.df))
    recombinations.df["kept_distance_M"]          <- rep(0, nrow(recombinations.df))
    recombinations.df["kept_recombinations"]      <- rep(0, nrow(recombinations.df))
    recombinations.df["recombination_rate"]       <- rep(0, nrow(recombinations.df))     ## crossovers per morgan
    recombinations.df["kept_informative_markers"] <- rep(0, nrow(recombinations.df))
    recombinations.df["kept_markers_adjusted"]    <- rep(0, nrow(recombinations.df))

    max_double_crossovers_below_threshold <- as.numeric(settings[6])
    fprint("OPTION: max double crossovers below threshold allowed = {max_double_crossovers_below_threshold}")
    print_blank()
    loop_start_time <- Sys.time()
    fprint("Processing chromosomes to keep, start time: {round(loop_start_time)}")
    for (i in 1:nrow(recombinations.df)) {
        # work out which chromosomes to keep
        kept_chromosomes <- NULL
        for (j in 1:n_chromosomes) {
            col_name1 <- fstring("num_below_threshold_chr{j}")
            all_distances_above_threshold <- recombinations.df[i,col_name1] <= max_double_crossovers_below_threshold
            if (do_perchr_informative_markers) {
                col_name2 <- fstring("informative_markers_chr{j}")
                not_lower_5pc_informative <- recombinations.df[i,col_name2] >= marker_thresholds[j]
                if (all_distances_above_threshold & not_lower_5pc_informative) {kept_chromosomes <- c(kept_chromosomes, j)}
            } else {if (all_distances_above_threshold) {kept_chromosomes <- c(kept_chromosomes, j)}}
        }
        # use kept chromosomes to calculate kept recombinations
        total_kept_recombinations <- 0
        total_kept_markers        <- 0
        for (j in kept_chromosomes) {
            col_name3 <- fstring("recombinations_chr{j}")
            total_kept_recombinations <- total_kept_recombinations + recombinations.df[i,col_name3]
            col_name4 <- fstring("informative_markers_chr{j}")
            total_kept_markers <- total_kept_markers + recombinations.df[i,col_name4]
        }
        # save kept chromosomes to dataframe
        recombinations.df[i,"kept_chromosomes"]         <- vector_to_string(kept_chromosomes)
        recombinations.df[i,"num_kept_chromosomes"]     <- length(kept_chromosomes)
        recombinations.df[i,"kept_distance_M"]          <- kept_distance <- sum(chromosome_distances[kept_chromosomes,]) / 100
        recombinations.df[i,"kept_recombinations"]      <- total_kept_recombinations
        recombinations.df[i,"recombination_rate"]       <- total_kept_recombinations / kept_distance
        recombinations.df[i,"kept_informative_markers"] <- total_kept_markers
        recombinations.df[i,"kept_markers_adjusted"]    <- total_kept_markers / kept_distance
        # print progress message
        if (i %% 5000 == 0) {fprint("chromosomes to keep processed for father-child pair {i} of {nrow(recombinations.df)} ({round(i/nrow(recombinations.df) * 100,2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")}
    }
    print_blank()
    
    total_chromosomes    <- n_chromosomes*nrow(recombinations.df)
    chromosomes_excluded <- total_chromosomes - sum(recombinations.df[,"num_kept_chromosomes"])
    fprint("TOTAL CHROMOSOMES EXCLUDED: {chromosomes_excluded} out of {total_chromosomes} ({round(chromosomes_excluded/total_chromosomes * 100,2)}%)")
    
    print_info(nrow_last_step, sires_last_step, chr_last_step, "CHROMOSOME EXCLUSIONS:")
    chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)
    
    
    if (do_perchr_informative_markers) {
        fprint("Chromosomes were excluded based on lowest {informative_marker_proportion*100}% per-chromsome informative markers.")
    } else {fprint("Chromosomes were not excluded based on informative markers per chromosome (per-chr filtering set to FALSE).")}


    # Filter to only individuals with adequately many chromosomes kept
    min_kept_chromosomes <- as.numeric(settings[12])
    fprint("OPTION: filtering to individuals with at least {min_kept_chromosomes} kept chromosomes.")
    recombinations.df <- recombinations.df[recombinations.df$num_kept_chromosomes >= min_kept_chromosomes,]
    print_info(nrow_last_step, sires_last_step, chr_last_step, fstring("AT LEAST {min_kept_chromosomes} KEPT CHROMOSOMES:"))
    nrow_last_step  <- nrow(recombinations.df)
    sires_last_step <- length(unique(recombinations.df$parent_id))
    chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)

    # If we did not filter informative markers per-chromosome, filter per-individual now
    do_perchr_informative_markers <- as.logical(settings[8])
    fprint("OPTION: filter informative markers per-chromosome = {do_perchr_informative_markers}")
    if (do_perchr_informative_markers == F) {
        informative_marker_proportion <- as.numeric(settings[10])
        fprint("Therefore, filtering out lowest {100*informative_marker_proportion}% informative markers per-individual.")
        recombinations.df <- recombinations.df[order(recombinations.df$total_informative_markers),]
        recombinations.df <- recombinations.df[round(informative_marker_proportion * nrow(recombinations.df)):nrow(recombinations.df),]
        recombinations.df <- recombinations.df[order(as.numeric(row.names(recombinations.df))),]
        print_info(nrow_last_step, sires_last_step, chr_last_step, fstring("TOP {100*(1-informative_marker_proportion)}% INFORMATIVE MARKERS:"))
    } else {fprint("Therefore, not filtering informative markers per-individual. {nrow(recombinations.df)} father-child pairs remain.")}
    nrow_last_step  <- nrow(recombinations.df)
    sires_last_step <- length(unique(recombinations.df$parent_id))
    chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)

    # save progress
    write.csv(recombinations.df, "files/recombinations_processing.csv", row.names = F)
    fprint("Progress saved in recombinations_processing.csv")

} else {
    fprint("OPTION: chromosome filtering skipped, reading data from recombinations_processing.csv instead")
    recombinations.df <- read.csv("files/recombinations_processing.csv")
}
##----------------------------------------------------------------------------------------------------------------------------------




## FILTERS FOR EXCESSIVE RECOMBINATION
##----------------------------------------------------------------------------------------------------------------------------------
# find max total recombinations allowed
if (length(settings) >= 18) {
    max_total_recombinations <- as.numeric(settings[20])
    max_recombinations_type  <- as.character(settings[18])
} else {
    max_total_recombinations <- 90
    max_recombinations_type  <- "number"
}
# remove individuals above max total recombinations threshold
if (max_recombinations_type == "proportion") {
    total_recombinations_message <- fstring("LOWEST {100*max_total_recombinations}% TOTAL RECOMBINATIONS:")
    max_recombinations_point <- round(max_total_recombinations * nrow(recombinations.df))
    recombinations.df <- recombinations.df[order(recombinations.df$total_recombinations),]
    recombinations.df <- recombinations.df[1:max_recombinations_point,]
    recombinations.df <- recombinations.df[order(as.numeric(row.names(recombinations.df))),]
} else {
    total_recombinations_message <- fstring("MAXIMUM OF {max_total_recombinations} TOTAL RECOMBINATIONS:")
    recombinations.df <- recombinations.df[recombinations.df$total_recombinations <= max_total_recombinations,]
}
print_info(nrow_last_step, sires_last_step, chr_last_step, total_recombinations_message)
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# Filter out very high recombinations rates (i.e. outlier), if we put a setting for that
if (length(settings) >= 24) {
    max_recombination_rate <- as.numeric(settings[24])
    recombinations.df <- recombinations.df[recombinations.df$recombination_rate <= max_recombination_rate,]
    print_info(nrow_last_step, sires_last_step, chr_last_step, fstring("MAX OF {max_recombination_rate} RECOMBINATIONS PER MORGAN:"))
    nrow_last_step  <- nrow(recombinations.df)
    sires_last_step <- length(unique(recombinations.df$parent_id))
    chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)
}
##----------------------------------------------------------------------------------------------------------------------------------




## PER-SIRE FILTERING
##----------------------------------------------------------------------------------------------------------------------------------
# read settings pertaining to per-sire filtering
min_halfsib_size <- as.numeric(settings[14])
min_unique_ages  <- as.numeric(settings[16])
if (length(settings) >= 28) {min_breeding_window <- as.numeric(settings[28])} else {min_breeding_window <- 0}
if (length(settings) >= 30) {discrete_method     <- as.character(settings[30])} else {discrete_method <- "round"}
if (length(settings) >= 32) {min_kids_in_year    <- as.numeric(settings[32])} else {min_kids_in_year <- 0}
if (length(settings) >= 34) {min_years_with_kids <- as.numeric(settings[34])} else {min_years_with_kids <- 2}
if (discrete_method == "floor") {
    discrete_years <- function(expression) {floor(expression)}
} else {
    discrete_years <- function(expression) {round(expression)}
}

# make variables for: remaining halfsibs in data, rounded age of sire, children per rounded age, number of rounded ages with children above a chosen threshold, and total time between first child and last child (breeding window)
loop_start_time <- Sys.time()
print_blank()
fprint("Processing per-sire information, start time: {round(loop_start_time)}")
sires <- unique(recombinations.df$parent_id)
i <- 1
recombinations.df$sire_age_years <- as.factor(discrete_years(recombinations.df$sire_age_at_collection / 365))
for (bull in sires) {
    his_offspring.df <- recombinations.df[recombinations.df$parent_id == bull,]
    remaining_halfsibs <- nrow(his_offspring.df)
    recombinations.df[recombinations.df$parent_id == bull,"halfsibs_in_data"] <- remaining_halfsibs
    
    his_ages <- unique(his_offspring.df$sire_age_at_collection)
    recombinations.df[recombinations.df$parent_id == bull,"sire_unique_ages"] <- length(his_ages)
    breeding_window <- max(his_ages) - min(his_ages)
    recombinations.df[recombinations.df$parent_id == bull,"sire_breeding_window"] <- breeding_window
    
    his_years <- discrete_years(his_offspring.df$sire_age_at_collection / 365)
    his_distinct_years <- sort(unique(his_years))
    recombinations.df[recombinations.df$parent_id == bull,"sire_rounded_ages"] <- vector_to_string(his_distinct_years)
    
    kids_each_year <- NULL
    for (age in his_distinct_years) {
        kids_that_year <- sum(his_years == age)
        kids_each_year <- c(kids_each_year, kids_that_year)
    }
    sire_ages_with_enough_kids <- sum(kids_each_year >= min_kids_in_year)
    recombinations.df[recombinations.df$parent_id == bull,"sire_children_per_age"] <- vector_to_string(kids_each_year)
    recombinations.df[recombinations.df$parent_id == bull,"sire_years_with_enough_children"] <- sire_ages_with_enough_kids
    
    if (i %% 500 == 0) {fprint("Per-sire information processed for sire {i} of {length(sires)} ({round(i/length(sires) * 100,2)}% complete). Time elapsed: {get_time_diff_min(loop_start_time, Sys.time())} min")}
    i = i+1
}
print_blank()


# remove halfsib families smaller than threshold
recombinations.df <- recombinations.df[recombinations.df$halfsibs_in_data >= min_halfsib_size,]
print_info(nrow_last_step, sires_last_step, chr_last_step, fstring("HALFSIB FAMILIES OF AT LEAST {min_halfsib_size}:"))
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)


# remove sires with too narrow a breeding window
recombinations.df <- recombinations.df[recombinations.df$sire_breeding_window >= min_breeding_window,]
print_info(nrow_last_step, sires_last_step, chr_last_step, fstring("SIRE BREEDING WINDOW OF AT LEAST {min_breeding_window} DAYS:"))
nrow_last_step  <- nrow(recombinations.df)
sires_last_step <- length(unique(recombinations.df$parent_id))
chr_last_step   <- sum(recombinations.df$num_kept_chromosomes)
##----------------------------------------------------------------------------------------------------------------------------------




# write output to csv - written to current directory, not "Cows"
print_blank()
fprint("Filtering complete. {nrow(recombinations.df)} father-child pairs, {length(unique(recombinations.df$parent_id))} sires, and {sum(recombinations.df$num_kept_chromosomes)} chromosomes remain.")
write.csv(recombinations.df, "files/recombinations_filtered.csv", row.names = F)
print_end_message(start_time, F)