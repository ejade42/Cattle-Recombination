# Runs in qc directory, integrates genetic distances from Shen linkage map into the LIC data


## Set up
start_time <- Sys.time()
source("../scripts/functions.r")   ## load common functions
fprint("Start time: {round(start_time)}")


## Read in data, create "combined" variable for each dataset
## "combined" gives a unique ID to each SNP based on both chromosome and position
our_data <- read.table("all_cows_01_map.map")
colnames(our_data) <- c("chromosome", "snp_id", "genetic_distance_cM", "base_position")
linkage_map <- read.table("../input/Shen_linkage_maps.rmap", header=T)
linkage_map$combined <- linkage_map$Chr * 10^9 + linkage_map$Location
our_data$combined    <- our_data$chromosome * 10^9 + our_data$base_position

## Get average SNP-to-SNP distance between Jersey and Holstein males (in Morgans)
linkage_map$cow_avg <- (linkage_map$HO_m + linkage_map$JE_m) / 2


## Get cumulative distance (in cM) by adding these up along each chromosome, x100 for M -> cM
linkage_map$cumulative_distance_cM <- rep(0, nrow(linkage_map))
distance_M <- 0
chromosome <- 1
for (i in 1:nrow(linkage_map)) {
    if (chromosome == linkage_map[i,"Chr"]) {
        distance_M <- distance_M + linkage_map[i,"cow_avg"]
    } else {
        chromosome <- linkage_map[i,"Chr"]
        distance_M <- 0
    }
    linkage_map[i,"cumulative_distance_cM"] <- distance_M * 100
    if (i %% 5000 == 0) {
        fprint("distance calculated for marker {i} of {nrow(linkage_map)} ({round(i/nrow(linkage_map) * 100,2)}% complete)")
    }
}

## Merge datasets, copy cumulative distance into Column 3 for .map format
## Filter down to just the four columns for .map format, and remove any SNPs with missing IDs
new_data <- merge(our_data, linkage_map, all.x = F)
new_data$genetic_distance_cM <- new_data$cumulative_distance_cM
new_data <- new_data[,2:5]
new_data <- new_data[new_data$snp_id != ".",]

## Write .map file, and create a new file with list of all SNP IDs in final .map file
write.table(new_data, "all_cows_01_edited.map", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(new_data$snp_id, "all_cows_01_snp_ids_to_keep", row.names = F, col.names = F, quote = F)
print_end_message(start_time)