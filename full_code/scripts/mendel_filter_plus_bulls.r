## Runs in qc
## Set up
start_time <- Sys.time()
source("../scripts/functions.r")   ## load common functions
fprint("Start time: {round(start_time)}")


## Load in fam files before and after Mendel filtering
cows_before_mendel.df <- read.table("all_cows_07.fam")[,1:2]
colnames(cows_before_mendel.df) <- c("family_id", "individual_id")
fprint("Before Mendel filter, we have {nrow(cows_before_mendel.df)} cattle.")

cows_after_mendel.df <- read.table("all_cows_09_after_mendel.fam")[,1:2]
colnames(cows_after_mendel.df) <- c("family_id", "individual_id")
fprint("After Mendel filter, we have {nrow(cows_after_mendel.df)} cattle.")



## Load in sires
file_to_use <- "bull_breed_mix.csv"
#file_to_use <- "sire_(3151_animals_of_interest)_locations.csv"
sires.df <- read.csv(fstring("../input/{file_to_use}"))
sires.df <- sires.df[,c("lic_animal_key", colnames(sires.df)[2])]
colnames(sires.df) <- c("long_id", "total")

plink_ids.df <- read.table("../files/newID4plink_01")[,c(4,1)]
colnames(plink_ids.df) <- c("individual_id", "long_id")

merged_sires.df <- merge(sires.df, plink_ids.df)
merged_sires.df$total <- NULL

fprint("{nrow(merged_sires.df)} sires loaded from file.")



## Subset sires to all those present before Mendel filter, and to all those present after Mendel filter
merged_sires_pre_filter.df <- merged_sires.df[merged_sires.df$individual_id %in% cows_before_mendel.df$individual_id,]
fprint("{nrow(merged_sires_pre_filter.df)} sires present in pre-Mendel PLINK data.")

merged_sires_post_filter.df <- merged_sires.df[merged_sires.df$individual_id %in% cows_after_mendel.df$individual_id,]
fprint("{nrow(merged_sires_post_filter.df)} sires present in post-Mendel PLINK data.")



## Keep cattle if they are one of our 3151 sires, or if they passed the Mendel filter.
cattle_to_keep.df <- cows_before_mendel.df[cows_before_mendel.df$individual_id %in% cows_after_mendel.df$individual_id |
                                           cows_before_mendel.df$individual_id %in% merged_sires_pre_filter.df$individual_id,]
fprint("Preserving all sires (even if excluded by Mendel filter), we will keep {nrow(cattle_to_keep.df)} total cattle.")



## Write list of cattle to keep to file
filename <- "all_cows_09_cattle_ids_to_keep"
filename2 <- "../files/sires_of_interest"
write.table(cattle_to_keep.df, filename, row.names = F, col.names = F, quote = F)
write.table(merged_sires.df[,c(2,1)], filename2, row.names = F, col.names = F, quote = F)
fprint("IDs of {nrow(cattle_to_keep.df)} cattle written to {filename}")
print(paste("IDs of", nrow(merged_sires.df), "cattle written to", filename2), quote = F)
print_end_message(start_time)
