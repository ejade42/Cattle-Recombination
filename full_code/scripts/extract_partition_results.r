## runs in main directory, accessing the genome_partitioning directory, to extract the results of the per-chromosome heritability partitioning analysis

## LOAD FUNCTIONS
source("scripts/functions.r")    ## load common functions

extract_values_independent <- function(results) {## function takes output of readLines(), and extracts the h2, se, and p value from the .hsq file
    line <- strsplit(results[5], "\t", fixed = T)
    h2 <- as.numeric(line[[1]][2])
    se <- as.numeric(line[[1]][3])

    line <- strsplit(results[10], "\t", fixed = T)
    p <- as.numeric(line[[1]][2])
    c(h2, se, p)
}

extract_values_together <- function(results, chr) {
    line <- strsplit(results[chr+32], "\t", fixed = T)
    h2 <- as.numeric(line[[1]][2])
    se <- as.numeric(line[[1]][3])
    c(h2, se)
}

extract_values_ldak <- function(results, chr) {
    line <- strsplit(results[chr+17], " ", fixed = T)
    h2 <- as.numeric(line[[1]][2])
    sd <- as.numeric(line[[1]][3])
    c(h2, sd)
}



## EXTRACT GCTA RESULTS FOR EACH CHROMOSOME
genome_partition_results.df <- NULL
for (filter in c(0, 365)) {
    for (chr in 1:29) {
        ## READ IN CHROMOSOME INFORMATION FROM MAP FILE
        map <- read.table(fstring("qc/chr_{chr}/chr_{chr}_01.map"))
        bp <- map[nrow(map),4] - map[1,4]
        cM <- map[nrow(map),3] - map[1,3]
        
        ## READ IN RESULTS FROM LDAK ANALYSIS (ALSO FITTED TO ALL CHROMOSOMES SIMULTANEOUSLY)
        mean_ldak_results  <- readLines(fstring("genome_partitioning/mean_results_filter_{filter}_allchrs.reml"))
        slope_ldak_results <- readLines(fstring("genome_partitioning/slope_results_filter_{filter}_allchrs.reml"))
        mean_ldak <- extract_values_ldak(mean_ldak_results, chr)[1]
        slope_ldak <- extract_values_ldak(slope_ldak_results, chr)[1]

        ## COMBINE AND ADD NEW ROW TO DATA
        row <- c(filter, chr, bp, cM, mean_ldak, slope_ldak)
        genome_partition_results.df <- rbind(genome_partition_results.df, row)
    }
}



## WRITE RESULTS TO FILE
colnames(genome_partition_results.df) <- c("filter", "chromosome", "bp", "cM", "mean_ldak_h2", "slope_ldak_h2")
write.csv(genome_partition_results.df, "files/genome_partition_results.csv", row.names = F)