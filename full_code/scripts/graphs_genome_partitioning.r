# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
source("scripts/functions.r")



## Reorganise data
partition_results <- read.csv("files/genome_partition_results.csv")
category <- rep(c("(a) Mean, filter 0", "(b) Mean, filter 365", "(c) Slope, filter 0", "(d) Slope, filter 365"), each = 29)
values   <- c(partition_results$mean_ldak_h2, partition_results$slope_ldak_h2)
chrs     <- rep(1:29, 4)
bp       <- rep(partition_results$bp, 2)
cM       <- rep(partition_results$cM, 2)
partition_recoded <- data.frame(chromosome = chrs, bp, cM, category, h2 = values)



## Make genome partitioning plot
colour <- "darkblue"
ggplot(partition_recoded, aes(x = bp / 10^6, y = h2)) +
    stat_smooth(geom = "line", method = "lm", se = F, col = "orange", alpha = 0.8) +
    geom_point(size = 3.5, col = colour, alpha = 0.9) +
    geom_text(aes(label = chromosome), col = "white", size = 2) +
    xlab("Chromosome size (Mbp between first and last SNP)") + ylab("Estimated h^2 per chromosome") +
    theme_bw() +
    facet_wrap("category", nrow = 2, ncol = 2)
ggsave("output/genome_partitioning.png", dpi = 600, width = 8, height = 8)