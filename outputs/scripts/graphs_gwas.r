# Load in packages/functions
start_time <- Sys.time()
library(tidyverse)  ## v2.0.0
library(topr)       ## v2.0.0
library(egg)        ## v0.4.5
source("scripts/functions.r")
print_blank()
fprint("Start time: {round(start_time)}")


## to be used for all annotations:
linewidth  <- 0.7
linecolour <- "orange"
linealpha  <- 0.9

## to be used for all manhattan plots:
sign_thresh_label_size <- 0
hide_chrticks_from     <- 23


## MEAN, FILTER 0 - FIGURE 5
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fprint("WORKING ON PHENOTYPE: MEAN, FILTER: 0 (FIGURE 5)")
gwas_all <- read.table("input/GWAS_mean_filter_0_breed_all.mlma", header = T)
gwas_hol <- read.table("input/GWAS_mean_filter_0_breed_holstein.mlma", header = T)
gwas_jer <- read.table("input/GWAS_mean_filter_0_breed_jersey.mlma", header = T)

confidence <- 0.05
n_snps_all <- nrow(gwas_all)
n_snps_hol <- nrow(gwas_hol)
n_snps_jer <- nrow(gwas_jer)
bonferroni_all <- confidence / n_snps_all
bonferroni_hol <- confidence / n_snps_hol
bonferroni_jer <- confidence / n_snps_jer
fprint("{n_snps_all} SNPs remain in the data for all cows, Bonferroni threshold = {signif(bonferroni_all, 4)}, equivalent to score of {round(-log10(bonferroni_all), 3)}")
fprint("{n_snps_hol} SNPs remain in the data for Holstein, Bonferroni threshold = {signif(bonferroni_hol, 4)}, equivalent to score of {round(-log10(bonferroni_hol), 3)}")
fprint("{n_snps_jer} SNPs remain in the data for Jersey, Bonferroni threshold = {signif(bonferroni_jer, 4)}, equivalent to score of {round(-log10(bonferroni_jer), 3)}")


colnames(gwas_all)[c(1,3)] <- c("chrom", "pos")
all <- topr::manhattan(gwas_all, sign_thresh = signif(bonferroni_all, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    geom_segment(aes(x = 1.146*10^9, xend = 1.146*10^9, y = -Inf, yend = 14.0), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## first chr 10
    geom_segment(aes(x = 1.211*10^9, xend = 1.211*10^9, y = -Inf, yend = 12.2), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## second chr 10
    ggtitle("(a) All cattle, n = 1976")
hol <- topr::manhattan(gwas_hol, sign_thresh = signif(bonferroni_hol, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    geom_segment(aes(x = 1.146*10^9, xend = 1.146*10^9, y = -Inf, yend = 7.45), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## first chr 10
    geom_segment(aes(x = 1.211*10^9, xend = 1.211*10^9, y = -Inf, yend = 7.66), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## second chr 10
    geom_segment(aes(x = 1.914*10^9, xend = 1.912*10^9, y = -Inf, yend = 5.47), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 18
    ggtitle("(b) Near-pure Holstein only, n = 895")
jer <- topr::manhattan(gwas_jer, sign_thresh = signif(bonferroni_jer, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    geom_segment(aes(x = 1.146*10^9, xend = 1.146*10^9, y = -Inf, yend = 6.25), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## first chr 10
    geom_segment(aes(x = 1.211*10^9, xend = 1.211*10^9, y = -Inf, yend = 3.75), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## second chr 10
    geom_segment(aes(x = 0.907*10^9, xend = 0.907*10^9, y = -Inf, yend = 5.60), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 8
    geom_segment(aes(x = 0.350*10^9, xend = 0.350*10^9, y = -Inf, yend = 4.72), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 3
    ggtitle("(c) Near-pure Jersey only, n = 572")
qqall <- topr::qqtopr(gwas_all, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
        scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqhol <- topr::qqtopr(gwas_hol, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqjer <- topr::qqtopr(gwas_jer, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))


plot <- ggarrange(all, qqall, hol, qqhol, jer, qqjer, nrow = 3, ncol = 2,
         widths = c(15,4), heights = rep(4,3))
ggsave("output/GWAS_filter0_mean.png", dpi = 600, plot = plot, width = 8.5, height = 7)
for (i in 1:5) {print_blank()}
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







## SLOPE, FILTER 365 - FIGURE 6
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fprint("WORKING ON PHENOTYPE: SLOPE, FILTER: 365 (FIGURE 6)")
gwas_all <- read.table("input/GWAS_slope_filter_365_breed_all.mlma", header = T)
gwas_hol <- read.table("input/GWAS_slope_filter_365_breed_holstein.mlma", header = T)
gwas_jer <- read.table("input/GWAS_slope_filter_365_breed_jersey.mlma", header = T)

confidence <- 0.05
n_snps_all <- nrow(gwas_all)
n_snps_hol <- nrow(gwas_hol)
n_snps_jer <- nrow(gwas_jer)
bonferroni_all <- confidence / n_snps_all
bonferroni_hol <- confidence / n_snps_hol
bonferroni_jer <- confidence / n_snps_jer
fprint("{n_snps_all} SNPs remain in the data for all cows, Bonferroni threshold = {signif(bonferroni_all, 4)}, equivalent to score of {round(-log10(bonferroni_all), 3)}")
fprint("{n_snps_hol} SNPs remain in the data for Holstein, Bonferroni threshold = {signif(bonferroni_hol, 4)}, equivalent to score of {round(-log10(bonferroni_hol), 3)}")
fprint("{n_snps_jer} SNPs remain in the data for Jersey, Bonferroni threshold = {signif(bonferroni_jer, 4)}, equivalent to score of {round(-log10(bonferroni_jer), 3)}")


colnames(gwas_all)[c(1,3)] <- c("chrom", "pos")
all <- topr::manhattan(gwas_all, sign_thresh = signif(bonferroni_all, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(a) All cattle, n = 444")
hol <- topr::manhattan(gwas_hol, sign_thresh = signif(bonferroni_hol, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    geom_segment(aes(x = 0.280*10^9, xend = 0.280*10^9, y = -Inf, yend = 4.33), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 2
    geom_segment(aes(x = 1.780*10^9, xend = 1.780*10^9, y = -Inf, yend = 4.15), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 17
    ggtitle("(b) Near-pure Holstein only, n = 203")
jer <- topr::manhattan(gwas_jer, sign_thresh = signif(bonferroni_jer, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    geom_segment(aes(x = 29524658, xend = 29524658, y = -Inf, yend = 4.84), col = linecolour, alpha = linealpha, linewidth = linewidth) +  ## chr 1
    ggtitle("(c) Near-pure Jersey only, n = 132")
qqall <- topr::qqtopr(gwas_all, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqhol <- topr::qqtopr(gwas_hol, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7), breaks = seq(0,8,2)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqjer <- topr::qqtopr(gwas_jer, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))


plot <- ggarrange(all, qqall, hol, qqhol, jer, qqjer, nrow = 3, ncol = 2,
         widths = c(15,4), heights = rep(4,3))
ggsave("output/GWAS_filter365_slope.png", dpi = 600, plot = plot, width = 8.5, height = 7)
for (i in 1:5) {print_blank()}
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







## MEAN, FILTER 365 - SUPPLEMENTARY FIGURE S3
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fprint("WORKING ON PHENOTYPE: MEAN, FILTER: 365 (SUPPLEMENTARY FIGURE S3)")
gwas_all <- read.table("input/GWAS_mean_filter_365_breed_all.mlma", header = T)
gwas_hol <- read.table("input/GWAS_mean_filter_365_breed_holstein.mlma", header = T)
gwas_jer <- read.table("input/GWAS_mean_filter_365_breed_jersey.mlma", header = T)

confidence <- 0.05
n_snps_all <- nrow(gwas_all)
n_snps_hol <- nrow(gwas_hol)
n_snps_jer <- nrow(gwas_jer)
bonferroni_all <- confidence / n_snps_all
bonferroni_hol <- confidence / n_snps_hol
bonferroni_jer <- confidence / n_snps_jer
fprint("{n_snps_all} SNPs remain in the data for all cows, Bonferroni threshold = {signif(bonferroni_all, 4)}, equivalent to score of {round(-log10(bonferroni_all), 3)}")
fprint("{n_snps_hol} SNPs remain in the data for Holstein, Bonferroni threshold = {signif(bonferroni_hol, 4)}, equivalent to score of {round(-log10(bonferroni_hol), 3)}")
fprint("{n_snps_jer} SNPs remain in the data for Jersey, Bonferroni threshold = {signif(bonferroni_jer, 4)}, equivalent to score of {round(-log10(bonferroni_jer), 3)}")


colnames(gwas_all)[c(1,3)] <- c("chrom", "pos")
all <- topr::manhattan(gwas_all, sign_thresh = signif(bonferroni_all, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(a) All cattle, n = 444")
hol <- topr::manhattan(gwas_hol, sign_thresh = signif(bonferroni_hol, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(b) Near-pure Holstein only, n = 203")
jer <- topr::manhattan(gwas_jer, sign_thresh = signif(bonferroni_jer, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = hide_chrticks_from, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(c) Near-pure Jersey only, n = 132")
qqall <- topr::qqtopr(gwas_all, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
        scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqhol <- topr::qqtopr(gwas_hol, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqjer <- topr::qqtopr(gwas_jer, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))


plot <- ggarrange(all, qqall, hol, qqhol, jer, qqjer, nrow = 3, ncol = 2,
         widths = c(15,4), heights = rep(4,3))
ggsave("output/GWAS_filter365_mean.png", dpi = 600, plot = plot, width = 8.5, height = 7)
for (i in 1:5) {print_blank()}
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







## SLOPE, FILTER 0 - SUPPLEMENTARY FIGURE S4
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fprint("WORKING ON PHENOTYPE: SLOPE, FILTER: 0 (SUPPLEMENTARY FIGURE S4)")
gwas_all <- read.table("input/GWAS_slope_filter_0_breed_all.mlma", header = T)
gwas_hol <- read.table("input/GWAS_slope_filter_0_breed_holstein.mlma", header = T)
gwas_jer <- read.table("input/GWAS_slope_filter_0_breed_jersey.mlma", header = T)

confidence <- 0.05
n_snps_all <- nrow(gwas_all)
n_snps_hol <- nrow(gwas_hol)
n_snps_jer <- nrow(gwas_jer)
bonferroni_all <- confidence / n_snps_all
bonferroni_hol <- confidence / n_snps_hol
bonferroni_jer <- confidence / n_snps_jer
fprint("{n_snps_all} SNPs remain in the data for all cows, Bonferroni threshold = {signif(bonferroni_all, 4)}, equivalent to score of {round(-log10(bonferroni_all), 3)}")
fprint("{n_snps_hol} SNPs remain in the data for Holstein, Bonferroni threshold = {signif(bonferroni_hol, 4)}, equivalent to score of {round(-log10(bonferroni_hol), 3)}")
fprint("{n_snps_jer} SNPs remain in the data for Jersey, Bonferroni threshold = {signif(bonferroni_jer, 4)}, equivalent to score of {round(-log10(bonferroni_jer), 3)}")


colnames(gwas_all)[c(1,3)] <- c("chrom", "pos")
all <- topr::manhattan(gwas_all, sign_thresh = signif(bonferroni_all, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = 23, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(a) All cattle, n = 1963")
hol <- topr::manhattan(gwas_hol, sign_thresh = signif(bonferroni_hol, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = 23, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(b) Near-pure Holstein only, n = 889")
jer <- topr::manhattan(gwas_jer, sign_thresh = signif(bonferroni_jer, 4), get_chr_lengths_from_data = T, ymax = 15, hide_chrticks_from = 23, sign_thresh_label_size = sign_thresh_label_size) + 
    scale_y_continuous(breaks = seq(0,14,2)) +
    theme_bw() + theme(panel.grid.major.x = element_blank(), 
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_line(linewidth=0.125),
                       panel.grid.minor.y = element_blank()) +
    ggtitle("(c) Near-pure Jersey only, n = 570")
qqall <- topr::qqtopr(gwas_all, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqhol <- topr::qqtopr(gwas_hol, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7), breaks = seq(0,8,2)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))
qqjer <- topr::qqtopr(gwas_jer, diagonal_line_color = "red") +
    scale_x_continuous(expand = c(0, 0), limits=c(0,5)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,7)) +
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.x = element_line(linewidth=0.2),
                       panel.grid.major.y = element_line(linewidth=0.2))


plot <- ggarrange(all, qqall, hol, qqhol, jer, qqjer, nrow = 3, ncol = 2,
         widths = c(15,4), heights = rep(4,3))
ggsave("output/GWAS_filter0_slope.png", dpi = 600, plot = plot, width = 8.5, height = 7)
for (i in 1:5) {print_blank()}
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







## GWAS RESULTS TABLE - ALL SNPS WITH P < 10^-4, FOR MEAN_0 AND SLOPE_365 PHENOTYPES - SUPPLEMENTARY TABLE S3
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fprint("EXTRACTING SNP-WISE RESULTS, FOR ALL SNPS WTIH P < 10^-4 (SUPPLEMENTARY TABLE S3)")
fprint("Reading in data")
mean_all <- read.table("input/GWAS_mean_filter_0_breed_all.mlma", header = T)
mean_hol <- read.table("input/GWAS_mean_filter_0_breed_holstein.mlma", header = T)
mean_jer <- read.table("input/GWAS_mean_filter_0_breed_jersey.mlma", header = T)
slope_all <- read.table("input/GWAS_slope_filter_365_breed_all.mlma", header = T)
slope_hol <- read.table("input/GWAS_slope_filter_365_breed_holstein.mlma", header = T)
slope_jer <- read.table("input/GWAS_slope_filter_365_breed_jersey.mlma", header = T)

fprint("Rounding p values")
mean_all$p <- signif(mean_all$p, 4)
mean_hol$p <- signif(mean_hol$p, 4)
mean_jer$p <- signif(mean_jer$p, 4)
slope_all$p <- signif(slope_all$p, 4)
slope_hol$p <- signif(slope_hol$p, 4)
slope_jer$p <- signif(slope_jer$p, 4)

fprint("Rounding -log10(p) scores")
mean_all$log_score <- round(-log10(mean_all$p), 3)
mean_hol$log_score <- round(-log10(mean_hol$p), 3)
mean_jer$log_score <- round(-log10(mean_jer$p), 3)
slope_all$log_score <- round(-log10(slope_all$p), 3)
slope_hol$log_score <- round(-log10(slope_hol$p), 3)
slope_jer$log_score <- round(-log10(slope_jer$p), 3)

fprint("Rounding allele frequencies")
mean_all$Freq <- round(mean_all$Freq, 3)
mean_hol$Freq <- round(mean_hol$Freq, 3)
mean_jer$Freq <- round(mean_jer$Freq, 3)
slope_all$Freq <- round(slope_all$Freq, 3)
slope_hol$Freq <- round(slope_hol$Freq, 3)
slope_jer$Freq <- round(slope_jer$Freq, 3)

fprint("Selecting *minor* allele frequencies")
for (i in 1:nrow(mean_all)) {
    mean_all[i,"Freq"] <- min(mean_all[i,"Freq"], 1 - mean_all[i,"Freq"])
}
for (i in 1:nrow(mean_hol)) {
    mean_hol[i,"Freq"] <- min(mean_hol[i,"Freq"], 1 - mean_hol[i,"Freq"])
}
for (i in 1:nrow(mean_jer)) {
    mean_jer[i,"Freq"] <- min(mean_jer[i,"Freq"], 1 - mean_jer[i,"Freq"])
}
for (i in 1:nrow(slope_all)) {
    slope_all[i,"Freq"] <- min(slope_all[i,"Freq"], 1 - slope_all[i,"Freq"])
}
for (i in 1:nrow(slope_hol)) {
    slope_hol[i,"Freq"] <- min(slope_hol[i,"Freq"], 1 - slope_hol[i,"Freq"])
}
for (i in 1:nrow(slope_jer)) {
    slope_jer[i,"Freq"] <- min(slope_jer[i,"Freq"], 1 - slope_jer[i,"Freq"])
}

fprint("Adding phenotype and breed columns")
mean_all$phenotype  <- mean_hol$phenotype  <- mean_jer$phenotype  <- "Mean"
slope_all$phenotype <- slope_hol$phenotype <- slope_jer$phenotype <- "Slope"
mean_all$breed <- slope_all$breed <- "All"
mean_hol$breed <- slope_hol$breed <- "Holstein"
mean_jer$breed <- slope_jer$breed <- "Jersey"

fprint("Subsetting to only useful columns")
mean_all <- mean_all %>% select(phenotype, breed, SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)
mean_hol <- mean_hol %>% select(phenotype, breed,  SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)
mean_jer <- mean_jer %>% select(phenotype, breed, SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)
slope_all <- slope_all %>% select(phenotype, breed, SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)
slope_hol <- slope_hol %>% select(phenotype, breed, SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)
slope_jer <- slope_jer %>% select(phenotype, breed, SNP, Chr, bp, Freq, p, log_score) %>% filter(., p < 10^-4)

fprint("Combining data")
gwas_results <- rbind(mean_all, mean_hol, mean_jer, slope_all, slope_hol, slope_jer)
colnames(gwas_results) <- c("phenotype", "breed", "SNP", "chromosome", "bp", "minor_allele_frequency", "p", "log_score")

fprint("Saving output")
write.csv(gwas_results, "output/gwas.csv", row.names = F)
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


print_blank()
print_end_message(start_time)
