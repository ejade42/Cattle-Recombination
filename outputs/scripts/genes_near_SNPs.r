## Runs in main directory, uses GWAS output and list of genes (input file) to find genes near significant SNPs


## Process SNPs from GWAS results
## -------------------------------------------------------------------------------------------
## Get significant SNP list

hol.list <- scan("files/gwas_file_holstein.list", what="character")
all.list <- scan("files/gwas_file_all.list", what="character")
jer.list <- scan("files/gwas_file_jersey.list", what="character")
# 4


## Process all cows

ngldp.cut <- 4

sig.dat <- NULL
for (i in all.list){
  skip_to_next <- FALSE
  tryCatch({
    dd <- read.table(paste0("input/", i), header = TRUE)
    dd$ngldp <- -log10(dd$p)
    sig <- subset(dd, dd$ngldp > ngldp.cut)
    f.nm <-  strsplit(i, "_")
    sig.nm <-  data.frame(sig, do.call(rbind, f.nm))
    sig1 <- sig.nm[,-c(11,13,15)]
    sig.dat  <- rbind(sig.dat, sig1)
  }, error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
}

names(sig.dat)[11:13] <- c("param", "days", "breed")

write.csv(sig.dat, file = "files/nlgp4_SNP_AWS_all.csv", quote = FALSE, row.names = FALSE)

dim(sig.dat) #20


## Process Holstein

ngldp.cut <- 4
  
sig.dat <- NULL
for (i in hol.list){
skip_to_next <- FALSE
tryCatch({
  dd <- read.table(paste0("input/", i), header = TRUE)
  dd$ngldp <- -log10(dd$p)
  sig <- subset(dd, dd$ngldp > ngldp.cut)
  f.nm <-  strsplit(i, "_")
  sig.nm <-  data.frame(sig, do.call(rbind, f.nm))
  sig1 <- sig.nm[,-c(11,13,15)]
  sig.dat  <- rbind(sig.dat, sig1)
}, error = function(e) { skip_to_next <<- TRUE})
if(skip_to_next) { next }
}

names(sig.dat)[11:13] <- c("param", "days", "breed")

write.csv(sig.dat, file = "files/nlgp4_SNP_AWS_holstein.csv", quote = FALSE, row.names = FALSE)
dim(sig.dat) #24


## Process Jersey

ngldp.cut <- 4

sig.dat <- NULL
for (i in jer.list){
  skip_to_next <- FALSE
  tryCatch({
    dd <- read.table(paste0("input/", i), header = TRUE)
    dd$ngldp <- -log10(dd$p)
    sig <- subset(dd, dd$ngldp > ngldp.cut)
    f.nm <-  strsplit(i, "_")
    sig.nm <-  data.frame(sig, do.call(rbind, f.nm))
    sig1 <- sig.nm[,-c(11,13,15)]
    sig.dat  <- rbind(sig.dat, sig1)
  }, error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
}

names(sig.dat)[11:13] <- c("param", "days", "breed")

write.csv(sig.dat, file = "files/nlgp4_SNP_AWS_jersey.csv", quote = FALSE, row.names = FALSE)

dim(sig.dat) #23



## Merge data

a1 <- read.csv("files/nlgp4_SNP_AWS_all.csv")
h1 <- read.csv("files/nlgp4_SNP_AWS_holstein.csv")
j1 <- read.csv("files/nlgp4_SNP_AWS_jersey.csv")
dim(h1)
dim(a1)
b1 <- rbind(h1,a1)
b1 <- rbind(b1,j1)
dim(b1) #67
write.csv(b1, file = "files/nlgp4_SNP_AWS.csv", quote = FALSE, row.names = FALSE)
## -------------------------------------------------------------------------------------------



# Cattle gwas results
# Need ARS-UCD1.2/bosTau9 
## -------------------------------------------------------------------------------------------

## Extract genes in significant intervals 

## Read gene name file
## comes from https://www.ncbi.nlm.nih.gov/gene/?term=ARS-UCD1.2/bosTau9
## this file was edited so it oculd be easily read in R
genes <- read.csv("input/gene_names_pos_forR.csv")
head(genes)
colnames(genes)[1] <- "tax_id"
dim(genes)# 31660    12


## Read SNP file

snp <- read.csv("files/nlgp4_SNP_AWS.csv")
dim(snp) # 67
head(snp)

newdat <- NULL
for (j in 1:length(snp[,1])) {
  skip_to_next <- FALSE
  tryCatch({
      dd  <- snp[j, ]
      chr <- dd$Chr
      gene.chr <- subset(genes, genes$chromosome==chr)
      row.names(gene.chr) <- seq(1,length(gene.chr$tax_id))
      ggb  <- gene.chr[(abs(dd$bp-gene.chr$start_position_on_the_genomic_accession))<250000,] 
      gggb <- gene.chr[(abs(dd$bp-gene.chr$end_position_on_the_genomic_accession))<250000,] 
      ggbb <- rbind(ggb, gggb)
      ggbb$SNP.param <- dd$param
      ggbb$SNP.breed <- dd$breed
      ggbb$SNP.days  <- dd$days
      ggbb$SNP       <- dd$SNP
      ggbb$SNP.pos   <- dd$bp
      ggbb$SNP.a1    <- dd$a1
      ggbb$SNP.a2    <- dd$a2  
      ggbb$SNP.freq  <- dd$Freq
      ggbb$SNP.b     <- dd$b
      ggbb$SNP.se    <- dd$se
      ggbb$SNP.p     <- dd$p
      ggbb$SNP.ngldp <- dd$ngldp
      newdat <- rbind(newdat, ggbb)
  }, error = function(e) {skip_to_next <<- TRUE})
  if (skip_to_next) {next}
}

dim(newdat) #2053

write.csv(newdat, file = "output/SNP_gene_info_500kb.csv", quote = FALSE, row.names = FALSE)
## -------------------------------------------------------------------------------------------
