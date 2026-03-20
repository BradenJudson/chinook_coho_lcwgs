#!/usr/bin/env Rscript

args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 3)
  stop("Usage: rda.R <ALLELEFREQS> <RDA> <AXES>")

# Install and load necessary packages ------------------------------------------

for (p in c("vegan","tidyverse","data.table")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
    suppressMessages(require(p, character.only = T))}
  rm(p)
}

# Allele frequencies -----------------------------------------------------------

# Issues with fread and rownames.
# Get rownames using awk '{print $1}' args[1] > 01_info_files/af_rownames.txt
freqs <- data.frame(fread(file = args[1]))
frqpops <- read.table("01_info_files/af_rownames.txt", col.names = "population")
rownames(freqs) <- frqpops[,1]
freqs <- freqs %>% dplyr::select(-c("population"))
dim(freqs)
rownames(freqs)

message("Read allele frequencies")

# RDA --------------------------------------------------------------------------

message("Reading RDA object")

biorda <- readRDS(args[2])

message("RDA object loaded")

# Outlier loci -----------------------------------------------------------------

message("Identifying outlier SNPs")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

axes <- as.numeric(args[3])

loadings <- scores(biorda, choices = c(1:axes), display = 'species')

message(paste0("Determining outliers across ", axes, " RDA axes"))

cand_list <- list()
for (i in 1:axes) {
  
  candSNP <- outliers(loadings[,i], z = 3)
  
  cand_list[[i]] <- cbind.data.frame(rep(paste0("RDA",i),
                                         times = length(candSNP)),
                                     names(candSNP), unname(candSNP))
  colnames(cand_list[[i]]) <- c("axis", "snp", "loading")
  
  if (i > 1) {
    cand_list[[i]] <- cand_list[[i]] %>%
      filter(!snp %in% cand_list[[i-1]]$snp)
  }
  
  print(paste0(length(candSNP), " outliers on RDA", i), quote = F)
  
}

all_cands <- do.call("rbind", cand_list) %>%
  separate(col = "snp", sep = 12,
           into = c("chromosome", "position"), remove = FALSE) %>%
  mutate(chromosome = as.factor(gsub("_", "", chromosome)))

head(all_cands)

nrow(all_cands)

write.csv(all_cands, "outlier_snps_rda.csv", row.names = FALSE)

# Outlier AF matrix ---------------------------------------------------------------------------

message("Creating outlier allele frequency matrix")

outlier_afs <- freqs[,c(colnames(freqs) %in% c(all_cands$snp))]

write.table(outlier_afs, "outlier_afs_matrix_afsGL4bio2PCA_n819.txt", row.names = FALSE, quote = FALSE, sep = "\t")

message("Done!")