#!/usr/bin/env Rscript

#SBATCH --account=grdi_genarcc

args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 4)
  stop("Usage: rda.R <ALLELEFREQS> <PCNMDF> <BIOCLIM> <BIOVARS>")

# Install and load necessary packages ------------------------------------------

#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

for (p in c("gradientForest","tidyverse","data.table")) {
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
#freqs <- freqs
head(freqs[,c(1:10)])
rownames(freqs)
dim(freqs)
sum(is.na(freqs))

message("Read allele frequencies")

# Bioclim ----------------------------------------------------------------------

bioclim <- read.csv(file = args[3]) %>%
  filter(period == "1970-2000") %>%
  column_to_rownames("Site") %>%
  dplyr::select(-c("period", "ssp"))

bio_vars <- read.delim(args[4], header = T, sep = "\t")

bioclim <- bioclim[match(rownames(freqs), rownames(bioclim)),]
message("read bioclim data")
sum(is.na(bioclim))
rownames(bioclim)

all.equal(rownames(bioclim), rownames(freqs))

# Run GF ----------------------------------------------------------------------

set.seed(95)

maxLevel <- floor(log2(0.368*nrow(bioclim)/2))

message("Running gradient forest model")

gradfor <- gradientForest(cbind(bioclim, freqs),
                          predictor.vars = colnames(bioclim),
                          response.vars  = colnames(freqs),
                          ntree = 250, trace = T,
                          maxLevel = maxLevel,
                          corr.threshold = 0.50)

message("Gradient forest model complete")

# Save for plotting, diagnostics and predictive purposes (run locally).
saveRDS(object = gradfor, "14_rda/gforest_4bio2PCA_4bio3rda3SD_afGT_n819.RDS")

message("Done!")