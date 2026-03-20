#!/usr/bin/env Rscript

args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 6)
  stop("Usage: rda.R <ALLELEFREQS> <PCADF> <PCAs> <BIOCLIM> <BIOVARS> <OUTPUT>")

# Install and load necessary packages ------------------------------------------

library(data.table); library(tidyverse); library(vegan)

#for (p in c("vegan","tidyverse","data.table")) {
#  if (!suppressMessages(require(p, character.only = T))) {
#    message(paste("Installing:", p))
#    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
#    suppressMessages(require(p, character.only = T))}
#  rm(p)
#}

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

# PCAS ------------------------------------------------------------------------

pcas <- read.csv(file = args[2], row.names = "Pop")

pcas_to_use <- c(read.delim(args[3])[,1])

pcas <- pcas[,c("PC1", "PC2")]
head(pcas)
rownames(pcas)

message("Read pca dataframe")

pcas <- pcas[match(rownames(freqs), rownames(pcas)),]

dim(pcas); head(pcas)

# Bioclim ----------------------------------------------------------------------

biovars <- c(read.delim(args[5])[,1])

bioclim <- read.csv(file = args[4]) %>%
  filter(period == "1970-2000") %>%
  column_to_rownames("Site") %>%
  select(c(biovars))

bioclim <- bioclim[match(rownames(freqs), rownames(bioclim)),]
message("Read bioclim data")

rownames(freqs); rownames(pcas); rownames(bioclim)

all.equal(rownames(bioclim), rownames(pcas))
all.equal(rownames(bioclim), rownames(freqs))

# Run RDA ----------------------------------------------------------------------

(rda_formula <- as.formula(paste0("freqs ~ ", paste0(colnames(bioclim), collapse = " + "), " + Condition(", paste0(colnames(pcas), collapse = " + "), ")")))

(biorda <- vegan::rda(formula = rda_formula, data = cbind(bioclim, pcas), scale = TRUE))

RsquareAdj(biorda)

outdir <- args[6]
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(args[6], "/")
}

saveRDS(object = biorda, paste0(outdir, "bio_rda_popAfsGL_4bio3RDA3SD_2PCA.RDS"))


# RDA values/outputs -----------------------------------------------------------

# Variance explained by each RDA axis.
(eigvsumm <- round(summary(eigenvals(biorda, model = "constrained")), 5))
write.table(eigvsumm, paste0(outdir, "rda_eigv.txt"), quote = F, row.names = F)

# Screeplot for variance explained by each axis.
png(paste0(outdir, "rda_screeplot.png"), width = 1200, height = 1000, res = 100)
screeplot(biorda); dev.off() # Variation explained by each RDA axis.

# Model VIF.
message("Writing VIF information")
(rda_vif <- round(vif.cca(biorda), 2)) # Check model variance inflation factors.
write.csv(rda_vif, paste0(outdir, "rda_vif.csv"), quote = F, row.names = T)

# RDA SNP loadings
message("Writing loading values")
loadings <- scores(biorda, choices = c(1:5), display = 'species')
png(paste0(outdir, "rda_loadings.png"), width = 1200, height = 1200, res = 100)
par(mfrow = c(2, 2))
hist(loadings[,1], main = "RDA1 Loadings", xlab = "RDA1")
hist(loadings[,2], main = "RDA2 Loadings", xlab = "RDA2")
hist(loadings[,3], main = "RDA3 Loadings", xlab = "RDA3"); dev.off()
snp_loadings <- as.data.frame(loadings)
write.csv(snp_loadings, paste0(outdir, "snp_loadings_full.csv"), row.names = TRUE, quote = F)

# RDA site scores (populations)
site_scores <- as.data.frame(scores(biorda, display = 'sites', choices = c(1:5), scaling = 3))
write.csv(site_scores, paste0(outdir, "rda_site_scores.csv"), row.names = T, quote = F)

# RDA environmental variable loadings
bio_scores <- as.data.frame(scores(biorda, display = 'bp', choices = c(1:5), scaling = 3))
write.csv(bio_scores, paste0(outdir, "rda_bio_scores.csv"), row.names = T, quote = F)

# Variance partitioning.
(vp <- varpart(freqs, bioclim, pcas))
saveRDS(object = vp,  paste0(outdir, "rda_varpart.RDS"))