#!/usr/bin/env Rscript

args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 3)
  stop("Usage: rda.R <ALLELEFREQS> <RDA> <AXES>")

# Install and load necessary packages ------------------------------------------

#install.packages("devtools",  repos = "https://mirror.its.dal.ca/cran", dependencies = T)
#library("devtools")
#install_github("jdstorey/qvalue")

for (p in c("vegan","tidyverse","data.table", "robust", "qvalue")) {
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

K <- as.numeric(args[3])

rdadapt <- function(rda, K) {
  
  zscores  <- rda$CCA$v[ ,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha  <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda   <- median(resmaha)/qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha/lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  
  return(data.frame(p.values = reschi2test,
                    q.values = q.values_rdadapt))
}


out <- rdadapt(biorda, K)
rownames(out) <- names(biorda$colsum)
write.csv(out, paste0("rda_mah_rdadapt_k", K, ".csv"), row.names = TRUE, quote = FALSE)

top001 <- out %>%
  mutate(out = as.factor(case_when(
    p.values < quantile(p.values, 0.001) ~ "Y",
    TRUE ~ "N"
  ))
  )

# Outlier AF matrix ---------------------------------------------------------------------------

message("Creating outlier allele frequency matrix")

nrow( out[out$p.values < 0.01, ])
summary(out$q.values)

outlier_afs <- freqs[,c(colnames(freqs) %in% rownames(top001[top001$out == "Y", ]))]

write.table(outlier_afs, 
            paste0("outlier_afs_matrix_3RDA4Bio_mahalanobis_k", K, "_top001.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

message("Done!")




