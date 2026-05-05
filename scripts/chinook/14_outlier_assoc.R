
library(tidyverse)

# -------------------------------------------------------------------------

# Read in all SNP loadings. 
loadings <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/snp_loadings_full.csv")

# Read in outlier frequency matrix.
outliers <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/outlier_afs_matrix_afsGT4bio2PCA_n819.txt")

# Isolate SNP loadings for outliers only and convert to long form.
outlier_loadings <- loadings[loadings$X %in% colnames(outliers),] %>% 
  dplyr::select(!starts_with(c("PC", "RDA4"))) %>%
  pivot_longer(cols = starts_with("RDA"),
               names_to = "axis", 
               values_to = "loading") %>% 
  dplyr::rename("snp" = "X") 

# z <- outlier_loadings %>% 
#   mutate(chrom = gsub("\\.1_.*", ".1", X),
#          pos = as.numeric(gsub(".*\\.1_", "", X))) %>% 
#   rowwise() %>% 
#   mutate(maxRDA = max(c_across(starts_with("RDA"))))

# Clear up some memory.
rm(loadings); gc()

# Isolate contemporary climate data.
conBio <- read.csv("data/chinook_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% 
  column_to_rownames("Site") %>% 
  dplyr::select(c("bio1", "bio5", "bio12", "bio15"))


# Assemble an empty matrix for the following for-loop.
bioCors <- matrix(nrow = length(unique(outlier_loadings$snp)), ncol = ncol(conBio))
colnames(bioCors) <- colnames(conBio)
rownames(bioCors) <- unique(outlier_loadings$snp)



# Calculate correlation between SNP frequencies and each bioclim variable.
for (i in 1:length(unique(outlier_loadings$snp))) {
  
  snpfrq <- outliers[,c(unique(outlier_loadings$snp)[i])]
  bioCors[i,] <- apply(conBio, 2, function(x) cor(x, snpfrq))

}

# Set a delta threshold for determining whether outliers are strongly
# associated with one or >1 predictor variables. Group bioclim
# into two broad groups for similarly associated SNPs.
delta <- 0.15
prec <- c("bio12", "bio15")
temp <- c("bio5",  "bio1")

# For each outlier SNP loading, find the correlation with each bioclim variable.
# Find the highest and second highest correlation and group into one of four categories:
# i) the bioclim variable most strongly associated
# ii) `precipitation` if the SNP is similarly associated with bio12 and bio15
# iii) `temperature` if the SNP is similarly associated with bio1 and bio5
# iv) `unclassified` if the SNP is similarly associated with >1 SNP not from the groups in ii or iii.
cand_df_cor <- as.data.frame(bioCors) %>% 
  rownames_to_column(var = "SNP") %>% 
  mutate(CHROM = gsub(".1\\_[[:digit:]]*.", ".1", SNP),
         POS   = as.numeric(gsub(".*\\.1_", "", SNP))) %>% 
  rowwise() %>% 
  mutate(top_cor = max(c_across(starts_with("bio"))),
         sec_cor = sort(c_across(starts_with("bio")), decreasing = T)[2],
         top_var = names(.[c(colnames(bioCors))])[which.max(c_across(c(colnames(bioCors))))],
         sec_var = names(.[c(colnames(bioCors))])[which(rank(c_across(c(colnames(bioCors)))*-1) == 2)],
         association = case_when(
           abs(top_cor - sec_cor) >  delta ~ top_var,
           abs(top_cor - sec_cor) <= delta & top_var %in% temp & sec_var %in% temp ~ "Temperature",
           abs(top_cor - sec_cor) <= delta & top_var %in% prec & sec_var %in% prec ~ "Precipitation",
           abs(top_cor - sec_cor) <= delta & ((top_var %in% prec & sec_var %in% temp) | (top_var %in% temp & sec_var %in% prec)) ~ "Unclassified"
         )) %>% dplyr::select(c("CHROM", "POS", starts_with("bio"), "association"))

# Number of SNP-bioclim associations.
cand_df_cor %>% group_by(association) %>% tally() %>% mutate(prop = sprintf(n/sum(n)*100, fmt = "%#.2f"))

write.csv(cand_df_cor, "data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/outlier_snp_corrs.csv", row.names = F)
