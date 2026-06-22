library(gradientForest); library(tidyverse); library(sf)
library(geodata); library(ggnewscale); library(cowplot)

set.seed(240)

# Putatively adaptive SNPs -----------------------------------------------------
# Allele frequencies -----------------------------------------------------------

pops <- read.table("data/coho_af_rownames_n83.txt")

outlier_afs <- read.table("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_19bioGF/outlier_afs_matrix_3RDA4Bio_mahalanobis_k3_top001.txt",
                          header = TRUE) %>% `rownames<-`(., pops[,1])

# GF ---------------------------------------------------------------------------

# Read gradient forest model.
gradfor <- readRDS("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/mah_topp001_outliers_GF_4bioclim.RDS")

# Assess predictor importance for the GF model.
png("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/importance_scores.png", 
    width = 2000, height = 1200, res = 250)
plot(gradfor)
dev.off()

sum(colnames(gradfor$Y) == colnames(outlier_afs)) == ncol(outlier_afs)
sum(rownames(gradfor$X) == rownames(outlier_afs)) == length(unique(pops[,1]))

fclim <- \(scenario) {
  
  bioclim <- read.csv("data/coho_bioclim.csv") %>% 
    filter(if(is.na(scenario)) is.na(ssp) else ssp == scenario) %>% 
    column_to_rownames("Site") %>%
    dplyr::select(colnames(gradfor$X))
  
  bioclim <- bioclim[match(rownames(outlier_afs), rownames(bioclim)),]
  
}

fpr <- fclim(NA)
f26 <- fclim(26)
f45 <- fclim(45)
f85 <- fclim(85)

pCurrent <- predict(gradfor, fpr)
pred26 <- predict(gradfor, f26)
pred45 <- predict(gradfor, f45)
pred85 <- predict(gradfor, f85)

go26 <- data.frame(go26 = sqrt(rowSums((pred26 - pCurrent)^2)))
go45 <- data.frame(go45 = sqrt(rowSums((pred45 - pCurrent)^2)))
go85 <- data.frame(go85 = sqrt(rowSums((pred85 - pCurrent)^2)))

genomic_offsets <- merge(go26, go45, by = 0) %>% 
  merge(., go85, by.x = "Row.names", by.y = 0) %>% 
  dplyr::rename("site" = Row.names)

write.csv(genomic_offsets, "data/GOs/coho_offset_19bioMahtop001_afsGTs_2pca_n650.csv", row.names = F)
write.csv(genomic_offsets, "chin_coho_shiny/coho_offset_19bio3SD_afsGTs_2pca_n650.csv", row.names = F)
