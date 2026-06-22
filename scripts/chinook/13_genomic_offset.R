
library(gradientForest); library(tidyverse); library(data.table)
library(sf); library(geodata)

set.seed(240)

# Putatively adaptive SNPs -----------------------------------------------------
# Allele frequencies -----------------------------------------------------------

pops <- read.table("data/chinook_af_rownames_n106.txt", col.names = "pop")

outlier_afs <- read.table("data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/outlier_afs_matrix_3RDA4Bio_mahalanobis_k3_top001.txt", 
                          header = T) %>% `rownames<-`(., pops$pop)

# GF ---------------------------------------------------------------------------

# Read gradient forest model.
gradfor <- readRDS("data/rdas/chinook_afsGT_scaled4bio_3RDA_topMahp001_n819_19bioGF/gforest_4bio_topMAHp001_glAFs_19bioGF.RDS")

# Assess predictor importance for the GF model.
png("plots/chinook_gradfor_GTafs_mahp001_importance_n819.png", width = 2000, height = 1200, res = 250)
plot(gradfor)
dev.off()

# Check that gradient forest and outlier matrix match (same pops and same SNPs in the same order).
sum(colnames(gradfor$Y) == colnames(outlier_afs)) == ncol(outlier_afs)
sum(rownames(gradfor$X) == rownames(outlier_afs)) == length(unique(pops$pop))


fclim <- \(scenario) {
  
  bioclim <- read.csv("data/chinook_bioclim.csv", row.names = 1) %>% 
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
  merge(., go85, by.y = 0, by.x = "Row.names") %>% 
 dplyr::rename("Site" = Row.names)

write.csv(genomic_offsets, "data/GOs/chinook_offset_4bio3SD_afsGT_2pca_n819.csv", row.names = F)
write.csv(genomic_offsets, "chin_coho_shiny/data/chinook_offset_19bio3SD_afsGT_2pca_n819.csv", row.names = F)

rm(gradfor); gc()
