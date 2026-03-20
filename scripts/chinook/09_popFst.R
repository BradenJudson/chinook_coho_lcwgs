
library(tidyverse); library(pheatmap)

# -------------------------------------------------------------------------

# Read in sample information.
popmap <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "chinook") %>% 
  filter(!fate %in% c("Discard (coverage too low)", "PCA outlier", "Drop population"))

# Write pop map in plink format for Fst calculations.
write.table(popmap[,c("final_bam", "final_bam", "population")],
            "data/pca_data/chinook_plink_popMap_n819.txt", sep = "\t", 
            quote = F, row.names = F, col.names = F)

# Calculate Hudson's Fst between all populations.
system("plink2 --vcf data/vcfs/chinook_lcwgs_maf005_n819_imputed.vcf.gz --double-id --aec --within data/pca_data/chinook_plink_popMap_n819.txt --fst CATPHENO method=hudson --out data/fst_data/chinook_imputed_n819fst")

# Read in Fst estimates and convert from long form to pairwise matrix.
fst <- read.delim("data/fst_data/chinook_imputed_n819fst.fst.summary")
fstmat <- ConGenFunctions::df2pmat(fst, var = "HUDSON_FST",
                                   grp1 = "X.POP1", grp2 = "POP2")
colnames(fstmat) <- rownames(fstmat) <- gsub("_.*", "", rownames(fstmat))

# Create and save pairwise Fst heatmap.
png(width = 1200, height = 1200, units = "px", filename = "plots/chinook_FstMat_imputed819.png")
(fst_heatmap <- pheatmap(mat = fstmat %>% 
                           `rownames<-`(trimws(str_replace_all(rownames(.), "([A-Z])", " \\1"))) %>% 
                           `colnames<-`(trimws(str_replace_all(colnames(.), "([A-Z])", " \\1"))),
         treeheight_row = 0, treeheight_col = 0,
         na_col = NA, legend = T,
         angle_col = 90)[[4]])
dev.off()
