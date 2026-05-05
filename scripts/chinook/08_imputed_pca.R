
library(tidyverse); library(scales)

# -------------------------------------------------------------------------

# Run PCA.
system("plink --vcf data/vcfs/chinook_lcwgs_maf005_n819_imputed.vcf.gz --aec --double-id --pca 1000 --out data/pca_data/chinook_imputed_n819")

#Read in site information and format correctly.
sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "chinook") %>% 
  group_by(region_revised) %>% 
  mutate(mLat = mean(latitude))

indvs <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "chinook") %>% 
  filter(!fate %in% c("Discard (coverage too low)", "PCA outlier", "Drop population"))

dat <- merge(indvs, sites, by.x = "population", by.y = "site") 

scree <- read_tsv("data/pca_data/chinook_imputed_n819.eigenval", 
                  col_names = "eigenval") %>% 
  mutate(axis = as.numeric(rownames(.)),
         PC = paste0("PC", axis)) %>% 
  arrange(axis) %>% 
  mutate(var_exp = round(eigenval/sum(eigenval), 3)) %>% 
  mutate(PC = factor(PC, levels = unique(PC)))

imp_pca <- read.table("data/pca_data/chinook_imputed_n819.eigenvec") %>% # Read in PCA data.
  .[,2:ncol(.)] %>%                                                      # Remove unnecessary column.
  `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>%              # Rename columns.
  mutate(id = gsub("Kand0309-", "Kand0309\\.", gsub("\\..*$","",
              gsub(".*IDT_i5_[0-9]{1,3}.", "", gsub("(.*/)", "", id))))) %>%        # Make label formatting consistent.
  merge(., dat, by = "id") %>% 
  mutate(region_revised = factor(region_revised))

# Arrange factors with respect to Latitude (orange = south to pink = north)
imp_pca$region_revised <- reorder(imp_pca$region_revised, imp_pca$mLat)

coldf <- imp_pca[,c("region_revised", "plot_col", "plot_shape", "mLat")] %>% 
  group_by(region_revised) %>% sample_n(1) %>% ungroup()

(pca_plot <- ggplot(data = imp_pca,
                    aes(x = PC1, y = PC2,
                        shape = plot_shape,
                        fill = plot_col)) +
    geom_point(size = 1.5) +
    scale_shape_identity(
      guide = "legend",
      breaks = coldf$plot_shape,
      labels = coldf$region_revised,
      name = NULL
    ) +
    scale_fill_identity(
      guide = "legend",
      breaks = coldf$plot_col,
      labels = coldf$region_revised,
      name = NULL
    ) + theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 10)) +
    guides(fill  = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)),
           shape = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2))) +
    labs(x = paste0("PC1 (", sprintf("%.1f", scree$var_exp[1]*100), "%)"),
         y = paste0("PC2 (", sprintf("%.1f", scree$var_exp[2]*100), "%)")))


ggsave("plots/chinook_pca12_imputed819.tiff", dpi = 300, width = 9, height = 6)
saveRDS(object = pca_plot, "data/pca_data/chinook_imputed_n819_PCAobj.RDS")
pca_dat <- pca_plot@data[,c("id", "population", "region_revised", "mLat", paste0("PC", 1:20))] 
write.csv(pca_dat, "data/pca_data/chinook_imputed_pca_scores_n819.csv", row.names = F)
write.csv(pca_dat, "chin_coho_shiny/data/chinook_imputed_pca_scores_n819.csv", row.names = F)
