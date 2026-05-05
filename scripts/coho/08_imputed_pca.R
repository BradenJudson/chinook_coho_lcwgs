library(tidyverse)

# Run PCA.
system("plink --vcf data/vcfs/coho_lcwgs_maf005_n650_imputed.vcf.gz --aec --double-id --pca 1000 --out data/pca_data/coho_imputed_pca_n650")

# Read in sample information.
indvs <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "coho") %>% 
  filter(!fate %in% c("Discard (coverage too low)", "Drop (PCA outlier)", "Drop population"))

#Read in site information and format correctly.
sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "coho") %>% 
  group_by(region_revised) %>% 
  mutate(mLat = mean(latitude))


unique(indvs[!indvs$population %in% sites$site, "population"])
unique(sites[!sites$site %in% indvs$population, "site"])

dat <- merge(x = indvs,
             y = sites[,c("site", "region_revised", "mLat", "plot_shape", "plot_col")],
             by.x = "population", by.y = "site")

imp_pca <- read.table("data/pca_data/coho_imputed_pca_n650.eigenvec") %>% 
  .[,2:ncol(.)] %>%  # Remove unnecessary column.
  `colnames<-`(., c("id", paste0("PC", 1:(ncol(.)-1)))) %>% 
  mutate(id = gsub(".*/", "", gsub("\\..*", "", id))) %>% 
  merge(., dat, by = "id") %>% 
  mutate(region_revised = as.factor(region_revised)) 

scree <- read_tsv("data/pca_data/coho_imputed_pca_n650.eigenval", 
                  col_names = "eigenval") %>% 
  mutate(axis = as.numeric(rownames(.)),
         PC = paste0("PC", axis)) %>% 
  arrange(axis) %>% 
  mutate(var_exp = round(eigenval/sum(eigenval), 3)) %>% 
  mutate(PC = factor(PC, levels = unique(PC)))

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


saveRDS(object = pca_plot, "data/pca_data/coho_imputed_n650_PCAobj.RDS")
ggsave("plots/coho_pca12_imputed650.tiff", dpi = 300, width = 9, height = 6)
coho_imp_pca_data <- imp_pca[,c("id", "population", "region_revised", "mLat", paste0("PC", c(1:20)))] 
write.csv(coho_imp_pca_data, "data/pca_data/coho_pca_data_n650_imputed.csv", row.names = F)
write.csv(coho_imp_pca_data, "chin_coho_shiny/data/coho_imputed_pca_scores_n650.csv", row.names = F)

