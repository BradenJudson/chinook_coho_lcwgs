library(tidyverse); library(ggpmisc)


# Assemble data -----------------------------------------------------------

chin_off <- read.csv("data/GOs/chinook_offset_19bio3SD_afsGT_2pca_n819.csv")[,c(1,4)]
chin_het <- read.csv("data/heterozygosity_files/chinook_imputed_popHet_n106.csv", row.names = 1)
chin_bio <- read.csv("data/chinook_bioclim.csv", row.names = 1) %>% 
  filter(is.na(ssp)) %>% dplyr::select(c(Site, bio1, bio5, bio12, bio15))
ch_sites <- read.csv("data/chin_site_info.csv")[,c(1,4,7)]
ch_sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "chinook") %>% dplyr::select(c(2,4,7))

chinook_data <- merge(merge(chin_off, chin_bio, by = "Site"), 
                      ch_sites, by.x = "Site", by.y = "site") %>% 
  merge(., chin_het, by.x = "Site", by.y = "population") %>% 
  mutate(Lineage = fct_reorder(region_revised, latitude, .fun = mean)) %>% 
  dplyr::select(-region_revised) %>% 
  pivot_longer(cols = -c(Site, go85, Lineage), names_to = "variable") 

coho_off <- read.csv("data/GOs/coho_offset_19bio3SD_afsGTs_2pca_n650.csv")[,c(1,4)]
coho_het <- read.csv("data/heterozygosity_files/coho_imputed_popHet_n83.csv", row.names = 1)
coho_bio <- read.csv("data/coho_bioclim.csv", row.names = 1) %>% 
  filter(is.na(ssp)) %>% dplyr::select(c(Site, bio1, bio5, bio12, bio15))
co_sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "coho") %>% dplyr::select(c(2,4,7))

coho_data <- merge(merge(coho_off, coho_bio, by.x = "site", by.y = "Site"),
                   co_sites, by = "site") %>% 
  merge(., coho_het, by.x = "site", by.y = "population") %>% 
  mutate(Lineage = fct_reorder(region_revised, latitude, .fun = mean)) %>% 
  dplyr::select(-region_revised) %>% 
  pivot_longer(cols = -c(site, go85, Lineage), names_to = "variable") 


# Plotting function and labels --------------------------------------------

labs <- c(
  "bio1" = "Annual mean temperature (°C) (bio1)",
  "bio5" = "Max temperature of the warmest month (°C) (bio5)",
  "bio12" = "Annual precipitation (mm) (bio12)",
  "bio15" = "Precipitation seasonality (%) (bio15)",
  "latitude" = "Latitude",
  "mHet" = "Population average heterozygosity"
)


plot_corrs <- \(df) {
  
  ggplot(data = df,
         aes(y = go85, x = value)) +
    geom_smooth(method = "lm", linetype = 2,
                colour = "gray25", linewidth = 1/2) +
    geom_point(aes(fill = Lineage),
               shape = 21, alpha = 4/5,
               size = 3/2) +
    facet_wrap(~variable, 
               scales = "free_x",
               labeller = as_labeller(labs)) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    stat_poly_eq(use_label(c("R2", "P"))) +
    labs(x = NULL, y = "Genomic offset (SSP8.5)") +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE, 
                               override.aes = list(size = 2))) +
    theme(legend.text = element_text(size = 10))
  
}

(chinook_corrs <- plot_corrs(chinook_data))
ggsave("plots/chin_go85_corrs.tiff", dpi = 300,
       width = 12, height = 6)

(coho_corrs <- plot_corrs(coho_data))
ggsave("plots/coho_go85_corrs.tiff", dpi = 300,
       width = 12, height = 6)



