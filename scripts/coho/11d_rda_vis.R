
library(tidyverse); library(viridis); library(ggrepel)

# ------------------------------------------------------------------------------

# Read in SNP loadings on each axis of the RDA.
snp_loadings <- read.csv("data/rdas/afsGT_4bioclim_3pcnms/snp_loadings_full.csv") 

# Visualize distribution of SNP loadings.
ggplot(data = snp_loadings %>% 
         pivot_longer(starts_with("RDA")), 
       aes(x = value)) + 
  geom_histogram(bins = 150) +
  facet_wrap(~name) + 
  theme_bw() +
  labs(x = "Loading", y = "Count")

ggsave("plots/rda_loadings_full.tiff", dpi = 300, width = 12, height = 8)


# Biplot -----------------------------------------------------------------------

# Assemble RDA biplot information for ggplot.
site_info   <- readxl::read_excel("data/sample_sites.xlsx") %>%
  filter(species == "coho") %>% rename("region" = "region_revised")
site_info$region <- reorder(site_info$region, site_info$latitude)
site_scores <- read.csv("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda_site_scores.csv", row.names = 1)
site_scores <- merge(site_scores, site_info, by.x = 0, by.y = "site")
site_scores$region <- reorder(site_scores$region, site_scores$latitude)
bio_scores <- read.csv("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda_bio_scores.csv")
eigv <- read.delim("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda_eigv.txt", sep = "")

# Function for visualizing RDA w/ ggplot + tidyverse functions. 
rda_biplot <- \(x, y) {
  
  scalar <- 40
  
  varx <- paste0(text = enquo(x), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(1) %>% .[2,])), "%)")[2]
  vary <- paste0(text = enquo(y), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(2) %>% .[2,])), "%)")[2]
  
  ggplot() +
    scale_shape_manual(values = c(rep(c(23,21), 10))) +
    # scale_fill_viridis(discrete = T) +
    geom_text_repel(data = bio_scores,
                    aes(x = {{x}}*scalar,
                        y = {{y}}*scalar,
                        label = X), 
                    seed = 240, nudge_x = 1,
                    segment.color = 'transparent') +
    theme_bw() +
    labs(x = varx, y = vary) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.background = element_blank()) +
    geom_segment(data = bio_scores,
                 aes(x = 0, y = 0,
                     xend = {{x}}*scalar,
                     yend = {{y}}*scalar),
                 colour = '#0868ac',
                 lineend = "round", 
                 linewidth = 1,
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_point(data = site_scores,
               aes(x = {{x}},
                   y = {{y}},
                   fill = factor(region),
                   text = Row.names,
                   shape = factor(region)),
               size = 3) 
  
}

# Plot biplot.
(rda12 <- rda_biplot(x = RDA1, y = RDA2))

ggsave("plots/coho_afsGT_bio43RDA3SDbiplot.tiff", 
       dpi = 300, width = 8, height = 6)

saveRDS(rda12, "data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda12_biplot.RDS")


# ------------------------------------------------------------------------------
