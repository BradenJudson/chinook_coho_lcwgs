
library(tidyverse); library(viridis); library(ggrepel)

# ------------------------------------------------------------------------------

# Read in SNP loadings on each axis of the RDA.
snp_loadings <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/snp_loadings_full.csv")

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
site_info   <- read.csv("data/chin_site_info.csv") %>%
  mutate(Site = gsub("ChilcotinLower", "Chilcotin", Site))
site_info$Lineage <- reorder(site_info$Lineage, site_info$Latitude)

site_scores <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda_site_scores.csv", row.names = 1)
site_scores <- merge(site_scores, site_info, by.x = 0, by.y = "Site")
site_scores$Lineage <- reorder(site_scores$Lineage, site_scores$Latitude)

bio_scores  <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda_bio_scores.csv")
eigv <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda_eigv.txt", sep = "")

# Function for visualizing RDA w/ ggplot + tidyverse functions. 
rda_biplot <- \(x, y) {
  
  scalar <- 33
  
  varx <- paste0(text = enquo(x), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(1) %>% .[2,])), "%)")[2]
  vary <- paste0(text = enquo(y), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(2) %>% .[2,])), "%)")[2]
  
  bio_scores <- bio_scores %>% mutate(
    label_x = RDA1 * (scalar + 1),
    label_y = RDA2 * (scalar + 1),
    angle = atan2(RDA2*scalar, RDA1*scalar)*180/pi) %>% 
    mutate(label_y = ifelse(angle > 90 | angle < -90, label_y + 1/2, label_y)) %>% 
    mutate(angle = ifelse(angle > 90 | angle < -90, angle + 180, angle)) 
     
  ggplot() +
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
                   fill  = factor(Lineage),
                   text  = Row.names,
                   shape = factor(Lineage)),
               size = 3) +
    scale_shape_manual(values = c(rep(c(21,23), 12))) +
    scale_fill_viridis(discrete = T) +
    geom_text(data = bio_scores,
              aes(x = label_x,
                  y = label_y,
                  label = X,
                  angle = angle), 
              hjust = -0.1, vjust = -0.1) +
    guides(fill  = guide_legend(ncol = 1, reverse = TRUE),
           shape = guide_legend(ncol = 1, reverse = TRUE))
  
}

# Plot biplot.
(rda12 <- rda_biplot(x = RDA1, y = RDA2) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))))

ggsave("plots/chinook_afsGT_indvbio43RDA3SDbiplot.tiff", 
       dpi = 300, width = 8, height = 6)

saveRDS(rda12, "data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda12_biplot.RDS")

