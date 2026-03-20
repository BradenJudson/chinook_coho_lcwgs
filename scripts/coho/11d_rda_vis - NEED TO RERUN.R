
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
site_info   <- read.csv("data/ch_site_info.csv") %>%
  mutate(Site = gsub("ChilcotinLower", "Chilcotin", Site))
site_info$Lineage <- reorder(site_info$Lineage, site_info$Latitude)
site_scores <- read.csv("data/rdas/afsGT_4bioclim_3pcnms/rda_site_scores.csv", row.names = 1)
site_scores <- merge(site_scores, site_info, by.x = 0, by.y = "Site")
site_scores$Lineage <- reorder(site_scores$Lineage, site_scores$Latitude)
bio_scores  <- read.csv("data/rdas/afsGT_4bioclim_3pcnms/rda_bio_scores.csv")
eigv <- read.delim("data/rdas/afsGT_4bioclim_3pcnms/rda_eigv_unconstrained.txt", sep = "")

# Function for visualizing RDA w/ ggplot + tidyverse functions. 
rda_biplot <- \(x, y) {
  
  scalar <- 40
  
  varx <- paste0(text = enquo(x), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(1) %>% .[2,])), "%)")[2]
  vary <- paste0(text = enquo(y), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(2) %>% .[2,])), "%)")[2]
  
  ggplot() +
    geom_point(data = site_scores,
               aes(x = {{x}},
                   y = {{y}},
                   fill = factor(Lineage),
                   text = Row.names,
                   shape = factor(Lineage)),
               size = 3) +
    scale_shape_manual(values = c(rep(c(21,23), 10))) +
    scale_fill_viridis(discrete = T) +
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
          # legend.position.inside = c(0.05, 0.11),
          legend.background = element_blank()) +
    geom_segment(data = bio_scores,
                 aes(x = 0, y = 0,
                     xend = {{x}}*scalar,
                     yend = {{y}}*scalar),
                 colour = '#0868ac',
                 lineend = "round", 
                 linewidth = 1,
                 arrow = arrow(length = unit(0.1, "inches"))) +
    # Manually label a few points of interest that fall farthest away from the centroid of the distribution.
    geom_text_repel(data = site_scores[site_scores$Pop %in% c("Imnaha", "Granite", "Okanagan", "Trinity"),],
                     direction = "x",point.padding = 5,
                     aes(x = {{x}},
                         y = {{y}},
                         label = Pop))
  
}

# Plot biplot.
(rda12 <- rda_biplot(x = RDA1, y = RDA2))

ggsave("plots/afsGT_indvbio43RDA3SDbiplot.tiff", 
       dpi = 300, width = 8, height = 6)

# ggsave("plots/afsGT_indvbio43RDA3SDbiplot2.tiff", 
#        dpi = 300, width = 8, height = 5.2)

# See which pop is which.
ggplotly(rda12, tooltip = "text")

# Visualize all axes of variation.
cowplot::plot_grid(plotlist = list(
  rda_biplot(x = RDA1, y = RDA2),
  rda_biplot(x = RDA2, y = RDA3),
  rda_biplot(x = RDA1, y = RDA3)
))

ggsave("data/rdas/afsGT_4bioclim_3pcnms/rda_biplots_multi.tiff",
       dpi = 300, width = 10, height = 10)

# ------------------------------------------------------------------------------
