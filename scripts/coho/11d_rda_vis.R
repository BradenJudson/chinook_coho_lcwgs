
library(tidyverse); library(viridis); library(ggrepel)

# ------------------------------------------------------------------------------

# Read in SNP loadings on each axis of the RDA.
snp_loadings <- read.csv("data/rdas/afsGT_4bioclim_3pcnms/snp_loadings_full.csv") 

# Visualize distribution of SNP loadings.
ggplot(data = snp_loadings %>% 
         pivot_longer(starts_with("RDA")), 
       aes(x = value)) + 
  geom_histogram(bins = 150) +
  facet_wrap(~ name) + 
  theme_bw() +
  labs(x = "Loading", y = "Count")

ggsave("plots/rda_loadings_full.tiff", dpi = 300, width = 12, height = 8)


# Biplot -----------------------------------------------------------------------

# Assemble RDA biplot information for ggplot.
# Assemble RDA biplot information for ggplot.
site_info <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "coho") %>% 
  group_by(region_revised) %>% 
  mutate(mLat = mean(latitude))

coldf <- site_info[,c("region_revised", "plot_col", "plot_shape", "mLat")] %>% 
  group_by(region_revised) %>% sample_n(1) %>% ungroup() %>% arrange(mLat)


site_scores <- read.csv("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/rda_site_scores.csv", row.names = 1)
site_scores <- merge(site_scores, site_info, by.x = 0, by.y = "site")
# site_scores$region <- reorder(site_scores$region, site_scores$latitude)
bio_scores <- read.csv("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/rda_bio_scores.csv")
eigv <- read.delim("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/rda_eigv.txt", sep = "")

# Function for visualizing RDA w/ ggplot + tidyverse functions. 
rda_biplot <- \(x, y) {
  
  scalar <- 15
  
  varx <- paste0(text = enquo(x), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(1) %>% .[2,])), "%)")[2]
  vary <- paste0(text = enquo(y), " (" , sprintf("%0.1f", 100*c(eigv %>% dplyr::select(2) %>% .[2,])), "%)")[2]
  
  bio_scores <- bio_scores %>% mutate(
    label_x = RDA1 * (scalar + 1),
    label_y = RDA2 * (scalar + 1))
  
  ggplot(data = site_scores,
         aes(x = {{x}},
             y = {{y}},
             fill  = plot_col,
             text  = Row.names,
             shape = plot_shape)) +
    theme_bw() +
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
    ) +
    labs(x = varx, y = vary) +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    geom_segment(data = bio_scores,
                 aes(x = 0, y = 0,
                     colour = X,
                     xend = {{x}}*scalar,
                     yend = {{y}}*scalar),
                 lineend = "round", 
                 linewidth = 1,
                 arrow = arrow(length = unit(0.1, "inches")),
                 inherit.aes = FALSE,
                 show.legend = FALSE) +
    scale_colour_manual(values = c('bio5' = 'firebrick3',
                                   'bio12' = "blue1",
                                   'bio1' = 'goldenrod1',
                                   'bio15' = "forestgreen")) +
    geom_point(size = 2) +
    geom_text(data = bio_scores,
              aes(x = label_x + c(0, -0.4, 0, -0.3),
                  y = label_y + c(-0.1, -0.2, 0, -0.1),
                  label = X), 
              hjust = -0.1, 
              vjust = -0.1,
              inherit.aes = FALSE) +
    guides(fill  = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)),
           shape = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)))
  
}

# Plot biplot.
(rda12 <- rda_biplot(x = RDA1, y = RDA2) +
    scale_x_continuous(expand = c(0, 1.5)) +
    scale_y_continuous(expand = c(0.01, 0.1)))

ggsave("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/coho_afsGT_bio43RDA3SDbiplot.tiff", 
       dpi = 300, width = 8, height = 6)

saveRDS(rda12, "data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/rda12_biplot.RDS")


# ------------------------------------------------------------------------------
