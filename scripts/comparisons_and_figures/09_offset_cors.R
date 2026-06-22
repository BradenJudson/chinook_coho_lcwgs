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


# -------------------------------------------------------------------------
full_sites <- readxl::read_excel("data/sample_sites.xlsx")[,c(7:9)] %>% 
  group_by(region_revised) %>% summarize(plot_col = unique(plot_col),
                                         plot_shape = unique(plot_shape))

chin_off_reg <- chinook_data %>% 
  mutate(rank_off = min_rank(go85))
  # summarise(m85_ch = mean(go85),
  #           s85_ch = sd(go85)/sqrt(length(region_revised)))
coho_off_reg <- coho_data %>%
  mutate(rank_off = min_rank(go85))
  # group_by(region_revised) %>% 
  # summarise(m85_co = mean(go85),
  #           s85_co = sd(go85)/sqrt(length(region_revised)))

off_reg <- merge(chin_off_reg, coho_off_reg)
off_can <- off_reg[!off_reg$region_revised %in% c("Alaska", "Coastal Washington", "Southeast Alaska", 
                                                 "Columbia River (Lower)", "Oregon", "Puget Sound", "Haida Gwaii"),] %>% 
  merge(., full_sites)




off_can <- off_reg[!off_reg$region_revised %in% c("Columbia River (Lower)", "Puget Sound", "Haida Gwaii", 
                                                  "Coastal Washington", "Oregon"),] %>% 
  merge(., full_sites)



ggplot(data = off_can,
       aes(x = m85_ch,
           y = m85_co)) +
  geom_segment(aes(x = m85_ch - s85_ch,
                   xend = m85_ch + s85_ch,
                   y = m85_co, yend = m85_co),
               inherit.aes = F, show.legend = F) +
  geom_segment(aes(x = m85_ch, xend = m85_ch,
                   y = m85_co - s85_co,
                   yend = m85_co + s85_co),
               inherit.aes = F, show.legend = F) +
  geom_smooth(method = "lm", colour = "skyblue",
              linetype = 2, alpha = 1/6) +
  geom_point(size = 2, 
             aes(fill = plot_col,
                 shape = plot_shape)) +
  labs(x = "Chinook regional offsets",
       y = "Coho regional offsets") +
  theme_bw() +
  stat_poly_eq(use_label(c("R2", "P"))) +
  scale_fill_identity(
    guide = "legend",
    breaks = off_can$plot_col,
    labels = off_can$region_revised,
    name = NULL
  ) +
  scale_shape_identity(
    guide = "legend",
    breaks = off_can$plot_shape,
    labels = off_can$region_revised,
    name = NULL
  ) +
  guides(fill  = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)),
         shape = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)))


ggsave("plots/genomic_offset_regional_corrs.tiff",
       width = 10, height = 6, bg = 'white', dpi = 300)

j <- rbind(coho_off_reg[,c(1,2,8,11)] %>% mutate(sp = "Coho"),
           chin_off_reg[,c(1,2,8,11)] %>% mutate(sp = "Chinook") %>% 
             dplyr::rename("site" = "Site")) 

rr <- c("Columbia River (Interior ocean-type)",
        "Columbia River (Lower)",
        "California", "Columbia River (Interior stream-type)",
        "Haida Gwaii", "Oregon", "Southeast Alaska")

ggplot(data = j[!j$region_revised %in% rr,],
               aes(x = region_revised,
                   y = rank_off,
                   fill = sp)) +
  geom_boxplot(whisker.linetype =2,
               outlier.fill = NULL,
               outlier.shape = 21,
               outlier.size = 1.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/2),
        legend.title = element_blank()) +
  labs(x = NULL) +
  geom_pwc(
    method = "wilcox.test",
    label = "p.format",
    paired = F
  ) 




chin_off_reg <- chinook_data %>% 
  group_by(region_revised) %>% 
  filter(!region_revised %in% c("Columbia River (Lower)", "Puget Sound", 
                                "Coastal Washington", "Oregon", "Alaska",
                                "Southeast Alaska")) %>% 
  summarise(m85_ch = mean(go85),
             s85_ch = sd(go85)/sqrt(length(region_revised))) %>% 
  mutate(rank_off_ch = scales::rescale(m85_ch))
coho_off_reg <- coho_data %>%
  group_by(region_revised) %>% 
  filter(!region_revised %in% c("Columbia River (Lower)", "Puget Sound", 
                                "Coastal Washington", "Oregon", "Alaska",
                                "Southeast Alaska")) %>% 
  summarise(m85_co = mean(go85),
            s85_co = sd(go85)/sqrt(length(region_revised))) %>% 
  mutate(rank_off_co = scales::rescale(m85_co))



off_reg <- merge(chin_off_reg, coho_off_reg) %>% 
  merge(., full_sites) %>% 
  filter(!region_revised %in% c("Columbia River (Lower)", "Puget Sound", 
                               "Coastal Washington", "Oregon", "Alaska"))

ggplot(data = off_reg,
       aes(x = rank_off_ch,
           y = rank_off_co)) +
  # geom_segment(aes(x = m85_ch - s85_ch,
  #                  xend = m85_ch + s85_ch,
  #                  y = m85_co, yend = m85_co),
  #              inherit.aes = F, show.legend = F) +
  # geom_segment(aes(x = m85_ch, xend = m85_ch,
  #                  y = m85_co - s85_co,
  #                  yend = m85_co + s85_co),
  #              inherit.aes = F, show.legend = F) +
  geom_smooth(method = "lm", colour = "skyblue",
              linetype = 2, alpha = 1/6) +
  geom_point(size = 3, 
             aes(fill = plot_col,
                 shape = plot_shape)) +
  labs(x = "Chinook regional offsets (scaled)",
       y = "Coho regional offsets (scaled)") +
  theme_bw() +
  stat_poly_eq(use_label(c("R2", "P"))) +
  scale_fill_identity(
    guide = "legend",
    breaks = off_reg$plot_col,
    labels = off_reg$region_revised,
    name = NULL
  ) +
  scale_shape_identity(
    guide = "legend",
    breaks = off_reg$plot_shape,
    labels = off_reg$region_revised,
    name = NULL
  ) +
  guides(fill  = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)),
         shape = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2)))


