library(tidyverse)

# -------------------------------------------------------------------------

# Read in heterozygosity information.
ch <- read.csv("data/heterozygosity_files/chinook_imputed_het_n819.csv") %>% 
  mutate(sp = "Chinook", data = "Imputed")
co <- read.csv("data/heterozygosity_files/coho_imputed_het_n650.csv") %>% 
  mutate(sp = "Coho", data = "Imputed")

# Average population heterozygosity for Chinook.
write.csv(ch %>% group_by(population) %>% 
  summarise(mHet = mean(pOHET)), "data/heterozygosity_files/chinook_imputed_popHet_n106.csv")

dat <- rbind(ch, co)

sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  group_by(region_revised) %>% mutate(mLat = mean(latitude)) %>% 
  rename("population" = "site") 

full_imputed <- left_join(dat, sites) %>% 
  mutate(region = fct_reorder(region_revised, mLat, .desc = FALSE))

coldf <- data.frame(
  region = unique(as.character(sites$region_revised)),
  alpha = c(rep(c(0.3, 1), length(unique(sites$region_revised))/2), 0.3),
  colour = hue_pal()(length(unique(sites$region_revised)))
)


# Plot heterozygosity by region by species.
(het_plot <- ggplot(data = full_imputed,
                    aes(x = region, y = pOHET,
                        fill = plot_col)) +
    geom_point(alpha = 1/2, shape = 21) +
    scale_fill_identity(guide = "legend",
                        labels = coldf$region_revised) +
    geom_boxplot(alpha = 4/5, outliers = FALSE, show.legend = F) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.margin = unit(c(1,1,1,2), "lines"),
          legend.position = "none",
          legend.title = element_blank()) +
    labs(x = NULL, y = "Observed Heterozygosity") +
    guides(fill  = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2, alpha = 1)),
           shape = guide_legend(ncol = 1, reverse = TRUE, override.aes = list(size = 2, alpha = 1))))


(impHe <- het_plot + facet_wrap(~sp, ncol = 1))
ggsave("plots/regional_hets_imputed.tiff", dpi = 300,
       width = 10, height = 6, bg = 'white')

# -------------------------------------------------------------------------

# Compare imputed and lcWGS heterozygosity estimates by species. 
gen_div <- do.call("rbind", list(
  read.csv("data/heterozygosity_files/chinook_imputed_popHet_n106.csv") %>% 
    mutate(species = "Chinook", genotypes = "Imputed genotypes"),
  read.csv("data/heterozygosity_files/coho_imputed_popHet_n83.csv") %>% 
    mutate(species = "Coho", genotypes = "Imputed genotypes"),
  read.csv("data/heterozygosity_files/chinook_lcwgs_hets_n819.csv") %>% 
    mutate(species = "Chinook", genotypes = "Genotype likelihoods"),
  read.csv("data/heterozygosity_files/coho_lcwgs_hets_n650.csv") %>% 
    mutate(species = "Coho", genotypes = "Genotype likelihoods")
)) %>% pivot_wider(names_from = "genotypes", values_from = "mHet")

ggplot(data = gen_div, 
       aes(y = `Imputed genotypes`,
           x = `Genotype likelihoods`,
           group = species)) +
  geom_abline(slope = 1, linetype = 2) +
  geom_smooth(method = "lm", alpha = 1/6,
              aes(colour = species),
              show.legend = FALSE) +
  geom_point(aes(fill = species),
             shape = 21) +
  theme_bw() +
  labs(y = "Observed heterozygosity (imputed genotypes)",
       x = "Observed heterozygosity (genotype likelihoods)") +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9),
        legend.background = element_blank(),
        legend.key = element_blank())

ggsave("plots/heterozygosity_comp.tiff", dpi = 300,
       width = 8, height = 6)
