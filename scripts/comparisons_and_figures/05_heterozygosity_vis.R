library(tidyverse)

ch <- read.csv("data/heterozygosity_files/chinook_imputed_het_n819.csv") %>% 
  mutate(sp = "Chinook", data = "Imputed")
co <- read.csv("data/heterozygosity_files/coho_imputed_het_n650.csv") %>% 
  mutate(sp = "Coho", data = "Imputed")

write.csv(ch %>% group_by(population) %>% 
  summarise(mHet = mean(pOHET)), "data/heterozygosity_files/chinook_imputed_popHet_n106.csv")

dat <- rbind(ch, co)

sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  group_by(region_revised) %>% mutate(mLat = mean(latitude)) %>% 
  rename("population" = "site") 


full_imputed <- left_join(dat, sites) %>% 
  mutate(region = fct_reorder(region_revised, pOHET, .desc = FALSE))


het_plot_reg <- \(df) {
  
  (het_plot <- ggplot(data = df,
                      aes(x = region, y = pOHET,
                          fill = mLat)) +
     geom_point(alpha = 1/2, shape = 21) +
     scale_fill_viridis(name = "Latitude") +
     geom_boxplot(alpha = 4/5, outliers = FALSE) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
           plot.margin = unit(c(1,1,1,2), "lines"),
           legend.position = "inside",
           legend.position.inside = c(0.045, 1/3),
           legend.background = element_rect(colour = 'black', linewidth = 1/5)) +
     labs(x = NULL, y = "Observed Heterozygosity"))
  
}

(chinook_imp <- het_plot_reg(full_imputed) + 
    facet_wrap(~sp, ncol = 1))
ggsave("plots/regional_hets_imputed.tiff", dpi = 300,
       width = 10, height = 8)


# -------------------------------------------------------------------------

gen_div <- do.call("rbind", list(
  read.csv("data/chinook_imputed_het_n106.csv") %>% 
    mutate(species = "Chinook", genotypes = "Imputed genotypes"),
  read.csv("data/coho_imputed_het_n83.csv") %>% 
    mutate(species = "Coho", genotypes = "Imputed genotypes"),
  read.csv("data/chinook_lcwgs_hets_n106.csv", row.names = 1) %>% 
    mutate(species = "Chinook", genotypes = "Genotype likelihoods"),
  read.csv("data/coho_lcwgs_hets_n83.csv", row.names = 1) %>% 
    mutate(species = "Coho", genotypes = "Genotype likelihoods")
)) %>% pivot_wider(names_from = "genotypes", values_from = "het")

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
