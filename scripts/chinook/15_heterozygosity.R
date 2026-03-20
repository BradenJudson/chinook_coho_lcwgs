
library(tidyverse); library(sf); library(ggpmisc); library(viridis)

# -------------------------------------------------------------------------

# Read in population map and lineage information.
popmap <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "chinook") %>% 
  filter(!fate %in% c("Discard (coverage too low)", "PCA outlier", "Drop population"))

#Read in site information and format correctly.
sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "chinook") %>% 
  group_by(region_revised) %>% 
  mutate(mLat = mean(latitude))

# Arrange lineages by average latitude.
lin <- sites %>% group_by(region_revised) %>%
  summarize(mLat = mean(latitude)) %>%
  arrange(mLat)

# sites <- merge(sites, lin) %>%
#   arrange(mLat, latitude) %>%
#   mutate(region = factor(region_revised, levels = lin$region_revised)) %>%
#   mutate(sitenum = as.numeric(sitenum))

# Calculate individual-level heterozygosity using VCFtools --het.
# Adjust IDs and add population/lineage information.
# Calculate the number and proportion of heterozygous genotypes.
hets <- read.delim("data/heterozygosity_files/chinook_imputed_het_n819.het", sep = "") %>% 
  mutate(id = gsub("Kand0309-", "Kand0309\\.", gsub("\\..*", "", 
              gsub(".*/", "", gsub(".*IDT_i5_[0-9]{1,3}.", "", INDV))))) %>% 
  merge(., popmap, by = "id") %>% 
  merge(., sites, by.x = "population", by.y = "site") %>% 
  mutate(OHET = N_SITES - O.HOM., pOHET = OHET/N_SITES) %>% 
  dplyr::rename("region" = "region_revised") %>% 
  mutate(region = fct_reorder(region, OHET, .desc = FALSE))

write.csv(hets[,c("id","population", "region", "pOHET")], 
          "data/heterozygosity_files/chinook_imputed_het_n819.csv", row.names = F)

write.csv(hets %>% group_by(population) %>% summarize(het = mean(pOHET)), 
          "chin_coho_shiny/data/chinook_imputed_het_n106.csv", row.names = F)

