library(tidyverse)

chin_off <- read.csv("data/GOs/chinook_offset_19bio3SD_afsGT_2pca_n819.csv")
chin_het <- read.csv("data/heterozygosity_files/chinook_imputed_popHet_n106.csv", row.names = 1)
chin_ohet <- read.csv("data/heterozygosity_files/chinook_imputed_het_n819_outliers.csv") %>% 
  group_by(population) %>% summarize(oHet = mean(pOHET))

het <- merge(chin_het, chin_ohet) %>% 
  merge(., chin_off, by.x = "population", by.y = "Site") %>% 
  pivot_longer(cols = c("oHet", "mHet")) 


j <- het %>% 
  merge(., chin_ohet, by= "population")

ggplot(data = j, aes(x = value, y = go85)) +
  geom_point(shape = 21, aes(fill = region)) + facet_wrap(~name, scales = "free") +
  geom_smooth(method = "lm")

summary(lm(formula = go85 ~ value, data = het[het$name == "oHet",]))
