library(tidyverse)

lcpi <- function(dir) paste0(dir,
                             list.files(pattern = ".pestPG",
                                        path = dir)) %>% 
  set_names(.) %>% 
  map_df(read_table, .id = "FileName") %>% 
  group_by(FileName) %>% 
  mutate(pi = tP/nSites) %>% 
  summarise(meanPi = mean(pi, na.rm = TRUE),
            sdPi   = sd(pi, na.rm = TRUE)) %>% 
  mutate(pop = as.factor(tools::toTitleCase(gsub("\\_"," ",
                         gsub(".thetas.*","", gsub(".*/", "", FileName))))))

ch_pi <- lcpi(dir = "data/pi/chinook_likelihoods/")

write.csv(ch_pi, "data/pi/chin_likelihoods_pi_n106.csv", row.names = FALSE)
write.csv(ch_pi[,c(4,2)], "chin_coho_shiny/data/chin_likelihoods_pi_n106.csv", row.names = F)
