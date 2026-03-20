library(tidyverse)

coho <- read.csv("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/outlier_snp_bio_corrs.csv")
chin <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/outlier_snp_bio_corrs.csv")

group_out <- \(x) x %>% group_by(association) %>% tally() %>% mutate(prop = sprintf(n/sum(n)*100, fmt = "%#.2f"))

(coho_grp <- group_out(coho))
(chin_grp <- group_out(chin))


