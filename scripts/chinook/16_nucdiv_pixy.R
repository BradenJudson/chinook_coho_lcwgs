library(tidyverse)


# Generate nucleotide diversity from full VCF 
# e.g., pixy --stats pi --vcf $VCF --populations $POPMAP --window_size 10000
read_pi <- \(x) read.delim(x) %>% 
  group_by(pop) %>% summarize(meanPi = mean(avg_pi, na.rm = TRUE)) 

ch_pi <- read_pi("data/pi/chinook_n819_pixy_pi_pi.txt")

summary(ch_pi$meanPi)

write.csv(ch_pi, "data/pi/chinook_imputed_pi_n106.csv", row.names = FALSE)
write.csv(ch_pi, "chin_coho_shiny/data/chinook_imputed_pi_n106.csv", row.names = F)
