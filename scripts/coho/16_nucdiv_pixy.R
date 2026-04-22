library(tidyverse)

# Generate nucleotide diversity from full VCF 
# e.g., pixy --stats pi --vcf $VCF --populations $POPMAP --window_size 10000
read_pi <- \(x) read.delim(x) %>% 
  group_by(pop) %>% summarize(meanPi = mean(avg_pi, na.rm = TRUE)) 

co_pi <- read_pi("data/pi/coho_n650_pixy_pi.txt")

summary(co_pi$meanPi)

write.csv(ch_pi, "data/pi/coho_imputed_pi_n83.csv", row.names = FALSE)
write.csv(ch_pi, "chin_coho_shiny/data/coho_imputed_pi_n83.csv", row.names = F)
