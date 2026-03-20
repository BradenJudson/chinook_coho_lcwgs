library(tidyverse)

(hwe_files <- list.files(path = "data/heterozygosity_files/coho_lcwgs_n650_hwe/",
                         pattern = "*hwe.gz", full.names = T))

het_est <- list()

for (i in hwe_files) {
  
  pop <- gsub(".hwe.gz", "", gsub(".*/", "", i))
  
  print(pop)
  
  snp_het <- read.delim(i)
  
  het_est[[i]] <- data.frame(
    pop = pop,
    het = mean(snp_het$hetFreq)
  )
  
  rm(snp_het); gc()
  
}

pop_hets <- do.call("rbind", het_est) 
rownames(pop_hets) <- NULL

write.csv(pop_hets, "data/heterozygosity_files/coho_lcwgs_hets_n650.csv")
write.csv(pop_hets, "chin_coho_shiny/data/coho_lwcgs_hets_n83.csv")
