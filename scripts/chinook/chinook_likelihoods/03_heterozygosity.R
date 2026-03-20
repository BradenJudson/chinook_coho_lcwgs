library(tidyverse)

(hwe_files <- list.files(path = "data/heterozygosity_files/chinook_lcwgs_n819_hwe/",
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

write.csv(pop_hets, "data/heterozygosity_files/chinook_lcwgs_hets_n819.csv")
