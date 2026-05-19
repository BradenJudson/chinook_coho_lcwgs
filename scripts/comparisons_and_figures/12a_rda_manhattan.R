library(tidyverse); library(karyoploteR); library(GenomicRanges) 
library(ggplotify); library(cowplot)

get_chrom_info <- \(seqreport, gsp) {
  
  df <- read.delim(seqreport) %>% 
    filter(!Chromosome.name %in% c("Un", "MT")) %>% 
    dplyr::select(c("Chromosome.name", "Seq.length", "RefSeq.seq.accession")) %>% 
    mutate(Chromosome.name = gsub("LG", gsp, Chromosome.name)) %>% 
    mutate(value1 = sample(letters, nrow(.), replace = T),
           value2 = sample(letters, nrow(.), replace = T),
           start = 0) %>% 
    dplyr::select(c(Chromosome.name, start, Seq.length, value1, value2, RefSeq.seq.accession))
  
}

coho_chroms <- get_chrom_info("../genomes/coho/okis_sequence_report.tsv", "Oki")
chin_chroms <- get_chrom_info("../genomes/chinook/sequence_report.tsv",   "Ots")

out_chin_hom <- read.csv("data/homology/chin_homSNPs.csv")

out_chin <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda_loadings_manhattan_input.csv") %>%
  filter(out == "Y") %>% 
  merge(., read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/outlier_snp_bio_corrs.csv")[,c(1,7)],
        by.x = "X", by.y = "CHROM") %>% 
  mutate(POS = as.numeric(gsub(".*_", "", X)),
         CHROM = gsub(".1_.*", "\\.1", X),
         hom = case_when(X %in% out_chin_hom$X ~ "Y")) 

out_coho_hom <- read.csv("data/homology/coho_homSNPs.csv")

out_coho <- read.csv("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda_loadings_manhattan_input.csv") %>% 
  filter(out == "Y") %>% 
  mutate(POS = as.numeric(gsub(".*_", "", X))) %>% 
  mutate(CHROM = gsub(".2_.*", "\\.2", X)) %>% 
  merge(., read.csv("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/outlier_snp_bio_corrs.csv")[,c(1,2,7)],
        by = c("CHROM", "POS")) %>% 
  mutate(hom = case_when(X %in% out_coho_hom$X ~ "Y"))

# (range <- range(c(out_chin$maxld), out_coho$maxld))

kp <- plotKaryotype(chin_chroms[,c(6,2:3)], plot.type = 2)
kpPoints(kp, chr = out_chin$CHR, x = out_chin$POS, 
         y = out_chin$maxld, ymax = range[2], 
         ymin = range[1], cex = 1/3, r0 = 1/10)

plot_kar <- \(chrom_info, snp_info, species) {
  
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$topmargin <- 800
  
  range <- range(snp_info$maxld)
  
  unique_snps <- snp_info[is.na(snp_info$hom), ]
  shared_snps <- snp_info[snp_info$hom == "Y", ]
  
  kp <- plotKaryotype(chrom_info[,c(6,2:3)], plot.type = 2,
                      plot.params = pp)
  kpPoints(kp, chr = unique_snps$CHR, x = unique_snps$POS,
           y = unique_snps$maxld, ymax = range[2],
           ymin = range[1], cex = 1/3, r0 = 1/10)
  kpPoints(kp, chr = shared_snps$CHR, x = shared_snps$POS,
           y = shared_snps$maxld, ymax = range[2],
           ymin = range[1], cex = 1/3, r0 = 1/10,
           col = "red")
  kpAddMainTitle(kp, main = species, cex = 2)
  
  
}

(p1 <- as.grob(~plot_kar(chin_chroms, out_chin, "Chinook")))
(p2 <- as.grob(~plot_kar(coho_chroms, out_coho, "Coho")))

cowplot::plot_grid(plotlist = list(p1, p2))
ggsave("plots/sp_karyo_manh.tiff", dpi = 300, width = 24, height = 12, bg = 'white')




