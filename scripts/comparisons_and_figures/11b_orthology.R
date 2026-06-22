library(tidyverse)

# -------------------------------------------------------------------------


# potentially filter type != "gene" as well from below?


# Read and format GFF for each genome.
get_gff_map <- \(x) {
  df <- read.delim(file = x, skip = 9, sep = "\t", header = FALSE) %>% 
    filter(V3 != "gene") %>% 
    .[,c(1,3:5,9)] %>% 
    `colnames<-`(., c("chrom", "type", "start_pos", "end_pos", "info")) %>% 
    filter(type != "region") %>% 
    mutate(transcript = gsub(";.*", "", gsub(".*Parent=", "", info)),
           gene = gsub(";.*", "", gsub(".*gene=", "", info))) %>% 
    filter(!str_detect(transcript, "=")) %>% 
    dplyr::select(c(1,3:4,6:7)) %>% 
    filter(!str_detect(chrom, "^#"))
  return(df)
}

coho_gff <- get_gff_map("../genomes/coho/GCF_002021735.2_Okis_V2_genomic.gff")
chin_gff <- get_gff_map("../genomes/chinook/GCF_018296145.1/GCF_018296145.1_Otsh_v2.0_genomic.gff")

# Read output file from running orthofinder. 
orthogroups <- read.delim("data/orthology/Orthogroups.txt", sep = "", header = F) %>% 
  dplyr::rename("orthogroup" = "V1") %>% 
  mutate(orthogroup = gsub(":", "", orthogroup)) %>% 
  pivot_longer(cols = -orthogroup) %>% 
  dplyr::select(-name) %>% 
  dplyr::select(c(2,1)) %>% 
  dplyr::rename("transcript" = "value")

# Join orthogroups to genes by linking transcripts.
chg <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/RDAoutlier_snp_genes_3RDA3SD2PCA.txt") %>%
  left_join(., chin_gff[,c("gene", "transcript")], by = "gene") %>% 
  left_join(., orthogroups, by = "transcript")


cog <- read.delim("data/rdas/coho_afsGT_4bioclim_3RDA_Mahp001_n650/RDAoutlier_snp_genes_3RDA3SD2PCA.txt") %>%
  left_join(., coho_gff[,c("gene", "transcript")], by = "gene") %>% 
  left_join(., orthogroups, by = "transcript")

# Unique orthogroups and orthogroups associated with outlier SNPs in both species.
comb_orth <- data.frame(ortho = c(unique(chg$orthogroup), unique(cog$orthogroup))) %>% 
  filter(!is.na(ortho)) %>% group_by(ortho) %>% tally()
nrow(comb_orth[comb_orth$n > 1,]); nrow(comb_orth[comb_orth$n == 1,])

chinook_outlier_ortho <- chg %>% group_by(orthogroup) %>% tally()
coho_outlier_ortho <- cog %>% group_by(orthogroup) %>% tally()

shared_ortho <- comb_orth[comb_orth$n > 1,]

total_ortho <- orthogroups[orthogroups$orthogroup %in% c(chg$orthogroup),]

# If we randomly select genes from both species (same number of SNPs as there are outliers), how many
# associate with the same orthogroups by chance? Iterate this 10k times.
output <- integer()
for (i in 1:10000) {
  
  ch_sample <- sample(x = unique(orthogroups$orthogroup), size = length(unique(chg$orthogroup)))
  co_sample <- sample(x = unique(orthogroups$orthogroup), size = length(unique(cog$orthogroup)))
  
  output[i] <- table(co_sample %in% ch_sample)[2] %>% replace_na(0)
  
}


# Is the observed amount of orthogroup "sharedness" greater than what we expect by chance?
mean(output); nrow(shared_ortho)
2*pnorm(-abs((nrow(shared_ortho) - mean(output))/sd(output)))

# Visualize the above statistical test and save.
(hist <- ggplot(data = data.frame(n = unlist(output)),
                aes(x = n)) +
    geom_vline(xintercept = nrow(shared_ortho),
               linetype = 2, colour = "red2", linewidth = 1) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = 'gray90', 
                   colour = 'gray40', 
                   alpha = 2/3,
                   binwidth = 1) +
    # geom_density(linewidth = 1/2) +
    theme_bw() + 
    labs(x = "Orthogroups", 
         y = "Density") )

ggsave("plots/shared_orthogroups.tiff",
       width = 8, height = 6, dpi = 300)

