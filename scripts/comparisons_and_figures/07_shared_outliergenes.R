library(tidyverse)
set.seed(3190)

coho <- read.delim("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/RDAoutlier_snp_genes_3RDA3SD2PCA.txt")
coho_outlier_genes <- coho %>% filter(!str_detect(gene, "^LOC")) %>% group_by(gene) %>% tally()

chin <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/RDAoutlier_snp_genes_3RDA3SD2PCA.txt")
chinook_outlier_genes <- chin %>% filter(!str_detect(gene, "^LOC")) %>% group_by(gene) %>% tally()

shared_outlier_genes <- data.frame(genes = c(chinook_outlier_genes$gene, coho_outlier_genes$gene)) %>% 
  group_by(genes) %>% tally() %>% filter(n > 1)

# Read in annotation file and reformat. 
coho_genes <- read.table("../genomes/coho/GCF_002021735.2_Okis_V2_genomic.gtf", 
                         header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
  filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
  mutate(gene_id = gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
  `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info", "gene_id")) %>% 
  mutate(genegrp = case_when(grepl("^LOC[1-9]", gene_id) ~ "Non-annotated",
                             TRUE ~ "Annotated"))

chinook_genes <- read.table("../genomes/chinook/GCF_018296145.1_Otsh_v2.0_genomic.gtf", 
                            header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
  filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
  mutate(gene_id = gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
  `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info", "gene_id")) %>% 
  mutate(genegrp = case_when(grepl("^LOC[1-9]", gene_id) ~ "Non-annotated",
                             TRUE ~ "Annotated"))

shared <- rbind(
  coho_genes[coho_genes$genegrp == "Annotated",],
  chinook_genes[chinook_genes$genegrp == "Annotated",]) %>% 
  group_by(gene_id) %>% tally() %>% filter(n > 1)

output <- integer()
for (i in 1:10000) {
  
  ch_sample <- sample(x = shared$gene_id, size = nrow(chinook_outlier_genes[chinook_outlier_genes$gene %in% shared$gene_id,]))
  co_sample <- sample(x = shared$gene_id, size = nrow(coho_outlier_genes[coho_outlier_genes$gene %in% shared$gene_id,]))
  output[i] <- table(co_sample %in% ch_sample)[2]
  
}

shapiro.test(sample(output, 5e3))
hist(output,main = NULL, breaks = 50)
abline(v = nrow(shared_outlier_genes))
mean(output)
sd(output)

2*(1-pnorm(nrow(shared_outlier_genes), mean = mean(output), sd = sd(output)))

(hist <- ggplot() +
  geom_vline(xintercept = nrow(shared_outlier_genes),
             linetype = 2, colour = "red2", linewidth = 1) +
  geom_histogram(data = data.frame(n = unlist(output)), 
                 aes(x = n),
                 fill = 'gray90', colour = 'black', alpha = 2/3) +
  theme_bw() + labs(x = "Genes", y = "Count") +
  geom_boxplot(data = data.frame(n = unlist(output)),
               aes(x = n, y = 1500), inherit.aes = F,
               width = 50, outlier.alpha = 1/4))


(hist <- ggplot(data = data.frame(n = unlist(output)),
                aes(x = n)) +
    geom_vline(xintercept = nrow(shared_outlier_genes),
               linetype = 2, colour = "red2", linewidth = 1) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = 'gray90', colour = 'gray40',
                   alpha = 2/3) +
    geom_density(linewidth = 1/2) +
    theme_bw() +
    labs(x = "Genes", y = "Density"))


ggsave("plots/shared_genes_test.tiff",
       width = 8, height = 6, dpi = 300)


