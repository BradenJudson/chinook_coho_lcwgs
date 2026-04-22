library(tidyverse); library(ggVennDiagram)

coho <- read.delim("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/RDAoutlier_snp_genes_3RDA3SD2PCA.txt")
chin <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/RDAoutlier_snp_genes_3RDA3SD2PCA.txt")

coho_genes <- unique(coho$gene)
chin_genes <- unique(chin$gene)

out_genes <- do.call("rbind", list(
  data.frame(gene = unique(coho$gene), sp = "coho"),
  data.frame(gene = unique(chin$gene), sp = "chin")
)) 

out_genes %>% group_by(gene) %>% tally() %>% filter(n > 1) %>% .[,1] -> shared

outgene_list <- list(coho = unique(coho$gene),
                     chin = unique(chin$gene))

ggVennDiagram(outgene_list, label_alpha = 0,
              category.names = c("Coho", "Chinook")) +
  scale_fill_gradient(low = "white", high = "gray90") +
  theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  coord_flip()

ggsave("plots/supp/rdaoutGT_venn_genes.tiff", dpi = 250, width = 6, height = 6, bg = 'white')


chin_chr <- read.delim("../genomes/chinook/sequence_report.tsv") %>% 
  filter(Molecule.type == "Linkage Group")

coho_chr <- read.delim("../genomes/coho/okis_sequence_report.tsv") %>% 
  filter(Molecule.type == "Linkage Group")

(chin_manh <- ggplot() +
  geom_segment(data = chin_chr, 
               aes(x = RefSeq.seq.accession, y = 0,
                   xend = RefSeq.seq.accession, yend = Seq.length),
               linewidth = 1) +
  # geom_point(data = chin[!chin$gene %in% shared$gene,],
  #            aes(x = CHROM, y = POS)) +
  geom_point(data = chin[chin$gene %in% shared$gene,],
             aes(x = CHROM, y = POS),
             shape = 21, colour = "red2", fill = "red3") +
  theme_bw() + ggtitle("Chinook") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1/2)) +
  labs(y = "Position", x = NULL) )

(coho_manh <- ggplot() +
  geom_segment(data = coho_chr, 
               aes(x = RefSeq.seq.accession, y = 0,
                   xend = RefSeq.seq.accession, yend = Seq.length),
               linewidth = 1) +
  geom_point(data = coho[coho$gene %in% shared$gene,],
             aes(x = CHROM, y = POS),
             shape = 21, colour = "red2", fill = "red3") +
  theme_bw() + ggtitle("Coho") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1/2)) +
  labs(y = "Position", x = NULL) )

cowplot::plot_grid(plotlist = list(chin_manh, coho_manh), ncol = 1)
ggsave("plots/supp/comp_outliers.tiff", dpi = 300, width = 12, height = 10)


# -------------------------------------------------------------------------

ch_shared_snp <- chin[chin$gene %in% shared$gene,]
co_shared_snp <- coho[coho$gene %in% shared$gene,]

# -------------------------------------------------------------------------
# Read in annotation file and reformat. 
coho_genes <- read.table("../genomes/coho/GCF_002021735.2_Okis_V2_genomic.gtf", 
                         header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
  filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
  `rownames<-`(., gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
  `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info")) %>% 
  rownames_to_column("gene_id") %>% 
  mutate(genegrp = case_when(grepl("^LOC[1-9]", gene_id) ~ "Non-annotated"))


chinook_genes <- read.table("../genomes/chinook/GCF_018296145.1_Otsh_v2.0_genomic.gtf", 
                            header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
  filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
  `rownames<-`(., gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
  `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info")) %>% 
  rownames_to_column("gene_id") %>% 
  mutate(genegrp = case_when(grepl("^LOC[1-9]", gene_id) ~ "Non-annotated"))

common_genes <- rbind(
  coho_genes[is.na(coho_genes$genegrp), c(1:3)],
  chinook_genes[is.na(chinook_genes$genegrp), c(1:3)]) %>% 
  group_by(gene_id) %>% tally() %>% filter(n > 1)

nano <- rbind(
  coho_genes[coho_genes$genegrp == "Non-annotated",],
  chinook_genes[chinook_genes$genegrp == "Non-annotated",]) %>% 
  group_by(gene_id) %>% tally() %>% filter(n > 1)


chinook_outlier_genes <- chin %>% filter(!str_detect(gene, "^LOC")) %>% 
  group_by(gene) %>% tally()
coho_outlier_genes <- coho %>% filter(!str_detect(gene, "^LOC")) %>% 
  group_by(gene) %>% tally()


# x is number of possible genes
output <- integer()
for (i in 1:10000) {
  
  ch_sample <- sample(x = common_genes$gene_id, size = nrow(chinook_outlier_genes))
  co_sample <- sample(x = common_genes$gene_id, size = nrow(coho_outlier_genes))
  output[i] <- table(co_sample %in% ch_sample)[2]
  
}

output[is.na(output)] <- 0
# shapiro.test(output)
hist(output,main = NULL)
mean(output)
sd(output)


# -------------------------------------------------------------------------

coho_go <- read.csv("data/outlier_genes/coho_afsGT_goterms.csv")
chin_go <- read.csv("data/outlier_genes/chinook_afsGT_goterms.csv")

shared <- data.frame(go = c(coho_go[coho_go$pvalues < 0.05, "GO.ID"],
                            chin_go[chin_go$pvalues < 0.05, "GO.ID"])) %>% 
  group_by(go) %>% tally()




