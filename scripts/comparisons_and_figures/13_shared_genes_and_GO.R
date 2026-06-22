library(tidyverse); library(topGO); library(GenomicRanges); library(eulerr)

chinGO <- read.csv("data/outlier_genes/chinook_afsGT_goterms.csv")
cohoGO <- read.csv("data/outlier_genes/coho_afsGT_goterms.csv")

# Genes ------------------------------------------------------------------------

coho <- read.delim("data/rdas/coho_afsGT_4bioclim_3RDA_Mahp001_n650/RDAoutlier_snp_genes_3RDA3SD2PCA.txt") %>% 
  mutate(SNP = paste0(CHROM, "_", POS)) 

chin <- read.delim("data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/RDAoutlier_snp_genes_3RDA3SD2PCA.txt") %>% 
  mutate(SNP = paste0(CHROM, "_", POS)) 

toRanges <- \(df) {
  
  gr <- GRanges(seqnames = df$CHROM,
                ranges = IRanges(start = df$POS,
                                 end   = df$POS)) %>% 
    `names<-`(., paste0(df$CHROM, "_", df$POS))
  
  wpr <- GenomicRanges::promoters(gr,
                                  upstream   = 10e3,
                                  downstream = 10e3,
                                  use.names  = TRUE)
  
  
}

chin_regs_o <- toRanges(chin)
coho_regs_o <- toRanges(coho)


gtf2genes <- \(gtf) {
  
  df <- read.table(gtf, header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
    filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
    mutate(gene_id = gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
    `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info", "gene_id")) %>% 
    mutate(genegrp = case_when(grepl("^LOC[1-9]", gene_id) ~ "Non-annotated", TRUE ~ "Annotated"))

  gene_ranges <- GRanges(seqnames = df$chrom,
                         ranges = IRanges(start = df$start_pos,
                                          end   = df$end_pos)) %>% 
    `names<-`(., df$gene_id)
  
  }

# Read in annotation file and reformat. 
coho_genes <- gtf2genes("../genomes/coho/GCF_002021735.2_Okis_V2_genomic.gtf")
chin_genes <- gtf2genes("../genomes/chinook/GCF_018296145.1_Otsh_v2.0_genomic.gtf")

snps2genesDF <- \(snps, gtf) {
  
  overlaps <- findOverlaps(snps, gtf)
  out_snp_genes <- cbind.data.frame(
    names(snps)[queryHits(overlaps)],
    names(gtf)[subjectHits(overlaps)]) %>% 
    `colnames<-`(., c("snp", "gene")) %>% 
    mutate(CHR = as.factor(gsub("^([^_]*_[^_]*).*", "\\1", snp)),
           POS = as.numeric(gsub(".*.[0-9]{1}_", "", snp)))
  
}

chin_o_genes <- snps2genesDF(chin_regs_o, chin_genes)
coho_o_genes <- snps2genesDF(coho_regs_o, coho_genes)

shared_o <- data.frame(gene = c(unique(chin_o_genes$gene), 
                                unique(coho_o_genes$gene))) %>% 
  group_by(gene) %>% tally() %>% filter(n > 1)

png("plots/outlier_overlap.png", res = 300, width = 12, height = 8, units = "in")
plot(euler(list(
  "Coho (all)" = unique(c(coho_o_genes$gene)),
  "Chinook (all)" = unique(c(chin_o_genes$gene))
)))
dev.off()

# GO ----------------------------------------------------------------------

ch_gomap <- readMappings("../genomes/chinook/go_map.txt")
ch_ref_genes <- names(ch_gomap)

co_gomap <- readMappings("../genomes/coho/go_map.txt")
co_ref_genes <- names(co_gomap)

doGO <- \(genes, gomap) {
  
    out_genes <- genes[genes %in% names(gomap)]
    full_vect <- factor(as.integer(names(gomap) %in% out_genes))
    names(full_vect) <- names(gomap)
  
    GOdata <- new("topGOdata",
                  description = names(gomap),
                  ontology = "BP",
                  allGenes = full_vect,
                  nodeSize = 5,
                  annot = annFUN.gene2GO,
                  gene2GO = gomap)
    
    fishers <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    
    # Assemble GO stats, etc. into a dataframe with numeric p values.
    GO_results <- as.data.frame(GenTable(GOdata, fishers, 
                                         topNodes = length(fishers@score),
                                         numChar = 1e6)) %>% 
      dplyr::rename(pvalues = result1)
    GO_results$pvalues <- as.numeric(GO_results$pvalues)
    
    return(GO_results)

}

chin_o_go <- doGO(chin_o_genes$gene, ch_gomap)
coho_o_go <- doGO(coho_o_genes$gene, co_gomap)


