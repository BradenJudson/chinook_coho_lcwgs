library(tidyverse); library(topGO); library(GenomicRanges)


# -------------------------------------------------------------------------

# Define Chinook GO map and make a vector of annotated genes.
gomap <- readMappings("../genomes/chinook/go_map.txt")
ref_genes <- names(gomap)


# Read in annotation file and reformat. 
gtf <- read.table("../genomes/chinook/GCF_018296145.1_Otsh_v2.0_genomic.gtf", 
                  header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% 
  filter(V3 == "gene") %>% .[ , c(1, 4, 5, 7, 9)] %>% 
  `rownames<-`(., gsub("gene_id ", "", gsub(";.*$", "", .$V9))) %>% 
  `colnames<-`(., c("chrom", "start_pos", "end_pos", "strand", "gene_info"))

# Define regions associated with genes.
# Include strand information. 
genegr <- GRanges(seqnames = gtf$chrom, 
                  ranges = IRanges(start = gtf$start_pos, 
                                   end = gtf$end_pos)) %>% 
  `names<-`(., rownames(gtf))
strand(genegr) <- gtf$strand

# Expand gene regions to include up/downstream flanking sequence.
gwpgr <- promoters(genegr, 
                   upstream = 10e3, 
                   downstream = 10e3, 
                   use.names = TRUE)


# -------------------------------------------------------------------------


rda_outliers <- read.csv("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/outlier_snp_bio_corrs.csv") %>% 
  mutate( POS = as.numeric(gsub("*.*_", "", CHROM)),
          CHROM = gsub(".1_.*", ".1", CHROM))

# Convert outlier SNPs to "regions" (only 1b wide). 
# Name as unique loci.
outlier_snps <- GRanges(seqnames = rda_outliers$CHROM,
                        ranges = IRanges(start = rda_outliers$POS,
                                         end   = rda_outliers$POS)) %>% 
  `names<-`(., paste0(rda_outliers$CHROM, "_", rda_outliers$POS))

# Find overlaps between outlier SNPs and annotated genes.
outlier_overlaps <- findOverlaps(outlier_snps, gwpgr)
out_snp_genes <- cbind.data.frame(
  names(outlier_snps)[queryHits(outlier_overlaps)],
  names(gwpgr)[subjectHits(outlier_overlaps)]
) %>% `colnames<-`(., c("snp", "gene")) %>% 
  mutate( POS = as.numeric(gsub("*.*_", "", snp)),
          CHROM = gsub(".1_.*", ".1", snp))

write.table(out_snp_genes[,c(4,3,2)], "data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/RDAoutlier_snp_genes_3RDA3SD2PCA.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

# -------------------------------------------------------------------------

# Identify outlier-associated genes and create vector of detections/refs.
out_genes <- out_snp_genes$gene[out_snp_genes$gene %in% ref_genes]
full_vect <- factor(as.integer(ref_genes %in% out_genes))
names(full_vect) <- ref_genes

# Create GO object with terms above.
GOdata <- new("topGOdata",
              description = out_genes,
              ontology = "BP",
              allGenes = full_vect,
              nodeSize = 5,
              annot = annFUN.gene2GO,
              gene2GO = gomap)

# Test significance of GO terms. 
fishers <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

# Assemble GO stats, etc. into a dataframe.
# Reformat p-values and adjust for multiple comparisons.
GO_results <- as.data.frame(GenTable(GOdata, fishers, 
                                     topNodes = length(fishers@score),
                                     numChar = 1e6)) %>% 
  dplyr::rename(pvalues = result1)
GO_results$pvalues <- as.numeric(GO_results$pvalues)

write.csv(GO_results, "data/outlier_genes/chinook_afsGT_goterms.csv", row.names = F)

