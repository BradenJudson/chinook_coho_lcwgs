library(tidyverse); library(cowplot); library(pcadapt); library(qvalue)

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2pc <- \(vcf_file) system(paste("plink.exe --vcf", vcf_file, 
                                   " --make-bed --aec --double-id --out", 
                                   gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files necessary for pcadapt.
vcf2pc("data/vcfs/chinook_lcwgs_maf005_n819_imputed.vcf.gz")
vcf2pc("data/vcfs/coho_lcwgs_maf005_n650_imputed.vcf.gz")
# vcf2pc("data/vcfs_n324/chinook_lcwgs_maf005_n324_imputed.vcf.gz")


# Function for conducting pcadapt analyses.
pcadapt2 <- \(vcf, K, q.alpha) {
  
  # Read in bed file.
  input <- read.pcadapt(gsub("\\.vcf.gz", "\\.bed", vcf), type = "bed")
  # Conduct PCadapt analysis.
  x <- pcadapt(input = input, K = K, min.maf = 0.00)
  # Output screeplot for determining best K value.
  plot(x, option = 'screeplot')
  
  # Build and format output dataframe.
  snps <- read.delim(gsub("\\.vcf.gz", "\\.bim", vcf), 
                     header = FALSE)[,c(1,4)] %>% 
    `colnames<-`(., c("CHROM", "POS")) %>% 
    mutate(pval = x$pvalues,
           chi2 = x$chi2.stat,
           data = basename(gsub("\\.vcf", "", vcf)),
           # Adjust p-values for multiple comparisons.
           qval = qvalue(pval)$qvalues,
           out = as.factor(case_when(
             pval < quantile(pval, 0.001) ~ "Y",
             TRUE ~ "N"
            )),
           loc = paste0(CHROM, "_", POS))
  
  print(summary(snps$out)); return(snps)
  
}

# Conduct pcadapt with specified parameters. Returns dataframe of all SNPs and their outlier scores/status.
chinook_imp <- pcadapt2("data/vcfs/chinook_lcwgs_maf005_n819_imputed.vcf.gz", K = 4, q.alpha = 1/1000)
write.csv(chinook_imp, "data/pca_data/chinook_imputed_n819_pcadapt_p001.csv", row.names = F)

coho_imp <- pcadapt2("data/vcfs/coho_lcwgs_maf005_n650_imputed.vcf.gz", K = 5, q.alpha = 1/1000)
write.csv(coho_imp, "data/pca_data/coho_imputed_n650_pcadapt_p001.csv", row.names = F)

genome_manhattan <- \(ds, sp) {
  
  df <- ds %>% arrange(CHROM, POS) %>% 
    mutate(cPos = cumsum(as.numeric(POS)))
  
  axisdf <- df %>% group_by(CHROM) %>% 
    summarize(center = 1/2*(max(cPos) - min(cPos)) + min(cPos)) %>% 
    mutate(LG = rownames(.)) %>% 
    filter(!LG %in% c(24, 27, 32))
  
  axisdf <- if(sp == "Chinook") axisdf[!axisdf$LG %in% c(24, 27, 32), ] else axisdf
  
  ggplot() + 
    geom_point(data = df[df$out == "N",],
               aes(x = cPos, y = -log10(pval),
                   colour = CHROM)) + 
    scale_colour_manual(values = c(rep(c("gray30", "gray90"),
                                       length(unique(df$CHROM))/2))) +
    theme_bw() + labs(x = "Position (Mbp)") +
    geom_point(data = df[df$out == "Y",],
               aes(x = cPos, y = -log10(pval)),
               colour = "red2", inherit.aes = F,
               alpha = 0.5) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          strip.text = element_text(face = "bold", 
                                    color = "black", 
                                    hjust = 0, 
                                    size = 15),
          strip.background = element_blank()) +
    scale_x_continuous(label = axisdf$LG,
                       breaks = axisdf$center) +
    labs(x = NULL) +
    ggtitle(sp)
  
  
}

(ch_manh <- genome_manhattan(chinook_imp, sp = "Chinook"))
ggsave("plots/chinook_imputed_pcadapt_pTop001.tiff", dpi = 200, width = 16, height = 6)

(co_manh <- genome_manhattan(coho_imp, sp = "Coho"))
ggsave("plots/coho_imputed_pcadapt_pTop001.tiff", dpi = 200, width = 16, height = 6)

cowplot::plot_grid(plotlist = list(
  ch_manh + ggtitle("Chinook"), co_manh + ggtitle("Coho")),
  ncol = 1)
ggsave("plots/imputed_pcadapt_pTop001_bothsp.tiff", dpi = 200, width = 16, height = 10)
