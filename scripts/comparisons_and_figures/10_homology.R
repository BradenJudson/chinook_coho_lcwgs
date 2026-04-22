library(tidyverse); library(GenomicRanges)

set.seed(0703)

# Read in liftover files from RLS. Positions of SNPs in their original genome, the
# other species' genome, and the positions lifted BACK to the original genome. We
# filter the ones that do not align correctly to their original chromosome.
get_lift <- \(x) readxl::read_excel(x) %>% 
  filter(!is.na(snplift_to_original_chrom)) %>% 
  dplyr::select(c("original_chr", "original_pos", "snplift_chrom", "snplift_pos")) 

chinook <- get_lift("data/orthology/ChinSNPs_COGenome_final.xlsx")
coho <- get_lift("data/orthology/CohoSNPs_CNGenome_final.xlsx")

# Outlier ranges ---------------------------------------------------------------

# Function for converting files to reduced genomic ranges objects.
# We use a +/- 10Kbp buffer on each locus and collapse overlapping regions.
toRanges <- \(df, chrom, pos) {

  gr <- GRanges(seqnames = df[[chrom]],
                ranges = IRanges(start = df[[pos]],
                                 end   = df[[pos]])) %>% 
    `names<-`(., paste0(df[[chrom]], "_", df[[pos]]))
  
  wpr <- GenomicRanges::promoters(gr,
                   upstream = 10e3,
                   downstream = 10e3,
                   use.names = TRUE)
  
  regs <- reduce(wpr)

}

chin_orig <- toRanges(chinook, "original_chr", "original_pos")
sum(width((chin_orig)))/1e6 
chin_lift <- toRanges(chinook, "snplift_chrom", "snplift_pos")
sum(width((chin_lift)))/1e6

coho_orig <- toRanges(coho, "original_chr", "original_pos")
sum(width((coho_orig)))/1e6
coho_lift <- toRanges(coho, "snplift_chrom", "snplift_pos")
sum(width((coho_lift)))/1e6


# Stats ------------------------------------------------------------------------

# Read in the full SNP list for each species.
coho_angsd <- read.delim("data/angsd_sites_gls_n650.txt", header = FALSE,
                         col.names = c("chrom", "pos", "major", "minor"))

chin_angsd <- read.delim("data/angsd_n850_sites.txt", header = FALSE,
                         col.names = c("chrom", "pos", "major", "minor"))

# Model the background distribution for the number of SNPs within homologous outlier
# regions of the opposite species genome. For example, how many Coho outlier SNPs fall
# within regions that are homologous to the regions lifted over the from Chinook dataset?
backgdist <- \(angsd, original_outliers, lift_regions) {
  
  output <- integer()
  for (i in 1:10000) {
    
    # Random subsample of SNPs.
    sub_df <- angsd[sample(nrow(angsd), nrow(original_outliers)), c("chrom", "pos")]
    sub_gr <- GRanges(seqnames = sub_df$chrom,
                      ranges = IRanges(start = sub_df$pos,
                                       end   = sub_df$pos)) %>% 
      `names<-`(., paste0(sub_df$chrom, "_", sub_df$pos))
    
    overlaps <- findOverlaps(sub_gr, lift_regions)
    
    # How many fall within this region?
    output[i] <- length(overlaps)
    
  }
  
  return(output)
  
}

# These are the null distributions.
coho_g <- backgdist(coho_angsd, coho, chin_lift)
chin_g <- backgdist(chin_angsd, chinook, coho_lift)

# How many of the above cases are actually observed?
obs_snpoverlap <- \(original, lift_range) {
  
  original_ranges <- GRanges(seqnames = original$original_chr,
                             ranges = IRanges(start = original$original_pos,
                                              end   = original$original_pos)) %>% 
    `names<-`(., paste0(original$original_chr, "_", original$original_pos))
  
  length(findOverlaps(original_ranges, lift_range))
  
}

(coho_obs <- obs_snpoverlap(coho, chin_lift))
(chin_obs <- obs_snpoverlap(chinook, coho_lift))

# Histogram and p-value for the likelihood that the observed value above is 
# greater than the null distribution by chance.
hist_p <- \(dist, obs) {
  
  hist <- ggplot(data = data.frame(n = unlist(dist)),
                 aes(x = n)) +
    geom_vline(xintercept = obs, 
               colour = 'red2', linetype = 2) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = "gray90",
                   colour = "gray40",
                   alpha = 2/3) +
    geom_density(linewidth = 1/2) +
    theme_bw() +
    labs(x = "Overlaps", y = "Density") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))
  
  print(paste0("p-value = ", 2*(1-pnorm(obs, mean = mean(dist), sd = sd(dist)))))
  
  hist
  
}

(coho_hist <- hist_p(coho_g, coho_obs))
ggsave("plots/coho_homology_overlap_hist_10kbp.tiff", dpi = 300, width = 8, height = 6)
(chin_hist <- hist_p(chin_g, chin_obs))
ggsave("plots/chin_homology_overlap_hist_10kbp.tiff", dpi = 300, width = 8, height = 6)

