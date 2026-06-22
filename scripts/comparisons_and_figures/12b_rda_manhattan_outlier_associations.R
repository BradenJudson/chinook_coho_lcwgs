library(tidyverse); library(cowplot)

# -------------------------------------------------------------------------

# Reformat SNP p-values and associations into a dataframe ready to plot w/ ggplot2.
make_manhattan_df <- \(loadings, seqreport, outlier_associations) {
  
  snp_loadings <- read.csv(loadings) 
  # Below separates CHROM from POS based on the position of the underscores.
  snp_loadings$CHROM <- as.factor(gsub("^([^_]*_[^_]*).*", "\\1", snp_loadings$X))
  snp_loadings$POS <- as.numeric(gsub(".*.[0-9]{1}_", "", snp_loadings$X))
  
  df <- snp_loadings[,c("CHROM", "POS", "p.values", "q.values")]
  
  chrom_info <- read_tsv(seqreport) %>% 
    dplyr::select(c(`RefSeq seq accession`, `Chromosome name`)) %>% 
    `colnames<-`(., c("CHROM", "LG")) %>% filter(!LG %in% c("Un", "MT")) 
  chrom_info <- chrom_info[order(chrom_info$LG),]
  
  rda_snp_cors <- read.csv(outlier_associations) 
  
  # Colours chromosomes dark/light for better visuals.
  (ev_chr <- c(unlist(chrom_info[, 2])[c(FALSE, TRUE)]))
  (od_chr <- c(unlist(chrom_info[, 2])[c(TRUE, FALSE)]))
  
  rda_df <- left_join(df, chrom_info, by = 'CHROM') %>%
    mutate(POS = as.numeric(POS)) %>%
    left_join(., rda_snp_cors,
              by = c("CHROM", "POS")) %>% 
    mutate(noutfill = case_when(is.na(association) & LG %in% ev_chr ~ "EvNon",
                                is.na(association) & LG %in% od_chr ~ "OdNon")) %>%
    mutate(association = replace_na(association, "Non-outlier")) %>%
    mutate(association = tools::toTitleCase(association)) %>%
    arrange(CHROM, POS) %>%
    mutate(cPos = cumsum(as.numeric(POS)))

  # Intentionally order factor levels.
  rda_df$association <- factor(rda_df$association,
                               levels = c("Unclassified","Bio1", "Bio5", "Bio12", "Bio15",
                                          "Temperature", "Precipitation",
                                          "Non-outlier"))

  return(rda_df)
  
  
  
}


# Assemble dataframe for both species. 
chinook_df <- make_manhattan_df(loadings = "data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/rda_mah_rdadapt_k3.csv",
                                seqreport = "../genomes/chinook/sequence_report.tsv",
                                outlier_associations = "data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/outlier_snp_corrs.csv")
write.csv(chinook_df, "data/rdas/chinook_afsGT_4bioclim_2PCA_3RD_Mahp001_n819/rda_manhattan_df.csv")



coho_df <- make_manhattan_df(loadings = "data/rdas/coho_afsGT_4bioclim_3RDA_Mahp001_n650/coho_rda_mah_k3_rdadapt.csv",
                             seqreport = "../genomes/coho/okis_sequence_report.tsv",
                             outlier_associations = "data/rdas/coho_afsGT_4bioclim_3RDA_Mahp001_n650/outlier_snp_bio_corrs.csv")
write.csv(coho_df, "data/rdas/coho_afsGT_4bioclim_3RDA_Mahp001_n650/rda_manhattan_df.csv")

# Outlier colours.
outCols <- c(
  "Bio1" = 'goldenrod1', "Bio5" = 'firebrick3',
  "Bio12" = "blue1", "Bio15" = "forestgreen",
  "Non-outier" = 'gray96', "Unclassified" = 'gray15',
  "Precipitation" = "turquoise3", "Temperature" = 'orange1')

# Dummy plot for stealing the legend only.
legend_plot <- ggplot() +
  geom_point(data = coho_df[is.na(coho_df$noutfill),],
             aes(x = cPos, 
                 y = -log10(p.values),
                 colour = association),
             alpha = 0.8) + theme_bw() +
  scale_colour_manual(values = outCols) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text  = element_text(size = 14)) +
  guides(color = guide_legend(nrow = 1))

# Custom function for turning RDA p-values into a Manhattan plot.
rda_manhattan <- \(rda_df, plot_title) {
  
  rda_df <- rda_df[rda_df$p.values < 0.05,]
  
  # Set up a DF for the axis values, so labels are positioned nicely.
  axis_vals <- rda_df %>% 
    group_by(CHROM) %>% 
    summarize(center = 1/2*(max(cPos) - min(cPos)) + min(cPos),
              min = min(cPos, na.rm = TRUE)) %>% 
    mutate(LG = rownames(.)) 
  
  (maxRDA_manh <- ggplot() +
      # First plot non-outlier SNPs.
      geom_point(data = rda_df[!is.na(rda_df$noutfill),],
                 aes(x = cPos, 
                     y = -log10(p.values),
                     colour = noutfill),
                 alpha = 0.8,
                 show.legend = FALSE) +
      scale_colour_manual(values = c("gray80", "gray90")) +
      theme_bw() +
      ggnewscale::new_scale_colour() +
      # Plot unclassified SNPs second because they are darker - helps see associations.
      geom_point(data = rda_df[rda_df$association == "Unclassified",],
                 aes(x = cPos, 
                     y = -log10(p.values),
                     colour = association),
                 alpha = 0.8) +
      scale_colour_manual(values = outCols) +
      ggnewscale::new_scale_colour() +
      # Plot everything else.
      geom_point(data = rda_df[rda_df$association %in% c("Precipitation", "Bio5", "Bio12", "Temperature", "Bio15", "Bio1"),],
                 aes(x = cPos, y = -log10(p.values), colour = association), alpha = 0.8) +
      scale_colour_manual(values = outCols) +
      theme(legend.position = "none",
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.text  = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16),
            panel.grid.minor.x = element_line(size = 0.25, 
                                              colour = "gray25",
                                              linetype = 2)) +
      # Specify axis values from DF created earlier.
      scale_x_continuous(label = axis_vals$LG,
                         breaks = axis_vals$center, 
                         minor_breaks = axis_vals$min,
                         guide = guide_axis(minor.ticks = F),
                         expand = c(0.01, 0.01)) +
      labs(x = NULL, y = expression(-log[10](p-value)))) +
    guides(color = guide_legend(nrow = 1)) 
  
  if(!is.na(plot_title)) {
    maxRDA_manh <- maxRDA_manh + ggtitle(plot_title)
  }
  
  return(maxRDA_manh)
  
  
}

# Below, plot each Manhattan independently and together. 

(ch_manh <- cowplot::plot_grid(plotlist = list(
  cowplot::get_legend(legend_plot),
  rda_manhattan(chinook_df, plot_title = NA)), 
  ncol = 1, rel_heights = c(0.1, 1)))
ggsave("plots/chinook_RDA_colorManhattan_scaled.tiff", dpi = 300,
       height = 8, width = 24, bg = 'white')


(co_manh <- cowplot::plot_grid(plotlist = list(
  cowplot::get_legend(legend_plot),
  rda_manhattan(coho_df, plot_title = NA)), 
  ncol = 1, rel_heights = c(0.1, 1)))
ggsave("plots/coho_RDA_colorManhattan_scaled.tiff", dpi = 300,
       height = 8, width = 24, bg = 'white')

rm(ch_manh); rm(co_manh); gc()

cowplot::plot_grid(plotlist = list(
  cowplot::get_legend(legend_plot),
  rda_manhattan(chinook_df, plot_title = "Chinook"),
  rda_manhattan(coho_df, plot_title = "Coho")
), ncol = 1, rel_heights = c(0.1, 1, 1))

ggsave("plots/rda_manhattans_assoc_scaled.tiff", dpi = 300,
       height = 12, width = 24, bg = 'white')


