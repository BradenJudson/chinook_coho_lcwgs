
library(tidyverse); library(scales); library(plotly)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "chinook") %>% 
  group_by(region_revised) %>% 
  mutate(mLat = mean(latitude)) %>% 
  arrange(mLat)

# Arrange lineages by average latitude.
lin <- sites %>% group_by(Lineage) %>%
  summarize(mLat = mean(Latitude)) %>%
  arrange(mLat)

sites <- merge(sites, lin) %>%
  arrange(mLat, Latitude) %>%
  mutate(Lineage = factor(Lineage, levels = lin$Lineage)) %>%
  rownames_to_column("sitenum") %>% 
  mutate(sitenum = as.numeric(sitenum))

# # Dataframe for manual fill/alpha values per lineage
# coldf <- data.frame(
#   Lineage = as.character(levels(sites$Lineage)),
#   alpha = c(rep(c(0.3, 1), length(unique(sites$Lineage))/2), 0.3),
#   colour= hue_pal()(length(levels(sites$Lineage)))
# )
# 
# coldf[11:15,3] <- coldf[15:11,3]
# 


# indvs <- readxl::read_excel("data/chinook_lcwgs_samples.xlsx") %>% 
#   filter(!fate %in% c("Discard (coverage too low)", "PCA outlier", "Drop population"))
# indvs$site <- gsub("_OT", "", gsub("_ST", "", indvs$site))



indvs <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "chinook") %>% filter(is.na(fate))

dat <- merge(indvs, sites, by.x = "population", by.y = "site") 

# sites[!sites$Site %in% dat$site, "site"]
# unique(indvs[!indvs$site %in% sites$Site, "site"])

lcwgs_pca <- \(cov_mat, bam_list, PCX, PCY, scree_pos) {
  
  # Read in file/sample names. These must be in the same order as the bam 
  # files used during the PCA calculations!
  # Requires some minor formatting.
  files <- read.table({{bam_list}}, col.names = c("file")) %>% 
    mutate(id = gsub("Kand0309-", "Kand0309\\.", gsub("\\..*$","",
                gsub(".*IDT_i5_[0-9]{1,3}.", "", gsub("(.*/)", "", file)))))
  
  cov    <- read.table({{cov_mat}})               # Read in covariance matrix.
  mmepca <- eigen(cov)                            # Extract eigenvalues.
  eigvec <- as.data.frame(mmepca$vectors) %>%     # Extract eigenvectors.
    `colnames<-`(., gsub("V", "PC", colnames(.))) # Rename columns.
  
  # % variation explained by each PC axis.
  pca.eigenval.sum <- sum(mmepca$values)
  PC1 <- sprintf("%0.1f", mmepca$values[1]/pca.eigenval.sum*100)
  PC2 <- sprintf("%0.1f", mmepca$values[2]/pca.eigenval.sum*100)
  
  scree <- data.frame(PC = paste0("PC", c(1:8)),
                      EV = mmepca$values[1:8]/pca.eigenval.sum*100) %>% 
    mutate(sel = case_when( PC %in% c({{PCY}}, {{PCX}}) ~ "Y",
                           !PC %in% c({{PCY}}, {{PCX}}) ~ "N"),
           PC = fct_reorder(PC, desc(EV))) %>% 
    ggplot(data = ., aes(x = PC, y = EV, fill = sel)) +
    scale_fill_manual(values = c("gray90", "gray20")) +
    geom_bar(stat = 'identity', colour = "black") + theme_bw() + 
    theme(legend.position = "none", rect = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = NULL, y = "Variance explained (%)")
  
    
  pca_dat <- cbind(files, eigvec) %>% # Bind it all together.
    merge(., dat, by = "id")   # Add populations.
  
  pca_avgs <- pca_dat %>% group_by(population) %>%
    summarize(mPCX = mean(!!sym(PCX)),
              mPCY = mean(!!sym(PCY)))

  pca_dat <- merge(pca_dat, pca_avgs, by = "population") 
    # mutate(sitenum = as.numeric(sitenum)) %>%
    # arrange(sitenum) %>%
    # mutate(region_revised = factor(region_revised, levels = levels(sites$region_revised), ordered = T))


  # Visualize PCA consistent with above "vcf-based" PCAs.
  pc_scatter <- ggplot(data  = pca_dat,
                        aes(x = !!sym(PCX),
                            y = !!sym(PCY),
                            fill  = region_revised)) +
      geom_segment(aes(x = mPCX, y = mPCY,
                       xend = !!sym(PCX), yend = !!sym(PCY),
                       colour = region_revised),
                   show.legend = FALSE) +
      # scale_color_manual(values = coldf$colour) +
      geom_point(size = 2, shape = 21,aes(fill = region_revised)) + theme_bw() +
      # scale_fill_manual(values = alpha(c(coldf$colour), coldf$alpha)) +
      guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
      theme(legend.key = element_rect(fill = NA, color = NA)) +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE),
           shape = guide_legend(ncol =1, reverse = T)) +
      labs(x = paste0({{PCX}}, " (", sprintf("%0.1f", scree$data[scree$data$PC == {{PCX}}, "EV"]), "%)"),
           y = paste0({{PCY}}, " (", sprintf("%0.1f", scree$data[scree$data$PC == {{PCY}}, "EV"]), "%)"))

}

(pop_pca12 <- lcwgs_pca(cov_mat = "data/pca_data/chinook_angsd_n819.cov",
                        bam_list = "data/chinook_bam_list_n819.txt",
                        PCX = "PC1", PCY = "PC2", scree_pos = "br") +
    scale_x_continuous(transform = "reverse"))

pca_dat <- pop_pca12@data[,c("id", "population", "region_revised", "mLat", paste0("PC", c(1:20)))] 
write.csv(pca_dat, "data/pca_data_pc12.csv", row.names = F)
write.csv(pca_dat, "chin_coho_shiny/data/chinook_likelihoods_pca_scores_n819.csv", row.names = F)

ggsave("plots/pca12.tiff", dpi = 300, width = 10, height = 8)

pop_pcs <- pca_dat %>% group_by(Pop) %>% 
  summarize(PC1 = round(mean(PC1), 5),
            PC2 = round(mean(PC2), 5))

write.csv(pop_pcs, "data/pop_pcas_n106.csv", row.names = F)

