library(tidyverse); library(scales)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  filter(species == "coho") %>% 
  mutate(site = gsub("FishingBranchWeir", "FishingBranch", site)) %>% 
  group_by(region_revised) %>% mutate(mLat = mean(latitude))
indvs <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(species == "coho") %>% filter(is.na(fate))

dat <- merge(indvs, sites, by.x = "population", by.y = "site")

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
    # mutate(Lineage = factor(Lineage, levels = levels(sites$Lineage), ordered = T))
  
  
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
    geom_point(size = 2, shape = 21, aes(fill = region_revised)) + theme_bw() +
    # scale_fill_manual(values = alpha(c(coldf$colour), coldf$alpha)) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(legend.key = element_rect(fill = NA, color = NA),
          legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE),
           shape = guide_legend(ncol =1, reverse = T)) +
    labs(x = paste0({{PCX}}, " (", sprintf("%0.1f", scree$data[scree$data$PC == {{PCX}}, "EV"]), "%)"),
         y = paste0({{PCY}}, " (", sprintf("%0.1f", scree$data[scree$data$PC == {{PCY}}, "EV"]), "%)"))
  
}

(pop_pca12 <- lcwgs_pca(cov_mat = "data/pca_data/coho_angsd_n650.cov",
                  bam_list = "data/coho_bam_list_n650.txt",
                  PCX = "PC1", PCY = "PC2", scree_pos = "br"))

ggsave("plots/coho_PCA12_GLs_n650.tiff", dpi = 300,
       width = 12, height = 9)

pca_dat <- pop_pca12@data[,c("id", "population", "region_revised", "mLat", paste0("PC", c(1:20)))] 
write.csv(pca_dat, "data/pca_data_n650_likelihoods.csv", row.names = F)
write.csv(pca_dat, "chin_coho_shiny/data/coho_likelihoods_pca_scores_n650.csv", row.names = F)

pop_pcas <- pca_dat %>% group_by(site) %>% 
  summarize(PC1 = mean(PC1), PC2 = mean(PC2),
            PC3 = mean(PC3), PC4 = mean(PC4))

write.csv(pop_pcas, "data/PCA/pop_pcas_afsGL_n650.csv", row.names = FALSE)

pc_sd <- pca_dat %>% group_by(site) %>% 
  summarise(sdPC1 = sd(PC1)) %>% 
  mutate(site = fct_reorder(site, sdPC1, mean,.desc = T))

ggplot(data = pc_sd, 
       aes(x = site, y = sdPC1)) + 
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/2)) +
  labs(x = NULL)



(pc_bounds <- pca_dat[,c("id", "site", "PC1", "PC2")] %>% 
  pivot_longer(cols = c("PC1", "PC2")) %>% 
  group_by(site, name) %>%
    summarize(
      meanPC = mean(value),
      sdPC   = sd(value),
      upper = meanPC + 2*sdPC,
      lower = meanPC - 2*sdPC)
  )
  
  
pc_out <- pca_dat[,c("id", "site", "PC1", "PC2")] %>% 
  pivot_longer(cols = c("PC1", "PC2")) %>% 
  merge(., pc_bounds, by = c("site", "name")) %>% 
  mutate(out = case_when(value > upper | value < lower ~ "Y"))

(pcsdplot <- ggplot() +
  geom_point(data = pc_out,
             aes(x = site,
                 y = value,
                 colour = out)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/2)) +
  labs(x = NULL) +
    facet_wrap(~name, ncol = 1, scales = "free_y"))
  
ggplotly(pcsdplot)





ggplotly(ggplot(data = imp_pca,
       aes(x = PC1, 
           y = PC2,
           fill = site)) +
  geom_point(aes(text = id),
             shape = 21) +
  stat_ellipse(level = 0.9))
  
