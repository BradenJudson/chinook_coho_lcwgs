library(tidyverse); library(cowplot)

(chinook <- readRDS("data/rdas/chinook_afsGT_4bioclim_2PCA_3RDA3SD_n819/rda12_biplot.RDS"))
(coho <- readRDS("data/rdas/coho_afsGT_4bioclim_3RDA3SD_n650/rda12_biplot.RDS"))

leg <- cowplot::get_legend(chinook)

(biplots2 <- cowplot::plot_grid(plotlist = list(
  chinook + theme(legend.position = "none") +
    ggtitle("Chinook"),
  coho + theme(legend.position = "none") +
    ggtitle("Coho"),
  leg), nrow = 1,
  rel_widths = c(1/3, 1/3, 1/6)))

ggsave("plots/chinook_coho_imputed_RDA.tiff", dpi = 300,
       height = 6, width = 14, bg = 'white')
