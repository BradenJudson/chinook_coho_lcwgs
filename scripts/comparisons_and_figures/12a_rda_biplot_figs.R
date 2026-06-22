library(tidyverse); library(cowplot)

# -------------------------------------------------------------------------

# Read in RDA biplot objects and format together with consistent aspect ratios and a shared legend.
(chinook <- readRDS("data/rdas/chinook_afsGT_scaled4bio_3RDA_topMahp001_n819_4bioGF/rda12_biplot.RDS"))
(coho <- readRDS("data/rdas/coho_afsGT_scaled4bio_3RDA_topMahp001_n650_4bioGF/rda12_biplot.RDS"))

leg <- cowplot::get_legend(chinook)

(biplots2 <- cowplot::plot_grid(plotlist = list(
  chinook + theme(legend.position = "none") +
    ggtitle("Chinook"),
  coho + theme(legend.position = "none") +
    ggtitle("Coho"),
  leg), nrow = 1,
  rel_widths = c(1/3, 1/3, 1/6)))

ggsave("plots/chinook_coho_imputed_RDA_Mahp001.tiff", dpi = 300,
       height = 6, width = 14, bg = 'white')

