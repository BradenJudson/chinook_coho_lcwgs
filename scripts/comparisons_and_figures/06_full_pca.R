library(tidyverse); library(cowplot)

ch <- read.csv("data/pca_data/chinook_imputed_pca_scores_n819.csv") %>%
  mutate(sp = "Chinook", dt = "Imputed")

co <- read.csv("data/pca_data/coho_pca_data_n650_imputed.csv") %>%
  mutate(sp = "Coho", dt = "Imputed")

(chin_imp <- readRDS("data/pca_data/chinook_imputed_n819_PCAobj.RDS"))
(coho_imp <- readRDS("data/pca_data/coho_imputed_n650_PCAobj.RDS"))

leg1 <- cowplot::get_legend(chin_imp)

(comb <- cowplot::plot_grid(plotlist = list(
  chin_imp + theme(legend.position = "none", 
                   plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")) +
    ggtitle( "Chinook"), 
  coho_imp + theme(legend.position = "none", 
                   plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")) +
    ggtitle("Coho"), leg1),
  nrow = 1, rel_widths = c(0.4, 0.4, 0.2)))

ggsave("plots/chinook_coho_imputed_pca.tiff", dpi = 300,
       height = 6, width = 16, bg = 'white')
