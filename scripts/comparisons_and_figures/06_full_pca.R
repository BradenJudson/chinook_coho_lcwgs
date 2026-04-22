library(tidyverse); library(cowplot); library(scales)

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
  nrow = 1, rel_widths = c(0.4, 0.4, 0.18)))

ggsave("plots/chinook_coho_imputed_pca.tiff", dpi = 300,
       height = 6, width = 16, bg = 'white')


# -------------------------------------------------------------------------

chsc <- read.delim("data/pca_data/chinook_imputed_n819.eigenval", 
                 col.names = "eigv", header = F) %>% 
  mutate(pvar = eigv/sum(eigv), sp = "Chinook",
         PC   = as.numeric(rownames(.)))

cosc <- read.delim("data/pca_data/coho_imputed_pca_n650.eigenval",
                   col.names = "eigv", header = F) %>% 
  mutate(pvar = eigv/sum(eigv), sp = "Coho",
         PC   = as.numeric(rownames(.)))

ggplot(data = rbind(chsc, cosc),
       aes(x = PC, y = pvar)) +
  geom_line() +
  geom_point() +
  facet_wrap(~sp, ncol = 1, scales = "free_y") +
  theme_bw() +
  scale_x_continuous(limits = c(1, 10),breaks = c(1:10),
                     labels = label_number(accuracy = 1)) +
  labs(x = "PC", y = "Percentage of explained variance") +
  scale_y_continuous(labels = scales::percent)

ggsave("plots/PCA_screes_imputed.tiff", dpi = 300,
       width = 8, height = 6)
