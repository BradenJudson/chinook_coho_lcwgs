library(tidyverse); library(cowplot)


data <- readxl::read_excel("data/coho_chinook_samples.xlsx") %>% 
  filter(!fate %in% c("Drop population", "Drop (PCA outlier)"))
data$final_reads <- as.numeric(data$final_reads)
data$final_coverage <- as.numeric(data$final_coverage)
data$species <- tools::toTitleCase(data$species)
data$seq_batch <- fct_infreq(as.factor(data$seq_batch))

data <- data %>% mutate(fate = case_when(
  downsample_fraction  < 1 & is.na(fate) ~ "Downsample",
  downsample_fraction == 1 ~ "Keep as-is",
  TRUE ~ fate
))

readcov <- \(sp) {
  
  ggplot(data = data[data$species == sp,],
         aes(y = final_coverage,
             x = final_reads/1e6)) +
    scale_shape_manual(values = c(25, 23, 21)) +
    geom_segment(aes(x = mean(data$final_reads, na.rm = T)/1e6,
                     xend = mean(data$final_reads, na.rm = T)/1e6,
                     y = -Inf, 
                     yend = Inf),
                 inherit.aes = FALSE,
                 linewidth = 1/2, 
                 linetype = 2) +
    geom_segment(aes(x = Inf,
                     xend = -Inf,
                     y = mean(data$final_coverage, na.rm = T),
                     yend = mean(data$final_coverage, na.rm = T)),
                 linewidth = 1/2, linetype = 2,
                 inherit.aes = FALSE) +
    geom_point(colour = "black",
               size = 2,
               alpha = 4/5,
               aes(fill = seq_batch,
                   shape = fate)) +
    labs(x = "Reads (millions)",
         y = "Average Individual Coverage") +
    guides(fill = guide_legend(override.aes = list(shape = 21),
                               nrow = 2), 
           shape = guide_legend(position = "right")) +
    theme_bw() +
    theme(panel.grid = element_line(color = "gray95"),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.background = element_rect(fill = "transparent", color = NA))
  
}

(co <- readcov("Coho") + guides(shape = 'none'))
(ch <- readcov("Chinook"))

cowplot::plot_grid(plotlist = list(co, ch), 
                   ncol = 2, rel_widths = c(0.85, 1),
                   labels = c("Coho", "Chinook"),
                   label_x = 0.03,
                   label_y = 0.85)

ggsave("plots/reads_cov.tiff", dpi = 300, width = 15,
       height = 7)

