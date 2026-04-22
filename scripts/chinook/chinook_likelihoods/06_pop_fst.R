
library(tidyverse); library(pheatmap)

# Read in Fst data from previous script ----------------------------------------

fst <- read.delim("data/fst_data/chinook_likelihoods_n819_fst.txt", sep = "") %>% 
  mutate(pops = gsub("Kitimat_", "Kitimat", gsub(".globalFst", "", 
                gsub("12c_pop_fst_n819/global_fst/", "", file)))) %>%
  mutate(pop1 = gsub(".*_", "", pops), pop2 = gsub("_.*", "", pops)) %>% 
  filter(pop1 != "KitimatOT" & pop2 != "KitimatOT") %>% 
  dplyr::select(c("pop1", "pop2", "FstUnweighted", "FstWeighted"))

# Convert distances into a pairwise matrix with a diagonal of zero.
matrix_names <- sort(unique(as.character(unlist(fst[c("pop1", "pop2")]))))
fst_mat <- matrix(0, length(matrix_names), length(matrix_names),
              dimnames = list(matrix_names, matrix_names))
fst_mat[as.matrix(fst[c("pop1", "pop2")])] <- fst$FstWeighted
fst_mat[as.matrix(fst[c("pop2", "pop1")])] <- fst$FstWeighted
diag(fst_mat) <- NA
class(fst_mat) <- "numeric"

# Visualize heatmap with pheatmap.
png(width = 1000, height = 1000, units = "px", filename = "plots/popfst.png")
pheatmap(fst_mat, treeheight_row = 0, treeheight_col = 0, na_col = NA)
dev.off()
