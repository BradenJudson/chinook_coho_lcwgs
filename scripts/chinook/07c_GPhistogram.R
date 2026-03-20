args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 2)
  stop("Usage: dbMEM.R <imputedGPs> <originalGPs>")

# Install and load necessary packages ------------------------------------------

for (p in c("sqldf","tidyverse","data.table")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
    suppressMessages(require(p, character.only = T))}
  rm(p)
}

# Chromosome information --------------------------------------------------------

# From: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018296145.1/
chrs <- read_tsv("../01_info_files/otsh_sequence_report.tsv", show_col_types = FALSE) %>%
  filter(!`Chromosome name` %in% c("Un", "MT")) # Remove MT and unassembled contigs.

# Count genotype probabilities per bin ------------------------------------------

# Have to leverage sql for this otherwise we do not have enough memory.
gl_hist <- \(csv) {
  
  # To populate from for-loop.
  chr_list <- list()
  
  # For each chromosome, do the following:
  for (LG in c(unique(chrs$`RefSeq seq accession`))) {
    
    # For parsing with sql.
    chromosome <- toString(sprintf("'%s'", LG))
    
    # Because this takes a while, print progress.
    print(paste("Processing", noquote(gsub("'", "", chromosome))))
    
    df <- fn$read.csv.sql(file = csv, eol = "\n",
                          sql = "select * from file where chr == ($chromosome) ") %>%
      pivot_longer(cols = -c("chr", "pos"), values_to = "GP") %>% select(-c(name))
    
    # Split the dataframe above and count genotype probabilities falling into 5% bins.
    # Convert the tabulated GP data and convert into a dataframe.
    bins <- as.data.frame(cut(df$GP, breaks = seq(0, 1, 1/20)) %>% table) %>%
      `colnames<-`(., c("GP_bin", gsub("'", "", chromosome)))
    
    # Reassign to position in the parent list.
    chr_list[[length(chr_list) + 1]] <- bins
    
    rm(df); gc()
    
  }
  
  # Join list components into a dataframe by common bin values.
  # Sum values across each chromosome for genome-wide frequencies.
  gp_bin_counts <- reduce(chr_list, full_join, by = "GP_bin") %>%
    mutate(total  = rowSums(select(., contains("NC"))),
           rightB = as.numeric(gsub("]", "", gsub(".*,", "", GP_bin)))) %>%
    relocate("rightB", 1)
  
  return(gp_bin_counts)
  
}

imputed  <- gl_hist(csv = args[1])
message("Imputed sample information binned for histogram.")
write.csv(imputed,  "imputed_gp_histogram_data.csv",  row.names = F)

original <- gl_hist(csv = args[2])
message("Nonimputed sample information binned for histogram.")
write.csv(original, "original_gp_histogram_data.csv", row.names = F)

# Isolate chromosome, position and genome-wide counts of genotypes per probability bin.
# Column 37 = sum across all chromosomes. Use ds = "dataset" as factor.
dat <- rbind(imputed[, c(1:2,37)] %>% mutate(ds = "Imputed"),
             original[,c(1:2,37)] %>% mutate(ds = "Original")) %>%
  mutate(ds = factor(ds, levels = c("Original", "Imputed")))

# Visualize and save histogram.
ggplot(data = dat[dat$total > 0,],
       aes(x = rightB, y = total/1e6, fill = ds)) +
  geom_col(position = "dodge",  just = 1,
           colour = "black", size = 1/5) +
  scale_x_continuous(breaks = seq(0, 1, 1/20)) +
  labs(y = "SNPs (millions)", x = "Genotype probability") +
  scale_fill_manual(values = c("grey90", "grey30")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9))
ggsave("gp_hist.tiff", dpi = 300, width = 8, height = 6)

