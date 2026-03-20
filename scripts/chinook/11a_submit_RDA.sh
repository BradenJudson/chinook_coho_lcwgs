#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --partition=standard
#SBATCH --mem=150GB
#SBATCH --account=grdi_genarcc

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../r_env

awk '{print $1}' 15_pop_afs_imputed/pop_af_matrix.txt | tail -n+2 > 01_info_files/af_rownames.txt

Rscript 00_scripts/20f_rda_PCAcond.R 15c_pop_afs_imputed_2MSNPs/pop_af_matrix.txt 01_info_files/pop_pcas_n106.csv optimal_pcas.txt 01_info_files/chinook_bioclim.csv 01_info_files/bioclim_vars4.txt 14_rda/afsGT_4bioclim_2pca_n819/