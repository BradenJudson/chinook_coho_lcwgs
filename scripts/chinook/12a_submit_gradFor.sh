#!/bin/bash

#SBATCH --time=15:00:00
#SBATCH --partition=standard
#SBATCH --mem=150GB
#SBATCH --account=grdi_genarcc
#SBATCH --job-name=lc_gf

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../r_env

Rscript 00_scripts/18_gForest.R 14_rda/afsGT_4bioclim_2pca_n819/outlier_afs_matrix_afsGL4bio2PCA_n819.txt 01_info_files/pcnms_n106pops.csv 01_info_files/chinook_bioclim.csv 01_info_files/bioclim_vars4.txt