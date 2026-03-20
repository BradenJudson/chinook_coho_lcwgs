#!/bin/bash

#SBATCH --time=15:00:00
#SBATCH --partition=standard
#SBATCH --mem=150GB
#SBATCH --account=grdi_genarcc
#SBATCH --job-name=lc_gf

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../r_env

Rscript 00_scripts/14b_gForest.R 12_rdas/afsGT_4bioclim_2PCA/outlier_afs_matrix_afsGT4bio2PCA_n650_2RDAs_3Msnps.txt 01_info_files/coho_bioclim_n83pops.csv 01_info_files/bioclim_vars4.txt
