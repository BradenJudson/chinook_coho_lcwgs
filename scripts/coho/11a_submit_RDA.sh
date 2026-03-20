#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --partition=standard
#SBATCH --mem=150GB
#SBATCH --account=grdi_genarcc

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../r_env

#awk '{print $1}' 11_rda_afs/pop_afgls_matrix.txt | tail -n+2 > 01_info_files/af_rownames.txt

Rscript 00_scripts/13b_runRDA_imputed.R 10a_pop_afs_n650_imputed/pop_af_matrix_imputed.txt 01_info_files/pop_pcas_afsGL_n650.csv 01_info_files/optimal_pcas.txt 01_info_files/coho_bioclim_n83pops.csv 01_info_files/bioclim_vars4.txt 12_rdas/afsGT_4bioclim_2PCA