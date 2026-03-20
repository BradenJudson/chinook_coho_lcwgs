#!/bin/bash

#SBATCH --job-name=lcwgs_pca
#SBATCH --account=grdi_genarcc
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5G
#SBATCH --time=72:00:00

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset/07d_angsd_n819_maf005
conda activate ../../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"

FULLBEAGLE="angsd_n819.beagle.gz"

pcangsd -b "$FULLBEAGLE" -t 64 -o ../08_pcangsd/angsd_n819_pcadapt --pcadapt --iter 500 --maf 0.05