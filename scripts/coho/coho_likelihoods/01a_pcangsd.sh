#!/bin/bash

#SBATCH --job-name=co_pcadapt
#SBATCH --account=grdi_genarcc
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5G
#SBATCH --time=72:00:00
#SBATCH --begin=now+48hours

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset/07c_angsd_n650
conda activate ../../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"

FULLBEAGLE="angsd_n650.beagle.gz"

pcangsd -b "$FULLBEAGLE" -t 64 -o ../08c_pcangsd_n650/angsd_n650_pcadapt --pcadapt --iter 500 --maf 0.05
