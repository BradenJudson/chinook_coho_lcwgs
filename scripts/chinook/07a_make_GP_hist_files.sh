#!/bin/bash

#SBATCH --job-name=geno_probs
#SBATCH --account=grdi_genarcc
#SBATCH --time=16:00:00
#SBATCH --mem=80GB

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset/17_gps
conda activate ../../lcwgs_env

Imp_VCF="../09_vcf/chinook_lcwgs_maf005_n676_imputed_concat.vcf.gz"

echo "Last array. Generating full maximum likelihood file."
# Combine outputs and write to csv.
paste genomic_positions.txt <(paste ./temp_array_imp/*.txt) | awk '{print NR, $0}' | sed 's/[[:space:]]\+/,/g' > maximum_imputed_GLs.csv
paste genomic_positions.txt <(paste ./temp_array_ori/*.txt) | awk '{print NR, $0}' | sed 's/[[:space:]]\+/,/g' > maximum_original_GLs.csv

# Define header: Row number, chromosome name, position, and individual ID. One for each VCF.
GP_HEADER=`echo -e "row\nchr\npos\n$(bcftools query -l "$Imp_VCF")" | tr '\n' ',' | sed '$ s/.$//'`

# Add header to existing csvs.
sed -i "1i $GP_HEADER" maximum_imputed_GLs.csv
sed -i "1i $GP_HEADER" maximum_original_GLs.csv

echo "Maximum likelihood files created. Running R script for creating the histogram."
conda activate ../../r_env
export TMPDIR=tmp
Rscript 09d_GPhistogram.R maximum_imputed_GLs.csv maximum_original_GLs.csv