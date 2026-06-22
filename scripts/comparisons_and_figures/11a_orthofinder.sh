#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=chain
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=64

source /gpfs/fs7/grdi/genarcc/wp3/judsonb/miniconda/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/genome_comp/
conda activate ../miniconda/envs/ortho_env/
  
CHIN="../../../common/genomes/Chinook_Salmon/GCF_018296145.1_Otsh_v2.0_genomic"
COHO="../../../common/genomes/Coho_Salmon/GCF_002021735.2_Okis_V2_genomic"

# Keep longest isoform only for each transcript
agat_sp_keep_longest_isoform.pl --gff "$CHIN".gff -o "$CHIN"_longestisoform.gff
agat_sp_keep_longest_isoform.pl --gff "$COHO".gff -o "$COHO"_longestisoform.gff

agat_sp_extract_sequences.pl -g "$CHIN"_longestisoform -f "$CHIN".fna -o chin_isoform_sequences.faa
agat_sp_extract_sequences.pl -g "$COHO"_longestisoform -f "$COHO".fna -o coho_isoform_sequences.faa

# Then run orthofinder on the two proteomes above. t = 64.