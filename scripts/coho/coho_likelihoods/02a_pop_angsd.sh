#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --job-name=afs_offco
#SBATCH --time=16:00:00
#SBATCH --mem=60GB
#SBATCH --array=1-83
#SBATCH --begin=now+48hours

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env

# Define variables and paths.
SITES="01_info_files/angsd_sites_gls_n650.txt"
BEAGLE="07c_angsd_n650/angsd_n650.beagle.gz"
DIR="10b_pop_afs_n650_hwe"
TEXT="txt_files"
CHROMS="01_info_files/coho_chroms.txt"
SAMPLE_INFO="01_info_files/pop_map_n650.txt"
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Coho_Salmon
GENOME="GCF_002021735.2_Okis_V2_genomic.fna"
POP=`awk '{print $3}' $SAMPLE_INFO | sort | uniq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

echo "$POP"

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] && [[ ! -f "$SITES" ]];
then

        # Define sites to use with predefined and consistent major/minor alleles for each locus.
        zcat "$BEAGLE" | awk '{OFS="\t"} NR>1 {split($1, a, "_"); print a[1]"_"a[2], a[3], $2, $3}' - > "$SITES"
        angsd sites index "$SITES"
fi

awk -v population="$POP" '{ if ($3 == population) {print $1}}' "$SAMPLE_INFO" > "$DIR"/"$TEXT"/"$POP"_bams.txt

angsd -b "$DIR"/"$TEXT"/"$POP"_bams.txt \
        -ref "$GENOMEFOLDER"/"$GENOME" \
        -anc "$GENOMEFOLDER"/"$GENOME" \
        -sites "$SITES" \
        -rf "$CHROMS" \
        -out "$DIR"/"$POP" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -minMapQ 30 \
        -minQ 30 \
        -only_proper_pairs 1 \
        -GL 1 -doMaf 1 \
        -doSaf 1 \
        -doMajorminor 3 \
        -minInd 5 \
        -doHWE 1 -minHetFreq 0 -maxHetFreq 1