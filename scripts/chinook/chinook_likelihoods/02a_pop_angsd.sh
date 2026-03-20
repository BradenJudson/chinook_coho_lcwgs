#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=chin_pop
#SBATCH --time=12:00:00
#SBATCH --mem=100GB
#SBATCH --array=1-106

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../lcwgs_env

# Define variables.
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
OUTFOLDER="10c_pop_angsd_n819"
SITES="01_info_files/sites_n819.txt"
TEXT="10c_pop_angsd_n819/txt_files"
SAMPLE_INFO="01_info_files/chinook_lcwgs_samples.csv"
CHROMS="01_info_files/chinook_chroms.txt"
BEAGLE="18_ld/pruned_snps_n819.beagle.gz"
POP=`awk -F, '{if($16 != "Drop population") print $6}' "$SAMPLE_INFO" | sort | uniq | head -n -1 | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] && [[ ! -f "$SITES" ]];
then

        # Define sites to use with predefined and consistent major/minor alleles for each locus.
        zcat "$BEAGLE" | awk '{OFS="\t"} NR>1 {split($1, a, "_"); print a[1]"_"a[2], a[3], $2, $3}' - > "$SITES"
        angsd sites index "$SITES"
fi

# Return site frequencies from genotype likelihoods for each population.
# Create a variable for each population that does not have spaces - makes file names nicer to work with (e.g., San Juan -> SanJuan).
POPNS=$(echo $POP | tr -d ' ')

# Creates a list of individuals to analyze for each population. As written, this filters out individuals that we do not want to consider.
awk -F, -v population="$POP" '{ if ($16 == "" && $6 == population) print $3}' "$SAMPLE_INFO" > "$TEXT"/"$POPNS"_bams.txt

# Run ANGSD. Maintain mapping filters, but site-specific filters rely on the provided sites file.
# Returns minor allele frequencies and SFSs for all SNPs present in at least five individuals.
angsd -bam "$TEXT"/"$POPNS"_bams.txt \
        -ref "$GENOMEFOLDER"/"$GENOME" \
        -anc "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/"$POPNS" \
        -rf  "$CHROMS" \
        -sites "$SITES" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -minMapQ 30 \
        -minQ 30 \
        -only_proper_pairs 1 \
        -GL 1 -doMaf 1 \
        -doMajorminor 3 \
        -doSaf 1 \
        -minInd 5 \
	-doHWE 1 -minHetFreq 0 -maxHetFreq 1