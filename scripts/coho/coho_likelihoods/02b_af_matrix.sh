#!/bin/bash

#SBATCH --job-name=GL_afmat
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../../healyt/envs/vcftools

# Define variables.
MAFFOLDER="10b_pop_afs_n650"
POPS=`ls -1 "$MAFFOLDER"/*.mafs.gz | awk -F"/" '{print $2}' | perl -pe 's/\.mafs\.gz//g'`
POPNUM=`ls -1 "$MAFFOLDER"/*.mafs.gz | wc -l`
TMP="temp"
COMMON_SITES="all_pops_common_snps_n650.txt"
OUTPUT="11_rda_afs"
AFS="common_afs"

# Find the SNPs common to all populations by removing SNPs that do not occur as many times as there are populations being analyzed.
sort <(zcat "$MAFFOLDER"/*.mafs.gz | tail -n+2 | awk '{print $1"_"$2}') -T "$TMP" | uniq -c | awk -v npop="$POPNUM" '$1==npop {print $2}' > 01_info_files/"$COMMON_SITES"

# Manipulate the maf files for each population to make an allele frequency matrix.
for file in $(ls -1 "$MAFFOLDER"/*.mafs.gz | perl -pe 's/\.mafs\.gz//g')
do
	name=$(basename $file)

	# Common site file for each population.
	awk 'NR==FNR {a[$1]; next} ($1"_"$2) in a' 01_info_files/"$COMMON_SITES" <(zcat "$MAFFOLDER"/"$name".mafs.gz) > "$OUTPUT"/"$AFS"/"$name".common_mafs.txt

	# Arranges the common maf files into a matrix with a header column in the form of: population name, locus 1, locus 2, locus 3, etc.
	if [[ $name == $(ls -1 "$MAFFOLDER"/*.mafs.gz | awk -F"/" '{print $2}' | perl -pe 's/\.mafs\.gz//g' | head -n 1) ]]
        then
        	cat <(awk '{print $1}' <(awk '{print $1"_"$2, 1-$7}' "$OUTPUT"/"$AFS"/"$name".common_mafs.txt) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk '{print "population\t", $0}') <(awk '{print $2}' <(awk '{print $1"_"$2, 1-$7}' "$OUTPUT"/"$AFS"/"$name".common_mafs.txt) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk -v p=$name '{print p"\t", $0}') > "$OUTPUT"/pop_af_matrix.txt

	else
		cat "$OUTPUT"/pop_af_matrix.txt <(awk '{print $2}' <(awk '{print $1"_"$2, 1-$7}' "$OUTPUT"/"$AFS"/"$name".common_mafs.txt) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk -v p=$name '{print p"\t", $0}') > "$OUTPUT"/temp.txt && mv "$OUTPUT"/temp.txt "$OUTPUT"/pop_af_matrix.txt

	fi

done

# Assemble the above into a matrix.
cut --complement -d$'\t' -f2 "$OUTPUT"/pop_af_matrix.txt > "$OUTPUT"/temp.txt && mv "$OUTPUT"/temp.txt "$OUTPUT"/pop_af_matrix.txt
