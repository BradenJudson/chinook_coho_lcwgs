#!/bin/bash

#SBATCH --job-name=vcffrqs
#SBATCH --account=grdi_genarcc
#SBATCH --mem=70G
#SBATCH --time=30:00:00

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env


INPUT="09_vcf/coho_lcwgs_maf005_n650_imputed.vcf.gz"
OUTPUT="10a_pop_afs_n650_imputed"
SAMPLEINFO="01_info_files/pop_map_n650.txt"

for pop in $(awk '{print $3}' "$SAMPLEINFO" | sort | uniq)
do
	echo $pop

	#awk -v population="$pop" '{ if ($3 == population) {print $1}}' "$SAMPLEINFO" > "$OUTPUT"/txt_files/"$pop"_indvs.txt

	#vcftools --gzvcf "$INPUT" --keep "$OUTPUT"/txt_files/"$pop"_indvs.txt --freq --out "$OUTPUT"/"$pop"_imputed_freqs


	if [[ $pop == $(awk '{print $3}' "$SAMPLEINFO" | sort | uniq | head -n 1) ]]
	then
		cat <(awk '{print $1}' <(awk '{ split ($5,a,":"); if (NR!=1) print $1"_"$2, a[2]}' "$OUTPUT"/"$pop"_imputed_freqs.frq) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk '{print "population\t", $0}') <(awk '{print $2}' <(awk '{ split ($5,a,":"); if (NR!=1) print $1"_"$2, a[2]}' "$OUTPUT"/"$pop"_imputed_freqs.frq) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk -v p=$pop '{print p"\t", $0}') > "$OUTPUT"/pop_af_matrix_imputed.txt
	else
		cat "$OUTPUT"/pop_af_matrix_imputed.txt <(awk '{print $2}' <(awk '{ split ($5,a,":"); if (NR!=1) print $1"_"$2, a[2]}' "$OUTPUT"/"$pop"_imputed_freqs.frq) | tr '\n' '\t' | perl -pe 's/\t$/\n/g' | awk -v p=$pop '{print p"\t", $0}') > "$OUTPUT"/temp.txt && mv "$OUTPUT"/temp.txt "$OUTPUT"/pop_af_matrix_imputed.txt

	fi

done

