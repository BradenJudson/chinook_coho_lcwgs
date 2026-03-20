#!/bin/bash

#SBATCH --partition=standard
#SBATCH --account=grdi_genarcc
#SBATCH --cpus-per-task=64
#SBATCH --time=06:00:00
#SBATCH --array=1-861

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env

file="01_info_files/coho_lcwgs_samples.csv"
# Predefined file with the bams and the fraction they need to be downsampled to to achieve the goal coverage.
BAM=`awk -F, '{ if( $10 != "Discard" || $14 != "PCA outlier" ) print $2 }' "$file"  | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FRAC=`awk -F, '{ if( $10 != "Discard" || $14 != "PCA outlier" ) print $12 }' "$file" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

OUT="06_downsampled_bams"
name=`echo "$BAM" |  awk -F '/' '{print $NF}' | perl -pe 's/\.bam//g' |  awk '{gsub(".*IDT_i5_[0-9]{1,3}.", "")}1'`

if (( $(echo "$FRAC < 1" | bc -l) )); then

        # Perform random downsampling. -s ensures that a read is never kept without its mate. Index bam.
        samtools view -b --subsample-seed 240 -s $FRAC ./"$BAM" > "$OUT"/"$name".downsampled.bam
        samtools index "$OUT"/"$name".downsampled.bam

else

        # Create symlink of unchanged original bam.
        ln -s "$BAM" "$OUT"/"$name".original.bam

fi

# Make a file list of original and downsampled files.
cat <(find "$OUT"/ -type f) <(readlink "$OUT"/*.bam) > 01_info_files/bam_list_n660sym_downSampled.txt
