#!/bin/bash

#SBATCH --job-name=lcwgs_dedup
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-860

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../lcwgs_env

INPUT="04b_deduplicated_bams_2425_2"
OUTPUT="05d_clipped_trimmed_bams_2526"
METRICS="metrics"
JAVA_OPTS="-Xmx80G"
TMPDIR="99_temp"

file=`ls "$INPUT"/*.bam | perl -pe 's/\.bam//g' | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Remove duplicate alignments
name=$(basename "$file")
echo "Deduplicating sample: $name"

name=$(basename "$file")
echo "Clipping sample: $name"

bam clipOverlap \
        --in "$INPUT"/"$name".bam \
        --out "$OUTPUT"/"$name".clip.bam \
        --stats

samtools index "$OUTPUT"/"$name".clip.bam

