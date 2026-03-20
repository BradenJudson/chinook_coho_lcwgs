#!/bin/bash

#SBATCH --job-name=offCo_clip
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-660
#SBATCH --dependency=afterok:5554224

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env

INPUT="05_deduplicated_bams"
OUTPUT="06_clipped_trimmed_bams"
METRICS="metrics"
JAVA_OPTS="-Xmx80G"
TMPDIR="temp"

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

