#!/bin/bash


#SBATCH --job-name=lcwgs_fastp
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --account=grdi_genarcc

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../lcwgs_env

LENGTH=100
QUAL=20
INPUT="raw_reads"
OUTPUT="02_trimmed_reads"
METRICS="metrics"
NCPU=16

# Load GNU parallel.
parallel=/fs/ssm/main/opt/parallel/parallel-20210922/parallel_20210922_ubuntu-20.04-amd64-64/bin/parallel

# Export for parallelization.
export fastp

ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/R[12]\.fastq\.gz//g' | \
        parallel -j "$NCPU" \
                fastp -i {}R1.fastq.gz -I {}R2.fastq.gz \
                -o "$OUTPUT"/{/}R1.fastq.gz \
                -O "$OUTPUT"/{/}R2.fastq.gz \
                --length_required="$LENGTH" \
                --qualified_quality_phred="$QUAL" \
                --correction \
                --trim_tail1=1 \
                --trim_tail2=1 \
                --json 11_metrics/{/}.json \
                --html 11_metrics/{/}.html \
                --report_title={/}.html
