#!/bin/bash


#SBATCH --job-name=coho_impute
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-30

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset/07c_angsd_n650

REGION=`cat ../01_info_files/coho_chroms.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
CHROM=`cat ../01_info_files/coho_chroms.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | awk -F":" '{print $1}'`
echo $CHROM

VCF="coho_lcwgs_maf005_n650"
FOLDER="../09_vcf"
JAVA_OPTS="-Xmx80G"
NCPU=64

conda activate ../../lcwgs_env

# Ouputs compressed vcf file from 10mb subset bcfs.
#bcftools concat *.bcf -O z -o "$FOLDER"/"$VCF".vcf.gz --threads $NCPU

# Have to switch conda environments because beagle v4 and bcftools are tricky to get to work together.
conda activate ../../../healyt/envs/beagle4

beagle "$JAVA_OPTS" \
        gl="$FOLDER"/"$VCF".vcf.gz \
        out="$FOLDER"/"$VCF"_"$CHROM"_imputed \
        gprobs=true \
        chrom="$REGION" \
        nthreads=$NCPU

conda activate ../../lcwgs_env

bcftools index "$FOLDER"/"$VCF"_"$CHROM"_imputed.vcf.gz