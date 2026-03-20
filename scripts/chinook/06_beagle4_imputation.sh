#!/bin/bash


#SBATCH --job-name=otsh_impute
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-34

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset/07_angsd_n676_maf005

REGION=`cat ../01_info_files/chinook_chroms.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
VCF="chinook_lcwgs_maf005_n676"
FOLDER="../09_vcf"
JAVA_OPTS="-Xmx80G"
NCPU=32
LAST=34

conda activate ../../lcwgs_env

# Ouputs compressed vcf file from 10mb subset bcfs.
bcftools concat *.bcf -O z -o "$FOLDER"/"$VCF".vcf.gz --threads $NCPU

# Have to switch conda environments because beagle v4 and bcftools are tricky to get to work together.
conda activate ../../../healyt/envs/beagle4

# Impute based on genotype probabilities one chromosome at a time.
beagle "$JAVA_OPTS" \
        gl="$FOLDER"/"$VCF".vcf.gz \
        out="$FOLDER"/"$VCF"_"$REGION"_imputed \
        gprobs=true \
        chrom="$REGION" \
        nthreads=$NCPU
    
# Switch back to use bcftools again.        
conda activate ../../lcwgs_env
        
bcftools index "$FOLDER"/"$VCF"_imputed.vcf.gz

# Then use bcftools concat to join each chromosome into a single VCF
if [[ $SLURM_ARRAY_TASK_ID == $LAST ]]
then

  bcftools concat -Oz -o "$FOLDER"/"$VCF"_imputed.vcf.gz "$FOLDER"/"$VCF"_*_imputed.vcf.gz 

fi