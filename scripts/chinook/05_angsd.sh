#!/bin/bash

#SBATCH --job-name=ch_angsd
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --time=32:00:00
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-217

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../lcwgs_env

# Define some variables.
REGION=`cat ../chinook_10mb_genome.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon/index_faidx
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
BAMFILES="bam_list_n676_downSampled.txt"
OUTFOLDER="07_angsd_n676_maf005"

#ls -1 06_downsampled_bams/*.bam > 01_info_files/"$BAMFILES"

# Run ANGSD.
angsd -b ./01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/chinook_angsd_"$REGION" \
        -nThreads 8 \
        -r "$REGION" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -SNP_pval 1e-10  \
        -minMapQ 30 \
        -minMaf 0.05 \
        -minQ 30 \
        -setMinDepth 676 \
        -setMaxDepth 2704 \
        -minInd 575 \
        -doMaf 1 \
        -doMajorMinor 1 \
        -doCounts 1 \
        -GL 1 -doGLF 2 \
        -only_proper_pairs 1 \
        -dumpCounts 2 \
        -rmTriallelic 1e-6 \
        -doBcf 1 -doPost 1 \
        -doGeno 5 -doQsDist 1