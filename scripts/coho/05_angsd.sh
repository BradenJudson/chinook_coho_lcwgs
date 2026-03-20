#!/bin/bash

#SBATCH --job-name=co_angsd
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --time=32:00:00
#SBATCH --mem=150GB
#SBATCH --array=1-30


source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env

REGION=`cat ../coho_chroms.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
BAMFILES="bam_list_n650.txt"
OUTFOLDER="07c_angsd_n650"
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Coho_Salmon/
GENOME="GCF_002021735.2_Okis_V2_genomic.fna"


angsd -b ./01_info_files/"$BAMFILES"  \
        -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/coho_angsd_"$REGION" \
        -nThreads 8 \
        -r "$REGION" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -SNP_pval 1e-10  \
        -minMapQ 30 \
        -minMaf 0.05 \
        -minQ 30 \
        -setMinDepth 650 \
        -setMaxDepth 2600 \
        -minInd 553 \
        -doMaf 1 \
        -doMajorMinor 1 \
        -doCounts 1 \
        -GL 1 -doGLF 2 \
        -only_proper_pairs 1 \
        -dumpCounts 2 \
        -rmTriallelic 1e-6 \
        -doBcf 1 -doPost 1 \
        -doGeno 5 -doQsDist 1