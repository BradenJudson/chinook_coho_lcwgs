#!/bin/bash


#SBATCH --job-name=lcwgs_bwa
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=24
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-660

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env


GENOME_FOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Coho_Salmon
GENOME="GCF_002021735.2_Okis_V2_genomic.fna"
TRIMMED_FOLDER="03_trimmed_reads"
ALIGNED_FOLDER="04_alignments"
NCPU=64
TEMP_FOLDER="temp"

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/coho_offset
conda activate ../lcwgs_env


GENOME_FOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Coho_Salmon
GENOME="GCF_002021735.2_Okis_V2_genomic.fna"
TRIMMED_FOLDER="03b_trimmed_reads_skeena"
ALIGNED_FOLDER="04_alignments"
NCPU=64
TEMP_FOLDER="temp"

file=`ls -1 "$TRIMMED_FOLDER"/*_R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | awk -F '/' '{print $NF}' | perl -pe 's/_R1.fastq.gz//'`

        base=$(basename $file | perl -pe 's/_R1.*//g' | awk '{gsub(".*IDT_i5_[0-9]{1,3}.", "")}1')
        echo "Aligning $base"

        ID=$(echo $base | perl -pe 's/\.[a-zA-Z0-9]{2,4}[_-][A-Z0-9].*$//g')
        PU=$(echo $base | perl -pe 's/NS\.X//g' | perl -pe 's/\.IDT_i7.*$//g')
        LB=$(echo $base | perl -pe 's/^.*[0-9]*\.IDT/IDT/g')
        SM=$(echo $base | perl -pe 's/NS.*i5.*_[0-9]*\.//g')

        bwa mem -t "$NCPU" \
                "$GENOME_FOLDER"/"$GENOME" \
                "$TRIMMED_FOLDER"/"$file"_R1.fastq.gz \
                "$TRIMMED_FOLDER"/"$file"_R2.fastq.gz \
                -R "@RG\tID:$ID\tPL:illumina\tPU:$PU\tSM:$SM\tLB:$LB" |
                tee >(samtools flagstat - > "$ALIGNED_FOLDER"/"$base".flagstat) |
                samtools view -b -q 20 - |
                samtools sort -T "$TEMP_FOLDER"/"$base" - > "$ALIGNED_FOLDER"/"$base".bam

        samtools index "$ALIGNED_FOLDER"/"$base".bam

        samtools idxstat "$ALIGNED_FOLDER"/"$base".bam > "$ALIGNED_FOLDER"/"$base".idxstat