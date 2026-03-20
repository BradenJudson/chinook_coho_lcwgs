#!/bin/bash


#SBATCH --job-name=lcwgs_bwa
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-195

# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset
conda activate ../lcwgs_env

GENOME_FOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Chinook_Salmon
GENOME="GCF_018296145.1_Otsh_v2.0_genomic.fna"
TRIMMED_FOLDER="02_trimmed_reads"
ALIGNED_FOLDER="03_alignments"
NCPU=64
TEMP_FOLDER="99_temp"

# Specify individual files for alignment.
file=`ls -1 ../../savager/chinook_lcwgs/raw_reads/*_R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | awk -F '/' '{print $NF}' | perl -pe 's/_R1.fastq.gz//'`
	
	# For tidier writing/naming downstream, only use the sample ID not the full fastq info.
        base=$(basename $file | perl -pe 's/_R1.*//g' | awk '{gsub(".*IDT_i5_[0-9]{1,3}.", "")}1')
        echo "Aligning $base"

	# BAM read group info required for downstream analyses. 
        ID=$(echo $base | perl -pe 's/\.[a-zA-Z0-9]{2,4}[_-][A-Z0-9].*$//g')
        PU=$(echo $base | perl -pe 's/NS\.X//g' | perl -pe 's/\.IDT_i7.*$//g')
        LB=$(echo $base | perl -pe 's/^.*[0-9]*\.IDT/IDT/g')
        SM=$(echo $base | perl -pe 's/NS.*i5.*_[0-9]*\.//g')

	# Conduct the alignment and print flagstat file for each sample.
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
