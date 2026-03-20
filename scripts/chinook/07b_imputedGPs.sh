#!/bin/bash

#SBATCH --job-name=geno_probs
#SBATCH --account=grdi_genarcc
#SBATCH --time=2:00:00
#SBATCH --mem=60GB
#SBATCH --array=1-676

# Total number of samples.
LAST=676 

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_offset/17_gps
conda activate ../../lcwgs_env

# Pre- and post-imputation VCFs.
Imp_VCF="../09_vcf/chinook_lcwgs_maf005_n676_imputed_concat.vcf.gz"
Ori_VCF="../09_vcf/chinook_lcwgs_maf005_n676.vcf.gz"
SNPs=`bcftools query -f "%POS\n" "$Imp_VCF" | wc -l`

# On the first sample we need to create the source files.
# I.e., the genotype probabilities for every sample and SNP and their positions.
# Only generate files if they do not already exist.
if [ ! -f "imputed_GLs.txt"  ]; then
        bcftools query -f "[\t%GP]\n" "$Imp_VCF" | sed "s/^[\t]*//" > imputed_GLs.txt
fi

if [ ! -f "original_GLs.txt" ]; then
        bcftools query -f "[\t%GP]\n" "$Ori_VCF" | sed "s/^[\t]*//" > original_GLs.txt
fi

# Positions for both VCFs are identical.
if [ ! -f "genomic_positions.txt" ]; then
        bcftools query -f "%CHROM\t%POS\n" "$Imp_VCF" > genomic_positions.txt
fi

# If the imputed and original likelihood files have information for every locus, then isolate
# the maximum likelihood per genotype per individual for both datasets. If the full likelihood
# files are not complete (i.e., not yet information for every SNP) then wait until those files
# are finalized and try again.
while true; do

       original_lines=$(wc -l < original_GLs.txt 2>/dev/null || echo 0)
       imputed_lines=$(wc  -l < imputed_GLs.txt  2>/dev/null || echo 0)

       if (( original_lines == $SNPs && imputed_lines == $SNPs )); then

        # Identify which individual to calculate genotype probabilities with.
       sample=`bcftools query -l $Imp_VCF | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | perl -pe 's/.dedup.clip.bam//g' | perl -pe 's/.dedup.clip.downsampled.bam//g' | awk -F'/' '{print $NF}'`
       echo $sample

       # Isolate chromosomes and genomic positions per SNP to add later on.
       cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' imputed_GLs.txt | tr "," "\t" | perl -lane '@a=sort @F;print join "\t", @a' | awk '{ print $3 }') > ./temp_array_imp/$sample.txt

       cat <(awk -F "\t" '{ print $'$SLURM_ARRAY_TASK_ID' }' original_GLs.txt | tr "," "\t" | awk '{for (i=1;i<=NF;i++) if ($i+0 == $i && $i ~ /e/) $i = sprintf("%.10f", $i) } 1' |\
                awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' | awk '{ print $3 }') | awk '{printf "%.2f\n",$1}' > ./temp_array_ori/$samp>#
               break
       else

       echo "Waiting for likelihood files to complete!"
       sleep 120

       fi
done

echo "Sample-specific maximum likelihood files written."

# When information is written for the last sample, we can combine information across
# all samples and use R to make our histogram of genotype probabilities.
# This last array takes significantly longer than the others, so we run it as a separate script.
if [[ $SLURM_ARRAY_TASK_ID == $LAST ]]
then

  sbatch 09c_make_GP_hist_files.sh 

fi
