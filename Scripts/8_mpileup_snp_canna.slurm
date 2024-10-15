#!/bin/bash
#SBATCH --job-name=8_Canna_mSNP # Name for job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00-24:00:00
#SBATCH --account=def-yeaman
#SBATCH --array=1-10
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure

module load StdEnv/2020
module load gcc/9.3.0


# Use local installation of bcftools 1.11
export PATH=/home/ahmccorm/projects/def-yeaman/ahmccorm/bcftools_1.11/bcftools-1.11/bin:$PATH


CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./chromosomes.txt | tr -d '\r\n')"

#tr -d '\r' < ./list.txt > ./cleaned_list.txt
#tr -d '\r' < ./ploidymap.txt > ./cleaned_ploidy.txt   #NEED TO KNOW LAYOUT HERE

bcftools mpileup -Ou -f /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna --bam-list ./list.txt -q 5 -r $CHROM -I -a FMT/AD | bcftools call -S ./ploidymap.txt -G - -f GQ -mv -Ov > ./vcf_output/$CHROM\.vcf
