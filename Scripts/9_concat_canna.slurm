#!/bin/bash
#SBATCH --job-name=9_canna_ccat # Name for job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=00-12:00:00
#SBATCH --account=X
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load htslib/1.11

export PATH=/home/ahmccorm/projects/def-yeaman/ahmccorm/bcftools_1.11/bcftools-1.11/bin:$PATH


tr -d '\r' < ./list.txt > ./cleaned_list.txt

bcftools concat -f ./cleaned_list.txt -Oz > ./concat/Cannabis_sativa_PRJNA734114.vcf.gz      #think list.txt is the chromosome VCFs in a list as thats the last output?
tabix -p vcf ./concat/Cannabis_sativa_PRJNA734114.vcf.gz
