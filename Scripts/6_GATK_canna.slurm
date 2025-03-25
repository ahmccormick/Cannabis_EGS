#!/bin/bash
#SBATCH --job-name=6_gatk_Canna # Name for job
#SBATCH --time=0-0:50
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --account=X
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


module load StdEnv/2020 samtools/1.12

samtools faidx /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna

module purge
module load picard java

#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=GCF_000331145.1_ASM33114v1_genomic.fna O=GCF_000331145.1_ASM33114v1_genomic.dict
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=GCF_900626175.2_cs10_genomic.fna O=GCF_900626175.2_cs10_genomic.dict
