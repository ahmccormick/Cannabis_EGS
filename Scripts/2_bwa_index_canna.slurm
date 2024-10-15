#!/bin/bash
#SBATCH --job-name=Canna_index # Name for job
#SBATCH --time=0-6:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --account=def-yeaman
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


module load bwa

bwa index -a bwtsw /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna
