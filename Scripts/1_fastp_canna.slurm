#!/bin/bash
#SBATCH --job-name=fastp_canna # Name for job
#SBATCH --time=0-7:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --account=X
#SBATCH --array=1-82
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


INPUT1="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list1.txt | tr -d '\r\n')"
INPUT2="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list2.txt | tr -d '\r\n')"
OUTPUT1="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list3.txt | tr -d '\r\n')"
OUTPUT2="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list4.txt | tr -d '\r\n')"

# Load the necessary environment module first
module load StdEnv/2020

# Load the specific version of fastp
module load fastp/0.23.4

# Check if fastp is available
which fastp

fastp -w  4 -i $INPUT1 -I $INPUT2 -o /home/ahmccorm/scratch/Cannabis/trimmed_fastq/$OUTPUT1\_trimmed.fastq.gz -O /home/ahmccorm/scratch/Cannabis/trimmed_fastq/$OUTPUT2\_trimmed.fastq.gz
