#!/bin/bash
#SBATCH --job-name=bwa_Canna # Name for job
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --account=def-yeaman
#SBATCH --array=1-82
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure

INPUT1="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list1.txt | tr -d '\r\n')"
INPUT2="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list2.txt | tr -d '\r\n')"
OUTPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list3.txt | tr -d '\r\n')"

module load bwa samtools

bwa mem -t 4 /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna "$INPUT1" "$INPUT2" > ./bwa_output/$OUTPUT\.sam
samtools view -Sb -q 10 ./bwa_output/$OUTPUT\.sam > ./bwa_output/$OUTPUT\.bam
rm ./bwa_output/$OUTPUT\.sam
samtools sort --threads  4 ./bwa_output/$OUTPUT\.bam > ./bwa_output/$OUTPUT\_sorted.bam
rm ./bwa_output/$OUTPUT\.bam
#samtools index ./bwa_output/$OUTPUT\_sorted.bam
# Create a CSI index for the sorted BAM file
samtools index -c ./bwa_output/$OUTPUT\_sorted.bam
