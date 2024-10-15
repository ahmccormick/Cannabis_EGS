#!/bin/bash
#SBATCH --job-name=RG_Canna # Name for job
#SBATCH --time=0-6:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --account=def-yeaman
#SBATCH --array=1-82
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure

INPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list1.txt | tr -d '\r\n')"
OUTPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list2.txt | tr -d '\r\n')"
NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list3.txt | tr -d '\r\n')"

module load picard java

#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$INPUT O=$OUTPUT\_RG.bam RGID=$NAME RGLB=$NAME\_LB RGPL=ILLUMINA RGPU=unit1 RGSM=$NAME

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$INPUT O=/home/ahmccorm/scratch/Cannabis/post_bwa/bwa_output/de_dup/add_RG/${OUTPUT}\_RG.bam RGID=$NAME RGLB=$NAME\_LB RGPL=ILLUMINA RGPU=unit1 RGSM=$NAME
