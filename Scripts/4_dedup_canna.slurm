#!/bin/bash
#SBATCH --job-name=Canna_dedup # Name for job
#SBATCH --time=0-6:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --account=X
#SBATCH --array=1-82
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


INPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list1.txt | tr -d '\r\n')"
OUTPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list2.txt | tr -d '\r\n')"


module load picard java


java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$INPUT OUTPUT=/home/ahmccorm/scratch/Cannabis/post_bwa/bwa_output/de_dup/$OUTPUT\_dedup.bam METRICS_FILE=/home/ahmccorm/scratch/Cannabis/post_bwa/bwa_output/de_dup/$OUTPUT\_DUP_metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
