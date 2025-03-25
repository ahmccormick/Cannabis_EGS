#!/bin/bash
#SBATCH --job-name=7_Canna_indel_ # Name for job
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --account=X
#SBATCH --array=1-82
#SBATCH --mail-user ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type ALL #events for email. ALL will alert you of job beginning/completion/failure


INPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list1.txt | tr -d '\r\n')"
OUTPUT="$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./list2.txt | tr -d '\r\n')"

module load StdEnv/2020 samtools/1.12
module list

samtools index $INPUT

module --force purge

# Initialize Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
# Activate your Conda environment
conda activate anna_conda

# Define the path to the GenomeAnalysisTK.jar file
GATK_JAR="/home/ahmccorm/miniconda3/envs/anna_conda/opt/gatk-3.8/GenomeAnalysisTK.jar"

# Now use GATK 3.8 commands
java -jar $GATK_JAR -T RealignerTargetCreator -R /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna -I $INPUT -o ./$OUTPUT\.intervals
java -jar $GATK_JAR -T IndelRealigner -R /home/ahmccorm/scratch/Reference_sequences/Cannabis_sativa/GCF_900626175.2_cs10_genomic.fna -I $INPUT -targetIntervals ./$OUTPUT\.intervals --consensusDeterminationModel USE_READS -o ./$OUTPUT\_realigned.bam

conda deactivate
