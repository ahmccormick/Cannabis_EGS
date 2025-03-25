#!/bin/bash
#SBATCH --job-name=10_LR_filter_canna # Name for job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00-6:00:00
#SBATCH --account=X
#SBATCH --mail-user=ahmccorm@hawaii.edu # email wish to be notified at
#SBATCH --mail-type=ALL #events for email. ALL will alert you of job beginning/completion/failure

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load htslib/1.11
module load vcftools/0.1.16

# Use local installation of bcftools 1.11
export PATH=/home/ahmccorm/projects/def-yeaman/ahmccorm/bcftools_1.11/bcftools-1.11/bin:$PATH

# Verify input files
echo "Listing input VCF files:"
ls -lh ./*vcf.gz

# Step 1: Filter with bcftools
echo "Running bcftools filter..."
bcftools filter -e 'MQ < 30' ./*vcf.gz -Oz -o ./tmp.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error in bcftools filter step."
    exit 1
fi

# Verify intermediate file
echo "Verifying intermediate VCF file:"
ls -lh ./tmp.vcf.gz

# Step 2: Filter with vcftools
echo "Running vcftools..."
vcftools --gzvcf ./tmp.vcf.gz --max-missing 0.7 --minQ 30 --minGQ 20 --minDP 5 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > ./Cannabis_sativa_PRJNA734114_filtered.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error in vcftools filtering step."
    exit 1
fi

# Step 3: Compress and index the final VCF file
echo "Compressing and indexing the final VCF file..."
tabix -p vcf ./Cannabis_sativa_PRJNA734114_filtered.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error in tabix indexing step."
    exit 1
fi

# Clean up intermediate file
echo "Cleaning up intermediate file..."
rm ./tmp.vcf.gz

echo "Script completed successfully."
