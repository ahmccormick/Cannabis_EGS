# Cannabis Environmental Genomic Selection

This repository provides the supporting data for a recently submitted manuscript entitled "X"

A pre-print of this work can be found at (to add shortly)

High resolution verions of main and supplemental figures are to be available on fig share once accepted https://figshare.com/authors/Anna_H_McCormick/17741367

# VCF Pipeline
Scripts used for .vcf generation are located in the scripts/vcf folder

# EGS Pipeline Steps
Scripts used for EGS are located in the scripts/EGS folder

Inputs needed to start
.vcf file with your genotype information and lat/long data for your samples in a .csv

1. Identification of core (1_step1_core.py)
   This script was used on the filtered vcf file in PyCharm to identify a core (in this case n=25, change accordingly for 30% of your dataset)
   
2. Convert vcf to format for rrBLUP (2_vcf_to_rrBLUP.py)
   In terminal locally with the below command structure:
   python vcf_to_rrBLUP.py -i Cannabis_sativa_PRJNA734114_filtered.vcf.gz -o Ren_PRJNA734114_rrBLUP_format

3. Extract relevant climate data from WorldClim with lat/long
   In R studio with the script 3_WorldClim_data_extraction.R)

4. Cross validation for climate variable prediction accuracy
   In R studio with the script Cannabis_cross_validation_2.R
   
6. Genomic selection with rrBLUP
   In R studio with the script Cannabis_GS_WORKING.R
 

   Figures for the manuscript can be found in scripts/EGS_Figures/Cannabis_EGS_Figures_manuscript.R 


