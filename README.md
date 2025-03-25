# Cannabis Environmental Genomic Selection

This repository provides the supporting data for a recently submitted manuscript entitled "Dissecting Genotype-Environment interactions with functional implications for parental selection in Cannabis Breeding"

A pre-print of this work can be found at: (to follow shortly)

High resolution verions of main and supplemental figures will be available on fig share https://figshare.com/authors/Anna_H_McCormick/17741367

# VCF Pipeline
Scripts used for .vcf generation are located in the Scripts folder 

# EGS Pipeline 

Inputs needed: .vcf file with your genotype information and lat/long data for your samples in a .csv.
The below steps detail the pipeline that was conducted on the Ren et al., 2021 dataset (PRJNA734114). This was similarly repeated on the fusion of the Ren et al., 2021 and Soorni et al., 2017 datasets in this manuscript.  

1. Identification of core (1_core_identification.py)
   This script was used on the generated vcf file (Cannabis_sativa_PRJNA734114_filtered.vcf.gz) in PyCharm to identify a core (change accordingly for 30% of your dataset)
   *You will need to edit this to your samples which have georeferences. If your .vcf is small an alternative is available in R and can be found in the Cannabis_EGS_Figures_manuscript.R file
   
2. Convert vcf to format for rrBLUP (2_vcf_to_rrBLUP.py)
   In terminal locally with the below command structure:
   python vcf_to_rrBLUP.py -i Cannabis_sativa_PRJNA734114_filtered.vcf.gz -o Ren_PRJNA734114_rrBLUP_format
   *may need edits if your vcf structure is different
   
3. Extract relevant climate data from WorldClim with lat/long
   In R studio with the script 3_WorldClim_data_extraction.R

4. Examining different models for climate variable prediction accuracy
   In R studio with the script 4_Cross_validation.R. If dataset is large recommend HPC compute. Requires import of xval_kfol_functions.R (Credit: Quinn Campbell)
   
5. Genomic selection with rrBLUP
   In R studio with the script 5_Genomic_Selection.R

R versions used in this work included 4.1.2 and 4.4.0
Figures made in R for this manuscript can be found in the Data folder under Cannabis_EGS_Figures_manuscript.R 


