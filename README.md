# Cannabis Environmental Genomic Selection

This repository provides the supporting data for a recently submitted manuscript entitled "X"

R script Cannabis_EGS_Figures_manuscript.R contains analysis for main and supplemental figures

A pre-print of this work can be found at (to add shortly)

High resolution verions of main and supplemental figures are to be available on fig share once accepted https://figshare.com/authors/Anna_H_McCormick/17741367





# Pipeline Steps

1. Identification of core (Step1_core.py)
   This script was used on the filtered vcf file in PyCharm to identify a core (n=25)
   
2. Convert vcf to format for rrBLUP (vcf_to_rrBLUP.py)
   In terminal locally with the below command structure:
   python vcf_to_rrBLUP.py -i Cannabis_sativa_PRJNA734114_filtered.vcf.gz -o Ren_PRJNA734114_rrBLUP_format


3. Extract relevant climate data from WorldClim with lat/long (R studio - WorldClim_data_extraction.R)
 Input
  Output
5. Cross validation for climate variable prediction accuracy (Cannabis_cross_validation_2.R)
  Input
  Output
6. Genomic selection with rrBLUP (R studio - Cannabis_GS_WORKING.R)
  Input
  Output
   
