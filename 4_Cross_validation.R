######################################################
# Step 4: Cross-validation of Models for Prediction Accuracy
######################################################
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")

# Read genotype data in rrBLUP format
data <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)

# Read environmental data from CSV
envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_all.csv', head = TRUE)  # full environmental dataset

# Read training dataset from CSV (only core environmental data)
trainingset <- read.csv("~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_core_n25.csv", head = TRUE)

# Set the row names for the environmental and training datasets
row.names(envdat) <- envdat$Sample_ID
row.names(trainingset) <- trainingset$Sample_ID

# Filter genotype data to include only columns matching Sample_ID in envdat
gd2 <- data[, c("rs.", "allele", "chrom", "pos", colnames(data)[colnames(data) %in% envdat$Sample_ID])]

# Set the SNP IDs as the row names for the genotype data
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1:4)]  # Isolate genotype data
g.in <- as.matrix(g.in)
g.in <- t(g.in)  # Transpose to align individuals as rows

### Remove SNPs with any missing data ###
g.in <- g.in[, colSums(is.na(g.in)) == 0]  # Remove SNPs with missing data

# Check dimensions of the cleaned genotype matrix
cat("Dimensions of g.in after removing SNPs with missing data:", dim(g.in), "\n")

### Load Environmental Data ###
# Select unique identifier and environmental data columns for training set
y.trainset <- trainingset[, c(1, 5:(ncol(trainingset)))]  # Select unique identifier and environmental data only

# Select unique identifier and environmental data columns for the full dataset
y.in <- envdat[, c(1, 5:ncol(envdat))]

# Convert environmental data to matrix format for RR-BLUP
y.trainset.mat <- as.matrix(y.trainset[, -1])  # Exclude the first column (Sample_ID)
y.in.mat <- as.matrix(y.in[, -1])  # Exclude the first column (Sample_ID)

# Ensure row alignment between genotype and environmental data
if (!all(row.names(g.in) %in% row.names(y.in.mat))) {
  stop("Mismatch between row names of g.in and y.in.mat")
}

# Ensure the same for the training set
if (!all(row.names(g.in) %in% row.names(y.trainset.mat))) {
  stop("Mismatch between row names of g.in and y.trainset.mat")
}

# Check dimensions before cross-validation
cat("Dimensions of g.in:", dim(g.in), "\n")
cat("Dimensions of y.in.mat:", dim(y.in.mat), "\n")
cat("Dimensions of y.trainset.mat:", dim(y.trainset.mat), "\n")

### Run k-fold Cross-validation for RR-BLUP ###
# k.xval function from external script, running 10-fold cross-validation with 50 repetitions
xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)

# Save the cross-validation results
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")
