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

# RR-BLUP ###
# k.xval function from external script, running 10-fold cross-validation with 50 repetitions
xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")

# Gaussian Kernel
K <- A.mat(g.in) #Compute relationship matrix using Gaussian Kernel
k_dist <- dist(K) #Calculate distance matrix from relationship matrix

#Run k-fold cross-validation for Gaussian Kernel
xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

#Exponential Kernel
#Compute relationship matrix using Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)

#Set row and column names for the Exponential Kernel matrix
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)

#Calculate distance matrix from Exponential Kernel relationship matrix
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

#Run k-fold cross-validation for Exponential Kernel
xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

#BayesCPi
#Run k-fold cross-validation for BayesCPi model
xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50, niter=3000,nburn=1200)
saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")


#####
#PLOTS
#####
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

# Gaussian Kernel - K-fold
#Load cross-validation results for Gaussian Kernel
gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

# Exponential Kernel - K-fold
#Load cross-validation results for Exponential Kernel
EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - K-fold
#Load cross-validation results for BayesCpi model
#bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
#bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
#bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"
gauss_kfold_10$model <- "Gaussian Kernel"
EXP_kfold_10$model <- "Exponential Kernel"
bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"
gauss_kfold_10$xval <- "Ten-Fold"
EXP_kfold_10$xval <- "Ten-Fold"
bayescpi_kfold_10$xval <- "Ten-Fold"


#Combine all model results into a single list
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

#Filter for specific traits of interest
all_bio <- all_models[all_models$trait %in% c('bio_01', 'bio_02', 'bio_03', 'bio_04', 'bio_05','bio_06',
                                              'bio_07','bio_08','bio_09','bio_010','bio_11',
                                              'bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))

ggplot(all_models, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)


