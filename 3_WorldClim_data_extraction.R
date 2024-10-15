########################################
#Step 1 - MOVED TO PYTHON (Step1_core.py)
#########################################
#Core collection identification from dataset
library(vcfR)
library(poppr)
library(adegenet)
library(corehunter)

vcf <- read.vcfR("~/R/cannabis_GEAV/Inputs/Cannabis_sativa_PRJNA734114_filtered.vcf.gz")

# Convert VCF to genind object and set ploidy
x <- vcfR2genind(vcf)
ploidy(x) <- 2

# Subset the genind object to include only the georeferenced 44 individuals
subset_ids <- c("SRR14708197", "SRR14708198", "SRR14708199", "SRR14708200", "SRR14708201", 
                "SRR14708202", "SRR14708203", "SRR14708204", "SRR14708205", "SRR14708206", 
                "SRR14708207", "SRR14708208", "SRR14708209", "SRR14708210", "SRR14708211", 
                "SRR14708212", "SRR14708213", "SRR14708217", "SRR14708222", "SRR14708228", 
                "SRR14708229", "SRR14708230", "SRR14708231", "SRR14708232", "SRR14708233", 
                "SRR14708234", "SRR14708235", "SRR14708236", "SRR14708237", "SRR14708240", 
                "SRR14708243", "SRR14708244", "SRR14708246", "SRR14708248", "SRR14708251", 
                "SRR14708254", "SRR14708258", "SRR14708259", "SRR14708272", "SRR14708273", 
                "SRR14708275", "SRR14708276", "SRR14708277", "SRR14708278")

x_subset <- x[indNames(x) %in% subset_ids, ]

# Compute genetic distance matrix for the subset
x.dist <- poppr::bitwise.dist(x_subset)

# Sample a core collection from the subset
my.data <- distances(as.matrix(x.dist))
core <- sampleCore(my.data, size = 25)

# Save core collection to CSV
write.csv(core, '~/R/cannabis_GEAV/Outputs/cannabis_core_n25.csv')

######################################################
#Step 2 (vcf_to_rrBLUP.py)
######################################################
#Converting .vcf to -1/0/1 format for rrBLUP

library(vcfR)
vcf_file <- "~/R/cannabis_GEAV/Inputs/Cannabis_sativa_PRJNA734114_filtered.vcf.gz"
vcf_data <- read.vcfR(vcf_file)

# Extract the genotype data
gt_matrix <- extract.gt(vcf_data)

# Custom function to convert genotypes to -1/0/1
convert_to_rrBLUP <- function(gt) {
  if (is.na(gt) || gt == "." || gt == "./.") {
    return(NA)  # Handle missing data
  } else if (gt == "0/0" || gt == "0|0") {
    return(-1)
  } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
    return(0)
  } else if (gt == "1/1" || gt == "1|1") {
    return(1)
  } else {
    return(NA)  # Return NA for any other unexpected format
  }
}

# Apply the conversion
geno_converted <- apply(gt_matrix, 1:2, convert_to_rrBLUP)

# Save the result
write.table(geno_converted, "~/R/cannabis_GEAV/Outputs/Ren_PRJNA734114_rrBLUP_format.txt", sep = "\t", col.names = NA, quote = FALSE)

######################################################
#Step 3
######################################################
#Extract relevant environmental data from WorldClim for latitude and longitudes 

library(raster)
worldclim_path <- "~/R/Data/climate_data/wc2.0_30s_bio/"
bio_files <- paste0(worldclim_path, "bio_", sprintf("%02d", 1:19), ".tif")
climate_layers <- stack(bio_files)

# Load your latitude and longitude dataset
data <- read.csv("~/R/cannabis_GEAV/Inputs/ren_lat_long.csv")
head(data)

# Prepare the coordinates (assuming Latitude and Longitude columns are present)
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data.csv", row.names = FALSE)

######################################################
#Step 4 
######################################################
# cross validation of models for prediction accuracy

library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

#load functions 
#source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") 
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R") 

#Read genotype data in rrBLUP format
data <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
colnames(data)

#Read environmental data from .csv
envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_all.csv', head = T) # full environmental dataset

#Read training dataset from .csv
trainingset <- read.csv("~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_core_n25.csv", head = T) #only have env data for core in this




# Filter the rrBLUP-formatted genotype data to include only columns that match Sample_ID in envdat
gd2 <- data[, c("rs.", "allele", "chrom", "pos", colnames(data)[colnames(data) %in% envdat$Sample_ID])]

# Set the rs. as the rownames of gd2 (use SNP IDs as row names)
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1:4)]  # This removes the first four columns to isolate genotype data

# Convert the genotype matrix to a matrix type
g.in <- as.matrix(g.in)
g.in <- t(g.in)



### Load Environmental Data ###
row.names(envdat) <- envdat$Sample_ID # set gen_id as rownames
row.names(trainingset) <- trainingset$Sample_ID



#Select unique identifier and environmental data columns for training set
y.trainset <- trainingset[,c(1,5:(ncol(trainingset)))] # select unique identifier and environmental data only 

#Select unique identifier and environmental data columns for the full dataset
y.in <- envdat[,c(1,5:ncol(envdat))] 





#RR-BLUP
#Prepare data for rrBLUP model
y.trainset.rr <- trainingset[,c(5:(ncol(trainingset)))]
y.in.rr <- envdat[,5:ncol(envdat)]

#Convert data to matrix format
y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)


#Run k-fold cross-validation for rrBLUP
dim(g.in)
dim(y.in.mat)
dim(y.trainset.mat)


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
setwd("~/R/barley_collab/parental_selection/cross_validation/")

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
bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)


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
all_bio <- all_models[all_models$trait %in% c('bio1', 'bio2', 'bio3', 'bio4', 'bio5','bio6',
                                              'bio7','bio8','bio9','bio10','bio11',
                                              'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))
ggplot(all_bio, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)






######################################################
#Step 5 
######################################################
# Genomic Selection with rrBLUP
library(rrBLUP)
library(dplyr)


gd1 <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
colnames(gd1)

#Read environmental data from .csv
envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_core_n25.csv', head = T) 
colnames(envdat)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.


# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1:4)]  # Removing 'rs.', 'allele', 'chrom', 'pos'
g.in <- as.matrix(gd3)  # Convert to matrix

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$Sample_ID


# Subset environmental data starting from the 5th column (bio_01 to bio_19)
y.in.rr <- envdat[, 5:ncol(envdat)]  # Environmental data starts at column 5

# Subset to core lines where 'Core' is TRUE
minicore_entries <- which(envdat$Core == TRUE)
y.in.rr <- y.in.rr[minicore_entries, ]  # Only select rows where 'Core' is TRUE

# Remove 'Core' column as it's not needed in the analysis
#y.in.rr <- y.in.rr[, -1]  # Remove 'Core' column

# Convert the environmental data to a matrix
y.in.mat <- as.matrix(y.in.rr)


# Set up training and genotype data
train <- row.names(y.in.mat) # Names of training lines
#g.train <- g.in[train,]  # Set training genotypes
#g.train <- g.in[, match(train, colnames(g.in))]
#g.train <- t(g.in[match(train, colnames(g.in)), ])  # Transpose the result to get rows as individuals
g.train <- g.in[, match(train, colnames(g.in))]  # Subset columns of g.in for training
g.train <- t(g.train)
dim(g.train) 


# List of traits to analyze
traits <- colnames(y.in.mat)

# Prediction setup
pred <- setdiff(row.names(g.in), train)  # Names of lines not in training set
#g.pred <- g.in[pred,]
#g.pred <- t(g.in[setdiff(colnames(g.in), train), ])  # Transpose for prediction
g.pred <- g.in[, setdiff(colnames(g.in), train)]  # Subset columns of g.in for prediction (samples not in training)

# Initialize objects for storing results
marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))  # Initialize with all lines (train + pred)
colnames(gebv_df) <- traits

# Check dimensions
dim(g.train)  # Should match the number of training samples
dim(y.train)
dim(g.pred)   # Should match the number of prediction samples




# RR-BLUP loop for each trait
for(t in 1:length(traits)) {
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train, trait])  # Training set for the trait
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training set
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Use match to ensure correct row alignment
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(row.names(g.train), row.names(g.in))  # Indices for train rows
  
  # Store the results in the combined dataframe
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Set row names for the result to match the line names from `g.in`
row.names(gebv_df) <- row.names(g.in)

write.csv(gebv_df, '~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis.csv')





######################################################
#Step 6 
######################################################
# Visualisation of the results
library(ggplot2)


data2 <-read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_plots.csv')

p1 <- ggplot(data2, aes(x = Cultivated_status, y = bio_01, fill = Cultivated_status)) +
  geom_boxplot(color = "black") +
  labs(title = "GEAV bio_01 - Annual Mean Temperature",
       x = "Cultivated Status",
       y = "bio_01") +
  scale_fill_manual(values = c("cultivar" = "green", 
                               "feral" = "orange", 
                               "landrace" = "blue")) +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")

print(p1)


data2 <-read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_plots.csv')
p1 <- ggplot(data2, aes(x = Phylogenetic_tree, y = bio_01, fill = Phylogenetic_tree)) +
  geom_boxplot(color = "black") +
  labs(title = "GEAV bio_01 - Annual Mean Temperature",
       x = "Cultivated Status",
       y = "bio_01") +
  scale_fill_manual(values = c("basal" = "orange", 
                               "hemp-type" = "green", 
                               "drug-type" = "red",
                               "drug-type feral" = "blue")) +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")

# Show plot
print(p1)
##Loop
# Load necessary libraries
library(ggplot2)

# Assuming your `bio` variables are named "bio_01" to "bio_19"
bio_vars <- paste0("bio_", sprintf("%02d", 1:19))  # Generate list of bio variables "bio_01" to "bio_19"

# Loop through each bio variable and generate a plot
for (bio_var in bio_vars) {
  
  # Create the ggplot for the current bio variable
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_var, fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = paste("GEAV", bio_var, "- Environmental Variable"),  # Dynamic title
         x = "Cultivated Status",
         y = bio_var) +
    scale_fill_manual(values = c("basal" = "orange", 
                                 "hemp-type" = "green", 
                                 "drug-type" = "red",
                                 "drug-type feral" = "blue")) +
    theme_minimal() +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Print each plot
  print(p)
  
  # Optionally save the plot to a file
  ggsave(paste0("~/R/cannabis_GEAV/Outputs/", bio_var, "_boxplot.pdf"), plot = p)
}



##########
# Load necessary libraries
library(ggplot2)

# Create a ggplot object
p <- ggplot(data2, aes(x = Phylogenetic_tree, y = bio_01, fill = Phylogenetic_tree)) +
  geom_boxplot(color = "black") +
  facet_grid(~factor(Cultivated_status), scales = "free_x", space = "free") +  # Separate by Cultivated Status with flexible space
  labs(title = "Phylogenetic Tree and Cultivated Status Comparison",
       x = NULL,
       y = "bio_01 - Annual Mean Temperature") +
  scale_fill_manual(values = c("basal" = "orange", 
                               "hemp-type" = "green", 
                               "drug-type" = "red", 
                               "drug-type feral" = "blue")) +  # Color mapping for Phylogenetic Tree
  theme_minimal() +
  theme(strip.background = element_blank(),  # Remove strip background for cleaner look
        panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
        axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        panel.grid.major = element_blank(),  # Optional: remove major grid lines
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = "longdash")

# Display the plot
print(p)




