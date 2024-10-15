######################################################
# Step 5: Genomic Selection with RR-BLUP
######################################################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_core_n25_greedy.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1:4)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$Sample_ID

# Subset environmental data for analysis (start from the 5th column)
y.in.rr <- envdat[, 5:ncol(envdat)]  
minicore_entries <- which(envdat$Core == TRUE)  # Subset to core lines where 'Core' is TRUE
y.in.rr <- y.in.rr[minicore_entries, ]  
y.in.mat <- as.matrix(y.in.rr)
print(y.in.mat [1:10, 1:10])  # Visualize 


# Training set
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train, ]  # Subset rows (individuals) directly
print(g.train[1:10, 1:10])  # Visualize part of g.train for debugging
dim(g.train)


# Remove SNPs with any missing data in g.train
g.train <- g.train[, colSums(is.na(g.train)) == 0]
dim(g.train)

# Check for remaining NA values in the genotype matrix
if (any(is.na(g.train))) stop("There are still missing values in g.train after SNP removal.")

cat("Dimensions of g.train after removing SNPs with missing data:", dim(g.train), "\n")

# Subset g.pred to include only individuals in `pred`
pred <- setdiff(row.names(g.in), train)  # This gives 57 individuals not in the training set
g.pred <- g.in[pred, ]

# Remove SNPs with missing data in g.pred
g.pred <- g.pred[, colSums(is.na(g.pred)) == 0]

# Find common SNPs between g.train and g.pred
common_snps <- intersect(colnames(g.train), colnames(g.pred))

# Subset both g.train and g.pred to keep only the common SNPs
g.train <- g.train[, common_snps]
g.pred <- g.pred[, common_snps]

# Check dimensions after SNP alignment
cat("Dimensions of g.train:", dim(g.train), "\n")
cat("Dimensions of g.pred after SNP alignment:", dim(g.pred), "\n")

# Ensure row alignment between g.pred and g.in
if (!all(row.names(g.pred) == pred)) stop("Row names of g.pred do not match prediction set.")

cat("Length of train:", length(train), "\n")
cat("Length of pred:", length(pred), "\n")
cat("Length of row names of g.pred:", length(row.names(g.pred)), "\n")

# Check if row names of g.pred match pred exactly
if (!all(row.names(g.pred) == pred)) {
  cat("Mismatch in row names between g.pred and pred:\n")
  mismatched_rows <- setdiff(row.names(g.pred), pred)
  print(mismatched_rows)
} else {
  cat("Row names of g.pred match pred.\n")
}

# List of traits to analyze
traits <- colnames(y.in.mat)

# Initialize object for storing results
gebv_df <- data.frame(matrix(nrow = nrow(g.in), ncol = length(traits)))
colnames(gebv_df) <- traits


# RR-BLUP loop for each trait
for (t in 1:length(traits)) {
  trait <- traits[t]
  
  # Set up training set for the trait
  y.train <- as.matrix(y.in.mat[train, trait])  
  
  # Run RR-BLUP model
  solve.out <- mixed.solve(y = y.train, Z = g.train, SE = FALSE, return.Hinv = FALSE)
  u.hat <- solve.out$u
  
  # Calculate GEBVs for both prediction and training sets
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  
  # Store results in the combined dataframe
  pred_rows <- match(row.names(g.pred), row.names(g.in))  # Indices for pred rows
  train_rows <- match(train, row.names(g.in))  # Indices for train rows
  
  gebv_df[pred_rows, t] <- GEBV  # Predictions for test lines
  gebv_df[train_rows, t] <- GEBV_train  # Predictions for training lines
}

# Final output: GEBVs for all individuals
row.names(gebv_df) <- row.names(g.in)
write.csv(gebv_df, '~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis.csv')






