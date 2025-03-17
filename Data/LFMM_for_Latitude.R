###########################
# LFMM for Latitude all
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

# Convert genotypes to numeric matrix (0 for ref, 1 for het, 2 for alt)
convert_to_numeric <- function(genotype, ref, alt) {
  if (genotype == paste(ref, ref, sep = "/")) return(0)
  else if (genotype == paste(ref, alt, sep = "/") || genotype == paste(alt, ref, sep = "/")) return(1)
  else if (genotype == paste(alt, alt, sep = "/")) return(2)
  else return(NA)
}

genotype_data <- t(apply(gt, 1, function(row, ref, alt) {
  sapply(1:length(row), function(i) convert_to_numeric(row[i], ref[i], alt[i]))
}, ref = ref_alleles, alt = alt_alleles))
genotype_data <- as.matrix(genotype_data)

# Impute missing values in genotype_data
genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
genotype_data <- t(genotype_data)

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Latitude_data <- as.matrix(env_data[ , "Latitude", drop = FALSE])
# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

Latitude_data[is.na(Latitude_data)] <- mean(Latitude_data, na.rm = TRUE)

#mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Latitude_data, K = 1)
mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Latitude_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Latitude_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Perform FDR adjustment on p-values
lfmm_results$fdr_p_value <- p.adjust(lfmm_results$p_value, method = "BH")

# Convert FDR-adjusted p-values to -log10 scale
lfmm_results$log_fdr_p <- -log10(lfmm_results$fdr_p_value)

# Define the FDR significance threshold 
fdr_threshold <- -log10(0.01)

# Plot with cumulative positions, 99th percentile line, FDR threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +  # 99th percentile
  geom_hline(yintercept = fdr_threshold, color = "purple", linetype = "dashed") +  # FDR threshold
  labs(
    title = "SNP Outliers for Latitude for All Groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p



###########################
# LFMM for Bio1 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

# Convert genotypes to numeric matrix (0 for ref, 1 for het, 2 for alt)
convert_to_numeric <- function(genotype, ref, alt) {
  if (genotype == paste(ref, ref, sep = "/")) return(0)
  else if (genotype == paste(ref, alt, sep = "/") || genotype == paste(alt, ref, sep = "/")) return(1)
  else if (genotype == paste(alt, alt, sep = "/")) return(2)
  else return(NA)
}

genotype_data <- t(apply(gt, 1, function(row, ref, alt) {
  sapply(1:length(row), function(i) convert_to_numeric(row[i], ref[i], alt[i]))
}, ref = ref_alleles, alt = alt_alleles))
genotype_data <- as.matrix(genotype_data)

# Impute missing values in genotype_data
genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
genotype_data <- t(genotype_data)


saveRDS(genotype_data, "~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")
#genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LRMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio1_data <- as.matrix(env_data[ , "bio_01", drop = FALSE])
Bio1_data[is.na(Bio1_data)] <- mean(Bio1_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio1_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio1_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 1 (Annual Mean Temperature) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio1_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



###########################
# LFMM for Bio2 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio2_data <- as.matrix(env_data[ , "bio_02", drop = FALSE])
Bio2_data[is.na(Bio2_data)] <- mean(Bio2_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio2_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio2_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 2 (Mean Diurnal Range) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio2_all.pdf", plot = p, device = "pdf", width = 15, height = 4)

###########################
# LFMM for Bio3 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio3_data <- as.matrix(env_data[ , "bio_03", drop = FALSE])
Bio3_data[is.na(Bio3_data)] <- mean(Bio3_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio3_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio3_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 3 (Isothermality) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio3_all.pdf", plot = p, device = "pdf", width = 15, height = 4)

################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio3.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio3.csv", row.names = FALSE)


###########################
# LFMM for Bio4 ALL
###########################




###########################
# LFMM for Bio5 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio5_data <- as.matrix(env_data[ , "bio_05", drop = FALSE])
Bio5_data[is.na(Bio5_data)] <- mean(Bio5_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio5_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio5_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 5 (Max Temperature of Warmest Month) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio5_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


###########################
# LFMM for Bio6 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio6_data <- as.matrix(env_data[ , "bio_06", drop = FALSE])
Bio6_data[is.na(Bio6_data)] <- mean(Bio6_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio6_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio6_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 6 (Min Temperature of Coldest Month) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
#ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio6_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
#write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio6.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################

library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
#write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio6.csv", row.names = FALSE)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio6_2.csv", row.names = FALSE)


###########################
# LFMM for Bio7 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio7_data <- as.matrix(env_data[ , "bio_07", drop = FALSE])
Bio7_data[is.na(Bio7_data)] <- mean(Bio7_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio7_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio7_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 7 (Temperature Annual Range) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio7_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


###########################
# LFMM for Bio8 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio8_data <- as.matrix(env_data[ , "bio_08", drop = FALSE])
Bio8_data[is.na(Bio8_data)] <- mean(Bio8_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio8_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio8_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 8 (Mean Temperature of Wettest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio8_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


###########################
# LFMM for Bio9 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

# Convert genotypes to numeric matrix (0 for ref, 1 for het, 2 for alt)
convert_to_numeric <- function(genotype, ref, alt) {
  if (genotype == paste(ref, ref, sep = "/")) return(0)
  else if (genotype == paste(ref, alt, sep = "/") || genotype == paste(alt, ref, sep = "/")) return(1)
  else if (genotype == paste(alt, alt, sep = "/")) return(2)
  else return(NA)
}

genotype_data <- t(apply(gt, 1, function(row, ref, alt) {
  sapply(1:length(row), function(i) convert_to_numeric(row[i], ref[i], alt[i]))
}, ref = ref_alleles, alt = alt_alleles))
genotype_data <- as.matrix(genotype_data)

# Impute missing values in genotype_data
genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
genotype_data <- t(genotype_data)


env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio9_data <- as.matrix(env_data[ , "bio_09", drop = FALSE])
Bio9_data[is.na(Bio9_data)] <- mean(Bio9_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio9_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio9_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 9 (Mean Temperature of the Driest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio9_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio9.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio9.csv", row.names = FALSE)


###########################
# LFMM for Bio10 ALL
###########################
env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio10_data <- as.matrix(env_data[ , "bio_10", drop = FALSE])
Bio10_data[is.na(Bio10_data)] <- mean(Bio10_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio10_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio10_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 10 (Mean Temperature of Warmest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio10_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



###########################
# LFMM for Bio11 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio11_data <- as.matrix(env_data[ , "bio_11", drop = FALSE])
Bio11_data[is.na(Bio11_data)] <- mean(Bio11_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio11_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio11_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 11 (Mean Temperature of Coldest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
#ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio11_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
#write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio11.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
#write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio11.csv", row.names = FALSE)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio11_2.csv", row.names = FALSE)


##################################################################################################################################################################
###########################
# LFMM for Bio12 ALL
###########################

library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio12_data <- as.matrix(env_data[ , "bio_12", drop = FALSE])
Bio12_data[is.na(Bio12_data)] <- mean(Bio12_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio12_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio12_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 12 (Annual Precipitation) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio12_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio12.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio12.csv", row.names = FALSE)


###########################
# LFMM for Bio13 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio13_data <- as.matrix(env_data[ , "bio_13", drop = FALSE])
Bio13_data[is.na(Bio13_data)] <- mean(Bio13_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio13_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio13_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 13 (Precipitation of the Wettest Month) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio13_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio13.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio13.csv", row.names = FALSE)



###########################
# LFMM for Bio14 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

# Convert genotypes to numeric matrix (0 for ref, 1 for het, 2 for alt)
convert_to_numeric <- function(genotype, ref, alt) {
  if (genotype == paste(ref, ref, sep = "/")) return(0)
  else if (genotype == paste(ref, alt, sep = "/") || genotype == paste(alt, ref, sep = "/")) return(1)
  else if (genotype == paste(alt, alt, sep = "/")) return(2)
  else return(NA)
}

genotype_data <- t(apply(gt, 1, function(row, ref, alt) {
  sapply(1:length(row), function(i) convert_to_numeric(row[i], ref[i], alt[i]))
}, ref = ref_alleles, alt = alt_alleles))
genotype_data <- as.matrix(genotype_data)

# Impute missing values in genotype_data
genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
genotype_data <- t(genotype_data)


env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio14_data <- as.matrix(env_data[ , "bio_14", drop = FALSE])
Bio14_data[is.na(Bio14_data)] <- mean(Bio14_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio14_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio14_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 14 (Precipitation of the Driest Month) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
#ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio14_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
#write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio14.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio14_2.csv", row.names = FALSE)


###########################
# LFMM for Bio15 ALL
###########################

###########################
# LFMM for Bio16 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio16_data <- as.matrix(env_data[ , "bio_16", drop = FALSE])
Bio16_data[is.na(Bio16_data)] <- mean(Bio16_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio16_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio16_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 16 (Precipitation of the Wettest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio16_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio16.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio16.csv", row.names = FALSE)


###########################
# LFMM for Bio17 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio17_data <- as.matrix(env_data[ , "bio_17", drop = FALSE])
Bio17_data[is.na(Bio17_data)] <- mean(Bio17_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio17_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio17_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 17 (Precipitation of the Driest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio17_all.pdf", plot = p, device = "pdf", width = 15, height = 4)



################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio17.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio17.csv", row.names = FALSE)

###########################
# LFMM for Bio18 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

# Convert genotypes to numeric matrix (0 for ref, 1 for het, 2 for alt)
convert_to_numeric <- function(genotype, ref, alt) {
  if (genotype == paste(ref, ref, sep = "/")) return(0)
  else if (genotype == paste(ref, alt, sep = "/") || genotype == paste(alt, ref, sep = "/")) return(1)
  else if (genotype == paste(alt, alt, sep = "/")) return(2)
  else return(NA)
}

genotype_data <- t(apply(gt, 1, function(row, ref, alt) {
  sapply(1:length(row), function(i) convert_to_numeric(row[i], ref[i], alt[i]))
}, ref = ref_alleles, alt = alt_alleles))
genotype_data <- as.matrix(genotype_data)

# Impute missing values in genotype_data
genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
genotype_data <- t(genotype_data)


env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio18_data <- as.matrix(env_data[ , "bio_18", drop = FALSE])
Bio18_data[is.na(Bio18_data)] <- mean(Bio18_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio18_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio18_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 18 (Precipitation of the Warmest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
#ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio18_all.pdf", plot = p, device = "pdf", width = 15, height = 4)


################################
# Filter for SNPs that pass the Bonferroni threshold
################################
significant_snps <- lfmm_results %>%
  filter(log_p >= bonferroni_threshold)

# Display significant SNPs for reference
print(significant_snps)

# Save the significant SNPs to a CSV file
#write.csv(significant_snps, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Significant_SNPs_Bonferroni_bio18.csv", row.names = FALSE)

################################
# Find where SNPs fall in annotation
################################
library(GenomicRanges)
library(rtracklayer)

# Convert significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = significant_snps$chromosome,
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos),
  mcols = significant_snps
)

# Import GFF file as GRanges object
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Find overlaps between SNPs and annotation
overlaps <- findOverlaps(snp_gr, annotation_gr)

# Extract overlapping SNPs and annotation features
overlapping_snps <- snp_gr[queryHits(overlaps)]
overlapping_features <- annotation_gr[subjectHits(overlaps)]

# Determine column names for SNP ID, feature ID, and feature type
snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]

# Combine information into a data frame
snp_annotation_df <- data.frame(
  SNP_ID = mcols(overlapping_snps)[[snp_col_name]],
  SNP_Chromosome = as.character(seqnames(overlapping_snps)),
  SNP_Position = start(overlapping_snps),
  Feature_ID = mcols(overlapping_features)[[feature_id_col]],
  Feature_Type = mcols(overlapping_features)[[feature_type_col]]
)

# View the first few rows of the result and save to CSV
head(snp_annotation_df)
write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_bio18_2.csv", row.names = FALSE)



###########################
# LFMM for Bio19 ALL
###########################
library(lfmm)
library(vcfR)
library(dplyr)
library(ggplot2)

# Load and Prepare Genotype Data from VCF
vcf <- read.vcfR("~/R/cannabis_GEAV/Ren_Soorni/GEA/merged_Ren_Soorni_snps.vcf.gz")
gt <- extract.gt(vcf, return.alleles = TRUE)
ref_alleles <- getREF(vcf)
alt_alleles <- getALT(vcf)

genotype_data <- readRDS("~/R/cannabis_GEAV/Ren_Soorni/GEA/LFMM_genotype_data.rds")

env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", header = TRUE, stringsAsFactors = FALSE)
Bio19_data <- as.matrix(env_data[ , "bio_19", drop = FALSE])
Bio19_data[is.na(Bio19_data)] <- mean(Bio19_data, na.rm = TRUE)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

mod.lfmm <- lfmm_ridge(Y = genotype_data, X = Bio19_data, K = 5)
pv <- lfmm_test(Y = genotype_data, X = Bio19_data, lfmm = mod.lfmm)
lfmm_results <- data.frame(SNP = 1:ncol(genotype_data), p_value = pv$pvalue[, 1], log_p = -log10(pv$pvalue[, 1]))

# Extract position and chromosome info
lfmm_results <- cbind(lfmm_results, chromosome = as.factor(vcf@fix[, "CHROM"]), pos = as.numeric(vcf@fix[, "POS"]))

# Define chromosome lengths based on reference for correct cumulative positioning
chromosome_lengths <- c(
  "NC_044371.1" = 101209240,
  "NC_044375.1" = 96346938,
  "NC_044372.1" = 94670641,
  "NC_044373.1" = 91913879,
  "NC_044374.1" = 88181582,
  "NC_044377.1" = 79335105,
  "NC_044378.1" = 71238074,
  "NC_044379.1" = 64622176,
  "NC_044376.1" = 61561104,
  "NC_044370.1" = 104987320
)
cumulative_starts <- c(0, cumsum(chromosome_lengths[-length(chromosome_lengths)]))

# Map chromosomes and compute cumulative positions
chromosome_mapping <- setNames(1:10, names(chromosome_lengths))
lfmm_results <- lfmm_results %>%
  mutate(
    chromosome_numeric = chromosome_mapping[as.character(chromosome)],
    cumulative_pos = pos + cumulative_starts[chromosome_numeric]
  ) %>%
  arrange(chromosome_numeric, cumulative_pos)

# Define significance thresholds
significance_threshold <- -log10(0.01)  # 99th percentile threshold
bonferroni_threshold <- -log10(0.05 / nrow(lfmm_results))  # Bonferroni threshold


#asthetics
library(RColorBrewer)
library(MetBrewer)

# Calculate midpoints for chromosome labels
midpoints <- cumulative_starts + chromosome_lengths / 2

# Define colors using the "VanGogh1" palette from MetBrewer
chromosome_colors <- met.brewer("VanGogh1", n = length(unique(lfmm_results$chromosome_numeric)))


# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
# Define updated chromosome labels with prefixes
chromosome_labels <- paste0(c(1:9, "10(X)"), "_", names(chromosome_mapping))

# Plot with cumulative positions, 99th percentile line, Bonferroni threshold, and customized labels
p <- ggplot(lfmm_results, aes(x = cumulative_pos, y = log_p, color = as.factor(chromosome_numeric))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = significance_threshold, color = "red", linetype = "dashed") +
  geom_hline(yintercept = bonferroni_threshold, color = "blue", linetype = "dashed") +
  labs(
    title = "SNP outliers for Bio 19 (Precipitation of the Coldest Quarter) for all groupings",
    x = "Cumulative Genomic Position",
    y = "-log10(p-value)",
    color = "Chromosome"
  ) +
  # Add lines to separate chromosomes
  geom_vline(xintercept = cumulative_starts[-1], color = "grey", linetype = "dotted", size = 0.5) + 
  # Adjust title size
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(size = 10)  # Adjust the title size here
  ) +
  scale_x_continuous(
    breaks = cumulative_starts + chromosome_lengths / 2,  # Center labels over chromosomes
    labels = chromosome_labels
  ) +
  scale_color_manual(values = chromosome_colors)

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Bio19_all.pdf", plot = p, device = "pdf", width = 15, height = 4)






