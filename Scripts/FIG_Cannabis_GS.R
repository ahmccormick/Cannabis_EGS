############################
#Figure 1
###########################
library(ggplot2)
library(raster)
library(viridis)

# Load species distribution .tiff file and reduce its resolution for faster plotting
species_distribution <- raster("~annamccormick/R/species_distribution/results/mean_soil_AND_worldclim.tif")
species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # Adjust 'fact' as needed

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

# Load your CSV data containing latitude and longitude
data <- read.csv("~/R/cannabis_GEAV/DataSetMerged/inputs/ren_lw_lat_long.csv")

# Define the limits for latitude and longitude based on the summary
lat_range <- c(10, 80)  # Adjust based on your data
lon_range <- c(-20, 150)  # Adjust based on your data

p4 <- ggplot() + 
  geom_raster(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") + 
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Dataset), size = 1.5, alpha = 0.8) +
  #borders("world", colour = "gray85", fill = "gray80") +  # Add a base map
  coord_fixed(ratio = 1.3, xlim = lon_range, ylim = lat_range) +  # Zoom into area of interest
  labs(title = "Species Distribution with Sample Points", x = "Longitude", y = "Latitude") +
  scale_color_manual(values = c("Ren et al 2021" = "#FF5733", "LeafWorks" = "#33FF57")) +  # Set custom colors for Dataset
  theme_minimal() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
ggsave("~/R/cannabis_GEAV/DataSetMerged/outputs/species_distribution_with_points.pdf", plot = p4, width = 12, height = 5.39, device = "pdf")

############################
#Figure 2
###########################
#SNPRelate
library(SNPRelate)
setwd("~/R/cannabis_GEAV/MERGE/")
vcf.fn <- "~/R/cannabis_GEAV/MERGE/merged_shared_snps.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "lw_ren.gds", method="copy.num.of.ref")
snpgdsSummary("lw_ren.gds")
genofile <- snpgdsOpen("lw_ren.gds")
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")


#FIGURE 2A
library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")

# Dataset
p<-ggplot(pca, aes(x = EV1, y = EV2, color = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, alpha = 0.5) +  # Adds ellipses with 95% confidence interval
  labs(title = "PCA Plot by Dataset", x = "PC1", y = "PC2") +
  theme_minimal() +  # Custom colors
  labs(title = "Dataset Origins", 
       x = "PC1 (8.25%)",  # Adding percentage for PC1
       y = "PC2 (4.46%)") +  # Adding percentage for PC2
  theme_minimal()+
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_2A.pdf", plot = p, device = "pdf", width = 8, height = 6)

#FIGURE 2B
#Core 
pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")

# Plot PCA with core indicated by color
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Core)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses if needed
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +  # Color based on Core
  labs(title = "Core ", 
       x = "PC1 (8.25%)", 
       y = "PC2 (4.46%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_2B.pdf", plot = p, device = "pdf", width = 8, height = 6)

#FIGURE 2C
#PHYLO 
pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")

# Phylogenetic_tree
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Phylogenetic_tree)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type feral" = "blue", 
                                "drug-type" = "red", 
                                "unknown" = "grey")) +  # Assign colors
  labs(title = "Phylogenetic Tree Categorisation", 
       x = "PC1 (8.25%)", 
       y = "PC2 (4.46%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_2C.pdf", plot = p, device = "pdf", width = 8, height = 6)

############################
#Figure 3
###########################

#Figure 3A
#PHYLO GROUP
data2 <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_greedy.csv')

library(gridExtra)
library(ggplot2)

bio_vars <- sprintf("bio_%02d", 1:19)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type"))


# Define the selected bio variables and their corresponding labels
selected_bio_vars <- c("bio_04", "bio_15", "bio_05", "bio_06", "bio_16", "bio_17")
selected_labels <- c(
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter"
)

# Initialize an empty list to hold the selected plots
selected_plots <- list()

# Loop through the selected bio variables and generate the plots
for (i in seq_along(selected_bio_vars)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = selected_bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = selected_labels[i], x = NULL, y = selected_bio_vars[i]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 8)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add the plot to the list
  selected_plots[[i]] <- p
}

# Save the selected plots into a PDF
pdf("~/R/cannabis_GEAV/Outputs/GEAV_selected_bio_plots_phylogenetic.pdf", width = 8, height = 10)
grid.arrange(grobs = selected_plots, ncol = 2)  # Arrange the plots in a grid
dev.off()


#Figure 3B
#PHYLO GROUP
# Load the new dataset
data2 <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/rrblup_GEBV_values_REN_LW.csv')

library(gridExtra)
library(ggplot2)

# Define the bio variables (you already have them in your dataset)
bio_vars <- sprintf("bio_%02d", 1:19)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type"))

# Define the selected bio variables and their corresponding labels from your new associations
selected_bio_vars <- c("bio_04", "bio_15", "bio_05", "bio_06", "bio_16", "bio_17")
selected_labels <- c(
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter"
)

# Initialize an empty list to hold the selected plots
selected_plots <- list()

# Loop through the selected bio variables and generate the plots
for (i in seq_along(selected_bio_vars)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = selected_bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = selected_labels[i], x = NULL, y = selected_bio_vars[i]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue", "unknown" = "grey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 8)) +  # Adjust title size to 16
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add the plot to the list
  selected_plots[[i]] <- p
}

# Save the selected plots into a PDF
pdf("~/R/cannabis_GEAV/DataSetMerged/outputs/GEAV_selected_bio_plots_REN_LW_phylogenetic.pdf", width = 8, height = 10)
grid.arrange(grobs = selected_plots, ncol = 2)  # Arrange the plots in a grid with 2 columns
dev.off()



############################
#Figure X
###########################
#MG Chemistry 













############################
#Supplemental Figure 1
###########################
#Recapitulation of Ren dataset

#SNPRelate
library(SNPRelate)
setwd("~/R/cannabis_GEAV/Outputs/")

vcf.fn <- "~/R/cannabis_GEAV/Inputs/Cannabis_sativa_PRJNA734114_filtered.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "ren.gds", method="copy.num.of.ref")
snpgdsSummary("ren.gds")
genofile <- snpgdsOpen("ren.gds")
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/cannabis_GEAV/Outputs/ren_PC.csv")

#####################
# PCA 
#####################
library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/Outputs/ren_PC.csv")

# Phylogenetic_tree
p<-ggplot(pca, aes(x = EV1, y = EV2, color = Phylogenetic_tree)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type feral" = "blue", 
                                "drug-type" = "red", 
                                "unknown" = "grey")) +  # Assign colors
  labs(title = "Phylogenetic Tree Categorisation", 
       x = "PC1 (3.33%)", 
       y = "PC2 (2.68%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )

ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S1A.pdf", plot = p, device = "pdf", width = 8, height = 6)

########Core 
library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/Outputs/ren_PC.csv")

# Plot PCA with core indicated by color
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Core)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses if needed
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +  # Color based on Core
  labs(title = "PCA Plot with Core Indicated", 
       x = "PC1 (3.33%)", 
       y = "PC2 (2.68%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S1B.pdf", plot = p, device = "pdf", width = 8, height = 6)

#  Cultivated status
p<-ggplot(pca, aes(x = EV1, y = EV2, color = Cultivated_status)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, alpha = 0.5) +  # Adds ellipses with 95% confidence interval
  scale_color_manual(values = c("feral" = "orange", 
                                "landrace" = "blue", 
                                "cultivar" = "green")) +  # Custom colors
  labs(title = "Cultivated Status Categorisation", 
       x = "PC1 (3.33%)",  # Adding percentage for PC1
       y = "PC2 (2.68%)") +  # Adding percentage for PC2
  theme_minimal()+
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S1C.pdf", plot = p, device = "pdf", width = 8, height = 6)



############################
#Supplemental Figure 2
###########################
#fastSTRUCTURE

library(MetBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggnewscale)

##############
#K2
##############
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 6)]  # Colors for K1 and K2

# Load your data
test <- read.csv("~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/canna_k2.csv")

# Sort by K1 and K2 percentages
test <- test %>%
  arrange(desc(K1), desc(K2))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type
type_colors <- c("basal" = "orange", 
                 "hemp_type" = "green", 
                 "drug_type_feral" = "blue", 
                 "drug_type" = "red", 
                 "unknown" = "grey")

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/fastSTRUCTURE_plot_k2.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)



##############
#K3
##############
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(MetBrewer)

hiroshige_palette <- met.brewer("Hiroshige", 8)
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Load your data
test <- read.csv("~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/canna_k3.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 ~ "K1",
    K2 >= K1 & K2 >= K3 ~ "K2",
    K3 >= K1 & K3 >= K2 ~ "K3"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("basal" = "orange", 
                 "hemp_type" = "green", 
                 "drug_type_feral" = "blue", 
                 "drug_type" = "red", 
                 "unknown" = "grey")

# Hiroshige palette for K categories
k_colors <- hiroshige_palette[c(1, 4, 6)]

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/fastSTRUCTURE_plot_k3.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)


##############
#K4
##############
#load fastSTRUCTURE data

hiroshige_palette <- met.brewer("Hiroshige", 8)
#selected_colors <- hiroshige_palette[c(1, 4,6, 8)]
k_colors <- hiroshige_palette[c(1,4,6,8)]  # Colors for K

test <- read.csv("~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/canna_k4.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 ~ "K4"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("basal" = "orange", 
                 "hemp_type" = "green", 
                 "drug_type_feral" = "blue", 
                 "drug_type" = "red", 
                 "unknown" = "grey")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 4, 6, 8)]  # Choose appropriate colors for K1 to K4

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/fastSTRUCTURE_plot_k4.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)



##############
#K5
##############

test <- read.csv("~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/canna_k5.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 ~ "K5"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("basal" = "orange", 
                 "hemp_type" = "green", 
                 "drug_type_feral" = "blue", 
                 "drug_type" = "red", 
                 "unknown" = "grey")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 8)]  # Choose appropriate colors for K1 to K5

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/fastSTRUCTURE_plot_k5.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)


##############
#K6
##############
#load fastSTRUCTURE data

# Load your data
test <- read.csv("~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/canna_k6.csv")

# Create a new column for the dominant K category
test <- test %>%
  mutate(dominant_K = case_when(
    K1 >= K2 & K1 >= K3 & K1 >= K4 & K1 >= K5 & K1 >= K6 ~ "K1",
    K2 >= K1 & K2 >= K3 & K2 >= K4 & K2 >= K5 & K2 >= K6 ~ "K2",
    K3 >= K1 & K3 >= K2 & K3 >= K4 & K3 >= K5 & K3 >= K6 ~ "K3",
    K4 >= K1 & K4 >= K2 & K4 >= K3 & K4 >= K5 & K4 >= K6 ~ "K4",
    K5 >= K1 & K5 >= K2 & K5 >= K3 & K5 >= K4 & K5 >= K6 ~ "K5",
    K6 >= K1 & K6 >= K2 & K6 >= K3 & K6 >= K4 & K6 >= K5 ~ "K6"
  ))

# Sort by the dominant K category and then by the actual values in descending order
test <- test %>%
  arrange(dominant_K, desc(K1), desc(K2), desc(K3), desc(K4), desc(K5), desc(K6))

# Convert to long format after sorting
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5", "K6"), names_to = "category", values_to = "value")

# Reorder the 'Name' factor based on the sorted data
long_data$Name <- factor(long_data$Name, levels = unique(test$Name))

# Assign specific colors to each Type and K category
type_colors <- c("basal" = "orange", 
                 "hemp_type" = "green", 
                 "drug_type_feral" = "blue", 
                 "drug_type" = "red", 
                 "unknown" = "grey")

# Hiroshige palette for K categories (adjust the palette if needed)
k_colors <- hiroshige_palette[c(1, 3, 5, 6, 7, 8)]  # Choose appropriate colors for K1 to K6

# Plot with the colored bar for Type underneath
final_plot <- ggplot() +
  # Plot the K1 and K2 bars with their own fill colors
  geom_bar(data = long_data, aes(x = Name, y = value, fill = category), stat = "identity") +
  scale_fill_manual(name = "K Categories", values = k_colors) +  # Colors for K1 and K2
  
  # Add the colored bars underneath for Type using a separate geom_tile
  new_scale_fill() +  # Necessary to separate the two scales
  geom_tile(data = test, aes(x = Name, y = -0.05, fill = Type), height = 0.1) +
  scale_fill_manual(name = "Type", values = type_colors) +  # Colors for Type
  
  labs(title = "", x = "", y = "Percent Identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 1),  # Make the text smaller
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),       # Smaller legend text
    legend.title = element_text(size = 8),      # Smaller legend title
    legend.key.size = unit(0.4, 'cm'),          # Smaller legend keys
    legend.spacing.y = unit(0.2, 'cm')          # Reduce vertical spacing between legend items
  )

print(final_plot)

# Save the plot as a PDF with reduced height
ggsave(filename = "~/R/cannabis_GEAV/DataSetMerged/fastSTRUCTURE/fastSTRUCTURE_plot_k6.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)


#########################################
#SUPPLEMENTARY FIGURE 2B
#########################################
library("FactoMineR")
library("factoextra")

pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")
#Select 3rd and 4th principal components
pca2 <- pca[3:4]
#pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a
ggsave("~/R/cannabis_GEAV/MERGE/LW_Ren_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b
ggsave("~/R/cannabis_GEAV/MERGE/LW_Ren_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")


pca <- read.csv("~/R/cannabis_GEAV/Outputs/ren_PC.csv")
pca2 <- pca[3:4]
#pca2 <- pca[3:6]

my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a
ggsave("Ren_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b
ggsave("Ren_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")



############################
#Supplemental Figure 3
###########################
#Ren + LW Phylogeny









############################
#Supplemental Figure 4
###########################
#Cross Validation Ren dataset
#ON HPC
library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

#load functions 
#source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") 
source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/xval_kfold_functions.R") 
#source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R") 

#Read genotype data in rrBLUP format
data <- read.delim("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
colnames(data)

#Read environmental data from .csv
envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/ren_n44_extracted_climate_data_all.csv', head = T) # full environmental dataset

#Read training dataset from .csv
trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/ren_n44_extracted_climate_data_core_n25.csv", head = T) #only have env data for core in this


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
#bayescpi_kfold_10$model <- "BayesCpi"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"
gauss_kfold_10$xval <- "Ten-Fold"
EXP_kfold_10$xval <- "Ten-Fold"
#bayescpi_kfold_10$xval <- "Ten-Fold"

#gauss_kfold_10 <- gauss_kfold_10[-1, ]  # Removing the first row
#EXP_kfold_10 <- EXP_kfold_10[-1, ]      # Adjust as needed, depending on where the extra row is


#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

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
                                              'bio_07','bio_08','bio_09','bio_10','bio_11',
                                              'bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))

p<- ggplot(all_models, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)
p
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S4.pdf", plot = p, device = "pdf", width = 18, height = 10)
rrblup_results <- all_bio[all_bio$model == "rrBLUP", ]

# Export the results to a CSV file
write.csv(rrblup_results, "~/R/cannabis_GEAV/MERGE/FIGS/rrblup_prediction_accuracy_ren.csv", row.names = FALSE)


############################
#Supplemental Figure 5
###########################
#Cross Validation Ren + LW dataset
#ON HPC
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")

# Read genotype data in rrBLUP format
data <- read.delim("~/R/cannabis_GEAV/DataSetMerged/outputs/Ren_LW_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(data)[1] <- "rs." #needed to rename the sample column here
colnames(data) <- gsub("^X", "", colnames(data)) #needed to remove the extra X added in 

# Read environmental data from CSV
envdat <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/ren_lw_n50_extracted_climate_data.csv', head = TRUE)  # full environmental dataset

# Read training dataset from CSV (only core environmental data)
trainingset <- read.csv("~/R/cannabis_GEAV/DataSetMerged/outputs/ren_lw_n50_extracted_climate_data.csv", head = TRUE)

# Set the row names for the environmental and training datasets
row.names(envdat) <- envdat$Sample_ID
row.names(trainingset) <- trainingset$Sample_ID


# Filter genotype data to include only columns matching Sample_ID in envdat
colnames(data)
gd2 <- data[, c("rs.", colnames(data)[colnames(data) %in% envdat$Sample_ID])]

# Set the SNP IDs as the row names for the genotype data
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1)]  # Isolate genotype data
g.in <- as.matrix(g.in)
g.in <- t(g.in)  # Transpose to align individuals as rows

### Remove SNPs with any missing data ###
#g.in <- g.in[, colSums(is.na(g.in)) == 0]  # Remove SNPs with missing data

# Impute missing values with the column mean, processing column by column
for (i in 1:ncol(g.in)) {
  g.in[is.na(g.in[, i]), i] <- mean(g.in[, i], na.rm = TRUE)
}

# Check dimensions after imputation
cat("Dimensions of g.in after imputation:", dim(g.in), "\n")


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
setwd("~/R/cannabis_GEAV/DataSetMerged/cross_validation_results/")

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

#gauss_kfold_10 <- gauss_kfold_10[-1, ]  # Removing the first row
#EXP_kfold_10 <- EXP_kfold_10[-1, ]      # Adjust as needed, depending on where the extra row is


#Combine all model results into a single list
model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10, bayescpi_kfold_10)
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10)

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
                                              'bio_07','bio_08','bio_09','bio_10','bio_11',
                                              'bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19'),]

#Plot 
all_bio$trait <- factor(all_bio$trait, levels = paste0("bio", 1:19))

p<-ggplot(all_models, aes(y = r.mean, x = model, color = model)) +
  theme_bw() +
  geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1)

p
ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S5.pdf", plot = p, device = "pdf", width = 18, height = 10)

# Filter the results for the rrBLUP model and the bio variables
rrblup_results <- all_bio[all_bio$model == "rrBLUP", ]

# Export the results to a CSV file
write.csv(rrblup_results, "~/R/cannabis_GEAV/MERGE/FIGS/rrblup_prediction_accuracy.csv", row.names = FALSE)


############################
#Supplemental Figure 6
###########################
#REN ALONE
#Figure 6A and B (Temp and ppt)

#PHYLO GROUP
data2 <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_greedy.csv')

library(gridExtra)
library(ggplot2)

bio_vars <- sprintf("bio_%02d", 1:19)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type"))


# Define temperature and precipitation associations separately
temperature_associations <- c(
  "BIO 1 - Annual Mean Temperature",
  "BIO 2 - Mean Diurnal Range (Max temp - Min temp)",
  "BIO 3 - Isothermality (BIO2/BIO7) (×100)",
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 7 - Temperature Annual Range (BIO5-BIO6)",
  "BIO 8 - Mean Temperature of Wettest Quarter",
  "BIO 9 - Mean Temperature of Driest Quarter",
  "BIO 10 - Mean Temperature of Warmest Quarter",
  "BIO 11 - Mean Temperature of Coldest Quarter"
)

precipitation_associations <- c(
  "BIO 12 - Annual Precipitation",
  "BIO 13 - Precipitation of Wettest Month",
  "BIO 14 - Precipitation of Driest Month",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter",
  "BIO 18 - Precipitation of Warmest Quarter",
  "BIO 19 - Precipitation of Coldest Quarter"
)

# Create a list to hold the temperature and precipitation plots separately
temperature_plots <- list()
precipitation_plots <- list()

# Loop through each bio variable to generate the plots for temperature variables
for (i in 1:length(temperature_associations)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = temperature_associations[i], x = NULL, y = bio_vars[i]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to temperature plots
  temperature_plots[[i]] <- p
}

# Loop through each bio variable to generate the plots for precipitation variables
for (i in 1:length(precipitation_associations)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_vars[i + length(temperature_associations)], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = precipitation_associations[i], x = NULL, y = bio_vars[i + length(temperature_associations)]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to precipitation plots
  precipitation_plots[[i]] <- p
}

# Save combined temperature plots into one PDF
pdf("~/R/cannabis_GEAV/Outputs/GEAV_combined_temperature_plots_phylogenetic_4x.pdf", width = 8, height = 10)
grid.arrange(grobs = temperature_plots, ncol = 3)  # Arrange the temperature plots in a grid
dev.off()

# Save combined precipitation plots into another PDF
pdf("~/R/cannabis_GEAV/Outputs/GEAV_combined_precipitation_plots_phylogenetic_4x.pdf", width = 8, height = 10)
grid.arrange(grobs = precipitation_plots, ncol = 3)  # Arrange the precipitation plots in a grid
dev.off()



############################
#Supplemental Figure 7
###########################
#REN ALONE

#Figure 7A and B (Temp and ppt)
#CULT GROUP 
data2 <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_greedy.csv')

library(gridExtra)
library(ggplot2)
bio_vars <- sprintf("bio_%02d", 1:19)
# Reorder the factor levels for Cultivated_status
data2$Cultivated_status <- factor(data2$Cultivated_status, levels = c("feral", "landrace", "cultivar"))

# Define temperature and precipitation associations separately
temperature_associations <- c(
  "BIO 1 - Annual Mean Temperature",
  "BIO 2 - Mean Diurnal Range (Max temp - Min temp)",
  "BIO 3 - Isothermality (BIO2/BIO7) (×100)",
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 7 - Temperature Annual Range (BIO5-BIO6)",
  "BIO 8 - Mean Temperature of Wettest Quarter",
  "BIO 9 - Mean Temperature of Driest Quarter",
  "BIO 10 - Mean Temperature of Warmest Quarter",
  "BIO 11 - Mean Temperature of Coldest Quarter"
)

precipitation_associations <- c(
  "BIO 12 - Annual Precipitation",
  "BIO 13 - Precipitation of Wettest Month",
  "BIO 14 - Precipitation of Driest Month",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter",
  "BIO 18 - Precipitation of Warmest Quarter",
  "BIO 19 - Precipitation of Coldest Quarter"
)

# Create a list to hold the temperature and precipitation plots separately
temperature_plots <- list()
precipitation_plots <- list()

# Loop through each bio variable to generate the plots for temperature variables
for (i in 1:length(temperature_associations)) {
  p <- ggplot(data2, aes_string(x = "Cultivated_status", y = bio_vars[i], fill = "Cultivated_status")) +
    geom_boxplot(color = "black") +
    labs(title = temperature_associations[i], x = NULL, y = bio_vars[i]) +
    scale_fill_manual(values = c("feral" = "orange", "cultivar" = "green", "landrace" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to temperature plots
  temperature_plots[[i]] <- p
}

# Loop through each bio variable to generate the plots for precipitation variables
for (i in 1:length(precipitation_associations)) {
  p <- ggplot(data2, aes_string(x = "Cultivated_status", y = bio_vars[i + length(temperature_associations)], fill = "Cultivated_status")) +
    geom_boxplot(color = "black") +
    labs(title = precipitation_associations[i], x = NULL, y = bio_vars[i + length(temperature_associations)]) +
    scale_fill_manual(values = c("feral" = "orange", "cultivar" = "green", "landrace" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to precipitation plots
  precipitation_plots[[i]] <- p
}

# Save combined temperature plots into one PDF
pdf("~/R/cannabis_GEAV/Outputs/GEAV_combined_temperature_plots.pdf", width = 8, height = 10)
grid.arrange(grobs = temperature_plots, ncol = 3)  # Arrange the temperature plots in a grid
dev.off()

# Save combined precipitation plots into another PDF
pdf("~/R/cannabis_GEAV/Outputs/GEAV_combined_precipitation_plots.pdf", width = 8, height = 10)
grid.arrange(grobs = precipitation_plots, ncol = 3)  # Arrange the precipitation plots in a grid
dev.off()




############################
#Supplemental Figure 8
###########################
#REN + LW

#PHYLO 
#Figure XA and B (Temp and ppt)

#PHYLO GROUP
data2 <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/rrblup_GEBV_values_REN_LW.csv')
bio_vars <- sprintf("bio_%02d", 1:19)
library(gridExtra)
library(ggplot2)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type"))

# Define temperature and precipitation associations separately
temperature_associations <- c(
  "BIO 1 - Annual Mean Temperature",
  "BIO 2 - Mean Diurnal Range (Max temp - Min temp)",
  "BIO 3 - Isothermality (BIO2/BIO7) (×100)",
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 7 - Temperature Annual Range (BIO5-BIO6)",
  "BIO 8 - Mean Temperature of Wettest Quarter",
  "BIO 9 - Mean Temperature of Driest Quarter",
  "BIO 10 - Mean Temperature of Warmest Quarter",
  "BIO 11 - Mean Temperature of Coldest Quarter"
)

precipitation_associations <- c(
  "BIO 12 - Annual Precipitation",
  "BIO 13 - Precipitation of Wettest Month",
  "BIO 14 - Precipitation of Driest Month",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter",
  "BIO 18 - Precipitation of Warmest Quarter",
  "BIO 19 - Precipitation of Coldest Quarter"
)

# Create a list to hold the temperature and precipitation plots separately
temperature_plots <- list()
precipitation_plots <- list()

# Loop through each bio variable to generate the plots for temperature variables
for (i in 1:length(temperature_associations)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = temperature_associations[i], x = NULL, y = bio_vars[i]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue", "unknown" = "grey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to temperature plots
  temperature_plots[[i]] <- p
}

# Loop through each bio variable to generate the plots for precipitation variables
for (i in 1:length(precipitation_associations)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_vars[i + length(temperature_associations)], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = precipitation_associations[i], x = NULL, y = bio_vars[i + length(temperature_associations)]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to precipitation plots
  precipitation_plots[[i]] <- p
}

# Save combined temperature plots into one PDF
pdf("~/R/cannabis_GEAV/DataSetMerged/outputs/GEAV_combined_LW_REN_temperature_plots_phylogenetic_4x.pdf", width = 8, height = 10)
grid.arrange(grobs = temperature_plots, ncol = 3)  # Arrange the temperature plots in a grid
dev.off()

# Save combined precipitation plots into another PDF
pdf("~/R/cannabis_GEAV/DataSetMerged/outputs/GEAV_combinedLW_REN_precipitation_plots_phylogenetic_4x.pdf", width = 8, height = 10)
grid.arrange(grobs = precipitation_plots, ncol = 3)  # Arrange the precipitation plots in a grid
dev.off()


############################
#Supplemental Figure 9
###########################
#REN + LW

#CULT
#Figure 9A and B (Temp and ppt)
data2 <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/rrblup_GEBV_values_REN_LW.csv')

library(gridExtra)
library(ggplot2)
bio_vars <- sprintf("bio_%02d", 1:19)
# Reorder the factor levels for Cultivated_status
data2$Cultivated_status <- factor(data2$Cultivated_status, levels = c("feral", "Landrace", "cultivar"))

# Define temperature and precipitation associations separately
temperature_associations <- c(
  "BIO 1 - Annual Mean Temperature",
  "BIO 2 - Mean Diurnal Range (Max temp - Min temp)",
  "BIO 3 - Isothermality (BIO2/BIO7) (×100)",
  "BIO 4 - Temperature Seasonality (std dev ×100)",
  "BIO 5 - Max Temperature of Warmest Month",
  "BIO 6 - Min Temperature of Coldest Month",
  "BIO 7 - Temperature Annual Range (BIO5-BIO6)",
  "BIO 8 - Mean Temperature of Wettest Quarter",
  "BIO 9 - Mean Temperature of Driest Quarter",
  "BIO 10 - Mean Temperature of Warmest Quarter",
  "BIO 11 - Mean Temperature of Coldest Quarter"
)

precipitation_associations <- c(
  "BIO 12 - Annual Precipitation",
  "BIO 13 - Precipitation of Wettest Month",
  "BIO 14 - Precipitation of Driest Month",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation)",
  "BIO 16 - Precipitation of Wettest Quarter",
  "BIO 17 - Precipitation of Driest Quarter",
  "BIO 18 - Precipitation of Warmest Quarter",
  "BIO 19 - Precipitation of Coldest Quarter"
)

# Create a list to hold the temperature and precipitation plots separately
temperature_plots <- list()
precipitation_plots <- list()

# Define bio_vars to correspond to the column names for temperature and precipitation variables
bio_vars <- c("bio_01", "bio_02", "bio_03", "bio_04", "bio_05", "bio_06", "bio_07", "bio_08", "bio_09", "bio_10", "bio_11", 
              "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")


# Loop through each bio variable to generate the plots for temperature variables
for (i in 1:length(temperature_associations)) {
  p <- ggplot(data2, aes_string(x = "Cultivated_status", y = bio_vars[i], fill = "Cultivated_status")) +
    geom_boxplot(color = "black") +
    labs(title = temperature_associations[i], x = NULL, y = bio_vars[i]) +
    scale_fill_manual(values = c("feral" = "orange", "cultivar" = "green", "Landrace" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to temperature plots
  temperature_plots[[i]] <- p
}

# Loop through each bio variable to generate the plots for precipitation variables
for (i in 1:length(precipitation_associations)) {
  p <- ggplot(data2, aes_string(x = "Cultivated_status", y = bio_vars[i + length(temperature_associations)], fill = "Cultivated_status")) +
    geom_boxplot(color = "black") +
    labs(title = precipitation_associations[i], x = NULL, y = bio_vars[i + length(temperature_associations)]) +
    scale_fill_manual(values = c("feral" = "orange", "cultivar" = "green", "Landrace" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to precipitation plots
  precipitation_plots[[i]] <- p
}

# Save combined temperature plots into one PDF
pdf("~/R/cannabis_GEAV/DataSetMerged/outputs/GEAV_combined_LW_REN_temperature_plots.pdf", width = 8, height = 10)
grid.arrange(grobs = temperature_plots, ncol = 3)  # Arrange the temperature plots in a grid
dev.off()

# Save combined precipitation plots into another PDF
pdf("~/R/cannabis_GEAV/DataSetMerged/outputs/GEAV_combined_LW_REN_precipitation_plots.pdf", width = 8, height = 10)
grid.arrange(grobs = precipitation_plots, ncol = 3)  # Arrange the precipitation plots in a grid
dev.off()

############################
#Supplemental Figure 10
###########################

library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")

# Add a condition to label a certain percentage of points, or sample a subset for labeling
set.seed(123)  # For reproducibility
pca$label <- ifelse(runif(nrow(pca)) < 0.99, pca$Name, NA)  # Labeling 99 points

# Plot with sample names distributed across all points
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Phylogenetic_tree)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type feral" = "blue", 
                                "drug-type" = "red", 
                                "unknown" = "grey")) +  # Assign colors
  geom_text(aes(label = label), size = 3, vjust = -1, hjust = -0.2, check_overlap = TRUE) +  # Add sample names, ensuring reasonable overlap management
  labs(title = "Phylogenetic Tree Categorisation", 
       x = "PC1 (8.25%)", 
       y = "PC2 (4.46%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )

ggsave("~/R/cannabis_GEAV/MERGE/FIGS/PCA_S10.pdf", plot = p, device = "pdf", width = 30, height = 20)



############################
#Supplemental Figure 11
###########################
#Prediction on 8X - Ren 
#HPC

######################################################
#Step 4 
######################################################
# cross validation of models for prediction accuracy
library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

#load functions 
#source("~/kantar_koastore/anna/Barley_collab/barley_parental/cross_validation_2/xval_kfold_functions.R") 
source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/xval_kfold_functions.R") 
#source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R") 

#Read genotype data in rrBLUP format
data <- read.delim("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation_8X/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
colnames(data)

#Read environmental data from .csv
envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation_8X/ren_n44_extracted_climate_data_8X_n44.csv', head = T) # full environmental dataset

#Read training dataset from .csv
trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation_8X/all_8_through_time.csv", head = T) #only have env data for core in this


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


#########
#PLOTS
#########







############################
#Supplemental Figure 12
###########################
#GS on 8X + plots - Ren 

#Extract relevant environmental data from WorldClim 8X
library(raster)
worldclim_path <- "~/R/Data/climate_data_8X/"

# Define the file patterns based on the types of climate variables
prec_files <- paste0(worldclim_path, "wc2.1_30s_prec_", sprintf("%02d", 1:12), ".tif")  # Precipitation
srad_files <- paste0(worldclim_path, "wc2.1_30s_srad_", sprintf("%02d", 1:12), ".tif")  # Solar radiation
tavg_files <- paste0(worldclim_path, "wc2.1_30s_tavg_", sprintf("%02d", 1:12), ".tif")  # Mean temperature
tmax_files <- paste0(worldclim_path, "wc2.1_30s_tmax_", sprintf("%02d", 1:12), ".tif")  # Maximum temperature
tmin_files <- paste0(worldclim_path, "wc2.1_30s_tmin_", sprintf("%02d", 1:12), ".tif")  # Minimum temperature
elev_file <- paste0(worldclim_path, "wc2.1_30s_elev.tif")


elev_layer <- raster(elev_file)
# Stack
prec_layers <- stack(prec_files)
srad_layers <- stack(srad_files)
tavg_layers <- stack(tavg_files)
tmax_layers <- stack(tmax_files)
tmin_layers <- stack(tmin_files)

climate_layers <- stack(prec_layers, srad_layers, tavg_layers, tmax_layers, tmin_layers, elev_layer)


# Load your latitude and longitude dataset
data <- read.csv("~/R/cannabis_GEAV/Inputs/ren_lat_long.csv")

# Prepare the coordinates (assuming Latitude and Longitude columns are present)
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_8X.csv", row.names = FALSE)


######################################################
# Step 5B: Genomic Selection on 8X 
######################################################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/Inputs/Ren_PRJNA734114_rrBLUP_clean.txt", sep = "\t", header = TRUE)
#envdat <- read.csv('~/R/cannabis_GEAV/Outputs/ren_n44_extracted_climate_data_core_n25_greedy.csv', head = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Outputs/all_8_through_time.csv', head = TRUE)

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
write.csv(gebv_df, '~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_8X_2.csv')


####################
#Plots for 8X traits
####################
library(ggplot2)
library(tidyr)
library(dplyr)

# Load your data (assuming your file is saved as 'climate_data.csv')
data <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_8X_2.csv')

# Reshape the data into a long format
data_long <- pivot_longer(data, 
                          cols = starts_with("wc2.1_30s"), 
                          names_to = c("variable", "month"), 
                          names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",
                          values_to = "value")

# Convert month to numeric so that ggplot orders the months correctly
data_long$month <- as.numeric(data_long$month)

# Now plot the results by Phylogenetic_tree and over months
ggplot(data_long, aes(x = month, y = value, color = Phylogenetic_tree, group = Phylogenetic_tree)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free_y") +  # Create separate panels for each variable (e.g., prec, srad, tavg, etc.)
  labs(title = "Climate Variables Over Time by Phylogenetic Tree", 
       x = "Month", 
       y = "Value") +
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type" = "red", 
                                "drug-type feral" = "blue")) +
  theme_minimal()

###################################
# Now plot the results using boxplots by Phylogenetic_tree and over months
ggplot(data_long, aes(x = factor(month), y = value, fill = Phylogenetic_tree)) + 
  geom_boxplot() +  # Use boxplots instead of lines
  facet_wrap(~variable, scales = "free_y") +  # Create separate panels for each variable (e.g., prec, srad, tavg, etc.)
  labs(title = "Climate Variables Over Time by Phylogenetic Tree", 
       x = "Month", 
       y = "Value") + 
  scale_fill_manual(values = c("basal" = "orange", 
                               "hemp-type" = "green", 
                               "drug-type" = "red", 
                               "drug-type feral" = "blue")) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for readability

####################
#median plot line #PHYLO STATUS
####################
library(ggplot2)
library(dplyr)
data <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_8X_2.csv')

# Reshape the data into a long format
data_long <- pivot_longer(data, 
                          cols = starts_with("wc2.1_30s"), 
                          names_to = c("variable", "month"), 
                          names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",
                          values_to = "value")

# Summarize the data by month and variable
data_summary <- data_long %>%
  group_by(Phylogenetic_tree, variable, month) %>%
  summarize(median_value = median(value), .groups = "drop")

# Set the order of the variables
data_summary$variable <- factor(data_summary$variable, 
                                levels = c("tavg", "tmin", "tmax", "prec", "vapr", "wind", "srad"))

# Plot using the summarized data with the custom variable order
ggplot(data_summary, aes(x = factor(month), y = median_value, color = Phylogenetic_tree, group = Phylogenetic_tree)) + 
  geom_line(size = 1.5) +  # Line plot for median values
  facet_wrap(~variable, scales = "free_y", ncol = 3) +  # More columns for clarity
  labs(title = "Median Climate Variables Over Time by Phylogenetic Tree", 
       x = "Month", 
       y = "Median Value") + 
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type" = "red", 
                                "drug-type feral" = "blue")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


############################
#Supplemental Figure 13
###########################
#Prediction on 8X - Ren + LW
#HPC

######################################################
#Step 4 
######################################################
# cross validation of models for prediction accuracy
library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

#load functions 
source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/cross_validation/xval_kfold_functions.R") 

#Read genotype data in rrBLUP format
data <- read.delim("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_LW/cross_validation_8x/Ren_LW_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(data)[1] <- "rs." #needed to rename the sample column here
colnames(data) <- gsub("^X", "", colnames(data)) #needed to remove the extra X added in 

#Read environmental data from .csv
#envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_LW/cross_validation/ren_lw_n50_extracted_climate_data.csv', head = T) # full environmental dataset
envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_LW/cross_validation_8x/ren_lw_n55_extracted_climate_data_8X.csv', head = T) # full environmental dataset

#Read training dataset from .csv
#trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_LW/cross_validation/ren_lw_n50_extracted_climate_data.csv", head = T) #only have env data for core in this
trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_LW/cross_validation_8x/ren_lw_n55_extracted_climate_data_8X.csv", head = T) #only have env data for core in this


# Set the row names for the environmental and training datasets
row.names(envdat) <- envdat$Sample_ID
row.names(trainingset) <- trainingset$Sample_ID

# Filter genotype data to include only columns matching Sample_ID in envdat
colnames(data)
gd2 <- data[, c("rs.", colnames(data)[colnames(data) %in% envdat$Sample_ID])]

# Set the SNP IDs as the row names for the genotype data
row.names(gd2) <- gd2$rs.

# Remove the first four columns (rs., allele, chrom, pos) to create the genotype matrix
g.in <- gd2[, -c(1)]  # Isolate genotype data
g.in <- as.matrix(g.in)
g.in <- t(g.in)  # Transpose to align individuals as rows

### Remove SNPs with any missing data ###
#g.in <- g.in[, colSums(is.na(g.in)) == 0]  # Remove SNPs with missing data



# Impute missing values with the column mean, processing column by column
for (i in 1:ncol(g.in)) {
  g.in[is.na(g.in[, i]), i] <- mean(g.in[, i], na.rm = TRUE)
}

# Check dimensions after imputation
cat("Dimensions of g.in after imputation:", dim(g.in), "\n")


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

#RR-BLUP 
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


#########
#PLOTS
#########








############################
#Supplemental Figure 14
###########################
#GS on 8X + plots - Ren + LW

#Extract relevant environmental data from WorldClim 8X
library(raster)
worldclim_path <- "~/R/Data/climate_data_8X/"

# Define the file patterns based on the types of climate variables
prec_files <- paste0(worldclim_path, "wc2.1_30s_prec_", sprintf("%02d", 1:12), ".tif")  # Precipitation
srad_files <- paste0(worldclim_path, "wc2.1_30s_srad_", sprintf("%02d", 1:12), ".tif")  # Solar radiation
tavg_files <- paste0(worldclim_path, "wc2.1_30s_tavg_", sprintf("%02d", 1:12), ".tif")  # Mean temperature
tmax_files <- paste0(worldclim_path, "wc2.1_30s_tmax_", sprintf("%02d", 1:12), ".tif")  # Maximum temperature
tmin_files <- paste0(worldclim_path, "wc2.1_30s_tmin_", sprintf("%02d", 1:12), ".tif")  # Minimum temperature
elev_file <- paste0(worldclim_path, "wc2.1_30s_elev.tif")


elev_layer <- raster(elev_file)
# Stack
prec_layers <- stack(prec_files)
srad_layers <- stack(srad_files)
tavg_layers <- stack(tavg_files)
tmax_layers <- stack(tmax_files)
tmin_layers <- stack(tmin_files)

climate_layers <- stack(prec_layers, srad_layers, tavg_layers, tmax_layers, tmin_layers, elev_layer)


# Load your latitude and longitude dataset
#data <- read.csv("~/R/cannabis_GEAV/Inputs/ren_lat_long.csv")
data <- read.csv("~/R/cannabis_GEAV/DataSetMerged/inputs/ren_lw_lat_long.csv")


# Prepare the coordinates (assuming Latitude and Longitude columns are present)
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "~/R/cannabis_GEAV/DataSetMerged/outputs/ren_lw_n55_extracted_climate_data_8X.csv", row.names = FALSE)


##################
# Step 5: Genomic Selection on 8X Ren + LW
####################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/DataSetMerged/outputs/Ren_LW_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(gd1)[1] <- "rs."

envdat <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/ren_lw_n55_extracted_climate_data_8X.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

row.names(g.in) <- gsub("^X", "", row.names(g.in)) #accidental X added
print(g.in [1:10, 1:10])  # Visualize 

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$Sample_ID

# Subset environmental data for analysis (start from the 5th column)
y.in.rr <- envdat[, 5:ncol(envdat)]  
minicore_entries <- which(envdat$Core == TRUE)  # Subset to core lines where 'Core' is TRUE
y.in.rr <- y.in.rr[minicore_entries, ]  
y.in.mat <- as.matrix(y.in.rr)
print(y.in.mat [1:10, 1:10])  # Visualize 

dim(y.in.mat)

# Training set
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train, ]  # Subset rows (individuals) directly

print(g.train[1:10, 1:10])  # Visualize part of g.train for debugging
dim(g.train)


# Remove SNPs with any missing data in g.train
#g.train <- g.train[, colSums(is.na(g.train)) == 0]
#g.train <- g.train[, colSums(is.na(g.train)) == 0]
#dim(g.train)

# Check for remaining NA values in the genotype matrix
#if (any(is.na(g.train))) stop("There are still missing values in g.train after SNP removal.")

#cat("Dimensions of g.train after removing SNPs with missing data:", dim(g.train), "\n")

# Subset g.pred to include only individuals in `pred`
pred <- setdiff(row.names(g.in), train)  # This gives 57 individuals not in the training set
g.pred <- g.in[pred, ]


##################################################
# Calculate the overall percentage of missing values in the entire matrix
total_missing <- sum(is.na(g.train))
total_values <- prod(dim(g.train))
percentage_missing <- (total_missing / total_values) * 100

cat("Overall percentage of missing data:", percentage_missing, "%\n")
#10.22%
##################################################

common_snps <- intersect(colnames(g.train), colnames(g.pred))

# Subset both g.train and g.pred to keep only the common SNPs
g.train <- g.train[, common_snps]
g.pred <- g.pred[, common_snps]

# Impute missing values in g.train and g.pred
# Impute missing values with the column mean
g.train[is.na(g.train)] <- apply(g.train, 2, function(x) mean(x, na.rm = TRUE)) #warning is ok
g.pred[is.na(g.pred)] <- apply(g.pred, 2, function(x) mean(x, na.rm = TRUE)) #warning is ok

# Check dimensions after SNP alignment and imputation
cat("Dimensions of g.train after imputation:", dim(g.train), "\n")
cat("Dimensions of g.pred after imputation:", dim(g.pred), "\n")

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
row.names(gebv_df) <- row.names(g.in) #row names from your input
write.csv(gebv_df, '~/R/cannabis_GEAV/DataSetMerged/outputs/rrblup_GEBV_values_8X_REN_LW.csv')

####################
#Plots for 8X traits
####################
library(ggplot2)
library(tidyr)
library(dplyr)

# Load your data (assuming your file is saved as 'climate_data.csv')
data <- read.csv('~/R/cannabis_GEAV/DataSetMerged/outputs/rrblup_GEBV_values_8X_REN_LW.csv')

# Reshape the data into a long format
data_long <- pivot_longer(data, 
                          cols = starts_with("wc2.1_30s"), 
                          names_to = c("variable", "month"), 
                          names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",
                          values_to = "value")

# Summarize the data by month and variable
data_summary <- data_long %>%
  group_by(Phylogenetic_tree, variable, month) %>%
  summarize(median_value = median(value), .groups = "drop")

# Set the order of the variables
data_summary$variable <- factor(data_summary$variable, 
                                levels = c("tavg", "tmin", "tmax", "prec", "vapr", "wind", "srad"))

# Plot using the summarized data with the custom variable order
p<-ggplot(data_summary, aes(x = factor(month), y = median_value, color = Phylogenetic_tree, group = Phylogenetic_tree)) + 
  geom_line(size = 1.5) +  # Line plot for median values
  facet_wrap(~variable, scales = "free_y", ncol = 3) +  # More columns for clarity
  labs(title = "Median Climate Variables Over Time by Phylogenetic Tree", 
       x = "Month", 
       y = "Median Value") + 
  scale_color_manual(values = c("basal" = "orange", 
                                "hemp-type" = "green", 
                                "drug-type" = "red", 
                                "drug-type feral" = "blue")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("~/R/cannabis_GEAV/DataSetMerged/outputs/Climate_Median_Values_Plot.pdf", plot = p, device = "pdf", width = 10, height = 6)


####################################################################################
# FIN
###################################################################################




