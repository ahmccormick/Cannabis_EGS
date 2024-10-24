##########################
#Figure 1
##########################
vcf4 <- read.vcfR("~/R/cannabis_GEAV/MERGE2/merged_Ren_Soorni_snps.vcf.gz")

#SNPRelate
library(SNPRelate)
setwd("~/R/cannabis_GEAV/MERGE2/")
vcf.fn <- "~/R/cannabis_GEAV/MERGE2/merged_Ren_Soorni_snps.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "soorni_ren.gds", method="copy.num.of.ref")
snpgdsSummary("soorni_ren.gds")
genofile <- snpgdsOpen("soorni_ren.gds")
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
write.csv(tab, "~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")

###########
# PCA 
###########
library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")

#Figure 1A
# Phylogenetic_tree
p<-ggplot(pca, aes(x = EV1, y = EV2, color = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("Ren et al 2021" = "purple3", 
                                "Soorni et al 2017" = "green3")) +  # Assign colors
  labs(title = "Dataset Origins", 
       x = "PC1 (4.11%)", 
       y = "PC2 (2.92%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )

ggsave("~/R/cannabis_GEAV/FIGS/PCA_1A.pdf", plot = p, device = "pdf", width = 8, height = 6)

#Figure 1B
# Core
library(ggplot2)
pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Core_greedy)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses if needed
  scale_color_manual(values = c("TRUE" = "green3", "FALSE" = "grey")) +  # Color based on Core
  labs(title = "Core greedy ", 
       x = "PC1 (4.11%)", 
       y = "PC2 (2.92%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
p
ggsave("~/R/cannabis_GEAV/FIGS/PCA_1B.pdf", plot = p, device = "pdf", width = 8, height = 6)

#Figure 1C
# Phylogenetic_tree
pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")
p<-ggplot(pca, aes(x = EV1, y = EV2, color = Phylogenetic_tree)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("basal" = "orange", 
                                "drug-type feral" = "blue", 
                                "hemp-type" = "green",
                                "drug-type" = "red",
                                "Population_1" = "#9370DB",
                                "Population_2" = "#008080")) +  # Assign colors
  labs(title = "Phylogenetic Tree Categorisation", 
       x = "PC1 (4.11%)", 
       y = "PC2 (2.92%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )
p
ggsave("~/R/cannabis_GEAV/FIGS/PCA_1C.pdf", plot = p, device = "pdf", width = 8, height = 6)


#Type/Chemistry
pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")
ggplot(pca, aes(x = EV1, y = EV2, color = Type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("Type I" = "red",
                                "Type II" = "orange",
                                "Type III" = "green", 
                                "unknown" = "grey")) +  # Assign colors
  labs(title = "Phylogenetic Tree Categorisation", 
       x = "PC1 (4.11%)", 
       y = "PC2 (2.92%)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),        # Adjust legend title size
    legend.text = element_text(size = 10),         # Adjust legend text size
    axis.title.x = element_text(size = 14),        # Adjust x-axis label size
    axis.title.y = element_text(size = 14),        # Adjust y-axis label size
    axis.text = element_text(size = 12)            # Adjust x and y axis text size
  )


#######################################
#Figure 2
#######################################
#Koppen 
library(raster)
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Ren_Soorni_Lat_Long_n111.csv")
head(data)  # This should show Longitude, Latitude, and Phylogenetic_tree columns

# Load the Köppen-Geiger climate classification raster (assuming you have the file)
climate_raster <- raster("~/R/cannabis_GEAV/Beck_KG_V1/Beck_KG_V1_present_0p083.tif")

# Convert raster to a data frame
climate_df <- as.data.frame(climate_raster, xy = TRUE)
colnames(climate_df) <- c("Longitude", "Latitude", "Climate_Class")


climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc -Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter ",
  `21` = "Dwa- Cold, dry winter, hot summer",
  `22` = "Dwb- Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

# Filter the dataset to include only the 14 classes you're interested in
climate_df <- climate_df %>%
  filter(Climate_Class %in% names(climate_labels))

p <- ggplot() +
  geom_raster(data = climate_df, aes(x = Longitude, y = Latitude, fill = factor(Climate_Class))) +
  scale_fill_viridis_d(option = "viridis", labels = climate_labels, name = "Climate Class") +  # Use 'viridis' palette
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Phylogenetic_tree),
             size = 4, alpha = 0.9, shape = 21, stroke = 1.5, fill = "white") +
  labs(title = "Köppen-Geiger Climate Classes with Phylogenetic Data Overlay", 
       x = "Longitude", 
       y = "Latitude") +
  scale_color_manual(values = c("basal" = "orange", 
                                "drug-type feral" = "blue", 
                                "hemp-type" = "green", 
                                "Population_1" = "#9370DB",
                                "Population_2" = "#008080",
                                "Population_3"= "grey"), name = "Phylogenetic Tree") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  coord_cartesian(xlim = c(-10, 150), ylim = c(20, 80))
p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/outputs/koppen_5.pdf", plot = p, width = 25, height = 10)

#########
#Extract info from points for boxplot
#########
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Ren_Soorni_Lat_Long_n111.csv")
climate_raster <- raster("~/R/cannabis_GEAV/Beck_KG_V1/Beck_KG_V1_present_0p083.tif")

# Extract the climate class for each point 
coordinates <- as.matrix(data[, c("Longitude", "Latitude")])
data$Climate_Class <- raster::extract(climate_raster, coordinates)

climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc -Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter ",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

data$Climate_Class_Description <- climate_labels[as.character(data$Climate_Class)]
head(data)

write.csv(data, "~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_Lat_Long_n111_climate_class.csv", row.names = FALSE)

#########
#BOXPLOT
#########
library(ggplot2)
library(dplyr)
library(viridis)
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_Lat_Long_n111_climate_class.csv")

# Group the data by Phylogenetic_tree and Climate_Label to get the count
phylo_climate_summary <- data %>%
  group_by(Phylogenetic_tree, Climate_Class_Description) %>%
  summarise(Count = n()) %>%
  ungroup()

# Create the 31 colors from the viridis palette
climate_colors_full <- viridis(30, option = "viridis")

# Create the full mapping of climate classes to colors using the color gradient
climate_labels <- c(
  `1` = "Af - Tropical, rainforest",
  `2` = "Am - Tropical, monsoon",
  `3` = "Aw - Tropical, Savanna",
  `4` = "BWh - Arid, desert, hot",
  `5` = "BWk - Arid, desert, cold",
  `6` = "BSh - Arid, steppe, hot",
  `7` = "BSk - Arid, steppe, cold",
  `8` = "Csa - Temperate, dry summer, hot summer",
  `9` = "Csb - Temperate, dry summer, warm summer",
  `10` = "Csc -Temperate, dry summer, cold summer",
  `11` = "Cwa - Temperate, dry winter, hot summer",
  `12` = "Cwb - Temperate, dry winter, warm summer",
  `13` = "Cwc - Temperate, dry winter, cold summer",
  `14` = "Cfa - Temperate, no dry season, hot summer",
  `15` = "Cfb - Temperate, no dry season, warm summer",
  `16` = "Cfc - Temperate, no dry season, cold summer",
  `17` = "Dsa - Cold, dry summer, hot summer",
  `18` = "Dsb - Cold, dry summer, warm summer",
  `19` = "Dsc - Cold, dry summer, cold summer",
  `20` = "Dsd - Cold, dry summer, very cold winter ",
  `21` = "Dwa - Cold, dry winter, hot summer",
  `22` = "Dwb - Cold, dry winter, warm summer",
  `23` = "Dwc - Cold, dry winter, cold summer",
  `24` = "Dwd - Cold, dry winter, very cold winter",
  `25` = "Dfa - Cold, no dry season, hot summer",
  `26` = "Dfb - Cold, no dry season, warm summer",
  `27` = "Dfc - Cold, no dry season, cold summer",
  `28` = "Dfd - Cold, no dry season, very cold winter",
  `29` = "ET - Polar, tundra",
  `30` = "EF - Polar, frost"
)

# Combine climate labels with the color gradient
climate_colors <- setNames(climate_colors_full, climate_labels)

# Display the full mapping of climate classes to colors
climate_colors
#then match

climate_colors <- c(
  "Af - Tropical, rainforest" = "#440154FF",
  "Am - Tropical, monsoon" = "#470E61FF",  # Updated
  "Aw - Tropical, Savanna" = "#481B6DFF",  # Updated
  "BWh - Arid, desert, hot" = "#482576FF", 
  "BWk - Arid, desert, cold" = "#46307EFF",  # Updated
  "BSh - Arid, steppe, hot" = "#443B84FF",  # Updated
  "BSk - Arid, steppe, cold" = "#404688FF",  # Updated
  "Csa - Temperate, dry summer, hot summer" = "#3C508BFF",  # Updated
  "Csb - Temperate, dry summer, warm summer" = "#38598CFF",  # Updated
  "Csc - Temperate, dry summer, cold summer" = "#33628DFF",  # Updated
  "Cwa - Temperate, dry winter, hot summer" = "#2F6B8EFF",  # Updated
  "Cwb - Temperate, dry winter, warm summer" = "#2C738EFF",  # Updated
  "Cwc - Temperate, dry winter, cold summer" = "#287C8EFF",  # Updated
  "Cfa - Temperate, no dry season, hot summer" = "#25838EFF",  # Updated
  "Cfb - Temperate, no dry season, warm summer" = "#228C8DFF",  # Updated
  "Cfc - Temperate, no dry season, cold summer" = "#1F948CFF",  # Updated
  "Dsa - Cold, dry summer, hot summer" = "#1E9D89FF",  # Updated
  "Dsb - Cold, dry summer, warm summer" = "#20A486FF",  # Updated
  "Dsc - Cold, dry summer, cold summer" = "#26AD81FF",  # Updated
  "Dsd - Cold, dry summer, very cold winter" = "#31B57BFF",  # Updated
  "Dwa - Cold, dry winter, hot summer" = "#3FBC73FF",  # Updated
  "Dwb - Cold, dry winter, warm summer" = "#4FC46AFF",  # Updated
  "Dwc - Cold, dry winter, cold summer" = "#61CB5FFF",  # Updated
  "Dwd - Cold, dry winter, very cold winter" = "#75D054FF",  # Updated
  "Dfa - Cold, no dry season, hot summer" = "#8BD646FF",  # Updated
  "Dfb - Cold, no dry season, warm summer" = "#A2DA37FF",  # Updated
  "Dfc - Cold, no dry season, cold summer" = "#B9DE28FF",  # Updated
  "Dfd - Cold, no dry season, very cold winter" = "#D1E11CFF",  # Updated
  "ET - Polar, tundra" = "#E8E419FF",  # Updated
  "EF - Polar, frost" = "#FDE725FF"
)


p <- ggplot(phylo_climate_summary, aes(x = Phylogenetic_tree, y = Count, fill = Climate_Class_Description)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = Climate_Class_Description), 
            position = position_fill(vjust = 0.5), 
            size = 1.8, color = "white") +  # White text for visibility
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Tendency of Samples Towards Climate Classes by Phylogenetic Tree Grouping",
       x = "Phylogenetic Tree Grouping",
       y = "Proportion of Samples",
       fill = "Climate Class") +
  scale_fill_manual(values = climate_colors) +  # Apply custom color palette
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
p
# Save the plot
ggsave("~/R/cannabis_GEAV/Ren_Soorni/outputs/boxplot_output_FIN.pdf", plot = p, width = 16, height = 8, units = "in")





#######################################
#Figure 3
#######################################
#Figure 3A 
#Ren GEAVs plotted for 6X env form WorldClim
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
#Soorni GEAVs plotted for 6X env form WorldClim
#PHYLO GROUP
data2 <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy.csv')

library(gridExtra)
library(ggplot2)

bio_vars <- sprintf("bio_%02d", 1:19)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type", "Population_1", "Population_2"))


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
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue", "Population_1" = "#9370DB","Population_2" = "#008080")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 8)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add the plot to the list
  selected_plots[[i]] <- p
}

# Save the selected plots into a PDF
pdf("~/R/cannabis_GEAV/Ren_Soorni/outputs/GEAV_selected_bio_plots_Ren_Soorni_n50_core_greedy.pdf", width = 8, height = 10)
grid.arrange(grobs = selected_plots, ncol = 2)  # Arrange the plots in a grid
dev.off()



#######################################
#Figure 4
#######################################
#Conceptual figure 


###################################################################################################################################################
#SUPPLEMENTAL FIGURES
#########################
#Supplemental Figure 1
##########################
########
#MAP
########
# Load necessary libraries
library(ggplot2)
library(maps)

# Load the world map data
world_map <- map_data("world")

# Load your data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_Lat_Long_n111_climate_class.csv")

# Create the plot with the world map, limit longitude and latitude, and color by Phylogenetic_tree
p<-ggplot() +
  # Plot the world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  # Add points for your latitude and longitude data, colored by Phylogenetic_tree
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Phylogenetic_tree), 
             size = 3, alpha = 0.7) +
  labs(title = "Geographical Distribution of Phylogenetic Categories", 
       x = "Longitude", 
       y = "Latitude") + 
  scale_color_manual(values = c("basal" = "orange", 
                                "drug-type feral" = "blue", 
                                "hemp-type" = "green", 
                                "Population_1" = "#9370DB",
                                "Population_2" = "#008080")) +  # Customize colors
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  # Limit the longitude and latitude range to focus on the region you want
  coord_cartesian(xlim = c(-10, 150), ylim = c(10, 80))
p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/outputs/tree_ranges.pdf", plot = p, width = 12, height = 5)

#####
#Boxplots of data ranges for the groups
#####
library(dplyr)
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_GEBV.csv")

# List of WorldClim variables
worldclim_vars <- c("bio_01", "bio_02", "bio_03", "bio_04", "bio_05", "bio_06", "bio_07", 
                    "bio_08", "bio_09", "bio_10", "bio_11", "bio_12", "bio_13", 
                    "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")

# Calculate range for each WorldClim variable grouped by 'Phylogenetic_tree'
env_range <- data %>%
  group_by(Phylogenetic_tree) %>%
  summarise(across(all_of(worldclim_vars), ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE), .names = "range_{col}"))

# View the result
print(env_range)

####
# Updated lookup table with descriptive labels and units for WorldClim variables
worldclim_labels <- c(
  bio_01 = "Annual Mean Temperature (°C)",
  bio_02 = "Mean Diurnal Range",
  bio_03 = "Isothermality",
  bio_04 = "Temperature Seasonality (Standard Deviation)",
  bio_05 = "Max Temperature of Warmest Month (°C)",
  bio_06 = "Min Temperature of Coldest Month (°C)",
  bio_07 = "Temperature Annual Range",
  bio_08 = "Mean Temperature of Wettest Quarter (°C)",
  bio_09 = "Mean Temperature of Driest Quarter (°C)",
  bio_10 = "Mean Temperature of Warmest Quarter (°C)",
  bio_11 = "Mean Temperature of Coldest Quarter (°C)",
  bio_12 = "Annual Precipitation (mm)",
  bio_13 = "Precipitation of Wettest Month (mm)",
  bio_14 = "Precipitation of Driest Month (mm)",
  bio_15 = "Precipitation Seasonality (Coefficient of Variation)",
  bio_16 = "Precipitation of Wettest Quarter (mm)",
  bio_17 = "Precipitation of Driest Quarter (mm)",
  bio_18 = "Precipitation of Warmest Quarter (mm)",
  bio_19 = "Precipitation of Coldest Quarter (mm)"
)

# Faceted boxplot with custom labels, colors, white background, grey gridlines, and units in titles
p<- ggplot(data_long, aes(x = Phylogenetic_tree, y = Value, fill = Phylogenetic_tree)) +
  geom_boxplot() +
  facet_wrap(~ WorldClim_variable, scales = "free_y", 
             labeller = as_labeller(worldclim_labels)) +
  scale_fill_manual(values = color_mapping) +  # Apply the custom color mapping
  theme_minimal() +  # Minimal theme with white background
  theme(
    panel.grid.major = element_line(color = "lightgrey"),  # Light grey major gridlines
    panel.grid.minor = element_line(color = "lightgrey"),  # Light grey minor gridlines
    panel.background = element_rect(fill = "white"),       # White background
    axis.text.x = element_text(angle = 90, hjust = 1)      # Rotate x-axis text
  ) +
  labs(title = "Boxplots of WorldClim Variables by Phylogenetic Tree",
       x = "Phylogenetic Tree",
       y = "Value",
       fill = "Phylogenetic Tree")

ggsave("~/R/cannabis_GEAV/Ren_Soorni/outputs/FS1B_worldclim_boxplots.pdf", plot = p, width = 20, height = 10)


#########################
#Supplemental Figure 2
##########################
#Ren alone PCAs


#done


#########################
#Supplemental Figure 3
##########################
#Phylogeny






#########################
#Supplemental Figure 4
##########################
#fastSTRUCTURE


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
test <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/canna_k2.csv")

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
                 "Population_1" = "#9370DB",
                 "Population_2" = "#008080")

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
ggsave(filename = "~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/fastSTRUCTURE_plot_k2.pdf", 
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
test <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/canna_k3.csv")

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
                 "Population_1" = "#9370DB",
                 "Population_2" = "#008080")

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
ggsave(filename = "~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/fastSTRUCTURE_plot_k3.pdf", 
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

test <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/canna_k4.csv")

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
                 "Population_1" = "#9370DB",
                 "Population_2" = "#008080")

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
ggsave(filename = "~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/fastSTRUCTURE_plot_k4.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)



##############
#K5
##############
hiroshige_palette <- met.brewer("Hiroshige", 8)

test <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/canna_k5.csv")

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
                 "Population_1" = "#9370DB",
                 "Population_2" = "#008080")

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
ggsave(filename = "~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/fastSTRUCTURE_plot_k5.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)


##############
#K6
##############
#load fastSTRUCTURE data

# Load your data
test <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/canna_k6.csv")

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
                 "Population_1" = "#9370DB",
                 "Population_2" = "#008080")

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
ggsave(filename = "~/R/cannabis_GEAV/Ren_Soorni/fastSTRUCTURE/fastSTRUCTURE_plot_k6.pdf", 
       plot = final_plot, 
       width = 8,    # Width of the plot (adjust as needed)
       height = 3)   # Reduced height (adjust as needed)





#########################################
#SUPPLEMENTARY FIGURE 4B
#########################################
library("FactoMineR")
library("factoextra")

#pca <- read.csv("~/R/cannabis_GEAV/MERGE/lw_ren_merged_PC.csv")
pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_ren_merged_PC.csv")
#Select 3rd and 4th principal components
#pca2 <- pca[3:4]
pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")
a
ggsave("~/R/cannabis_GEAV/MERGE/Soorni_Ren_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")
b
ggsave("~/R/cannabis_GEAV/MERGE/Soorni_Ren_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")


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



#########################
#Supplemental Figure 5
##########################
#Cross validation (WorldClim 19) for Ren dataset


#done





#########################
#Supplemental Figure 6
##########################
#Cross validation (WorldClim 19) for Ren + Soorni dataset







#########################
#Supplemental Figure 7
##########################
#Predictions for (WorldClim 19) for Ren (phylogenetic groupings)


#done









#########################
#Supplemental Figure 8
##########################
#Predictions for (WorldClim 19) for Ren (cultivated status groupings)


#done





#########################
#Supplemental Figure 9
##########################
#Predictions for (WorldClim 19) for Ren + Soorni (PCA groupings)







#########################
#Supplemental Figure 10
##########################
#Cross validation (WorldClim 8X) for Ren dataset

#done



#########################
#Supplemental Figure 11
##########################
#Predictions for (WorldClim 8X) for Ren (phylogenetic groupings)

#done





#########################
#Supplemental Figure 12
##########################
#Cross validation (WorldClim 8X) for Ren + Soorni dataset








#########################
#Supplemental Figure 13
##########################
#Predictions for (WorldClim 8X) for Ren + soorni datasets (PCA groupings)







