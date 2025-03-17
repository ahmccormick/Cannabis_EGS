#############################################
#Figure 1
#############################################
#Figure 1A
################
#MAP
########
library(ggplot2)
library(maps)
library(dplyr)

# Load the world map data
world_map <- map_data("world")

# Load your data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_Lat_Long_n111_climate_class.csv")

# Replace underscores with dashes in Phylogenetic_tree and capitalize first letters
data$Phylogenetic_tree <- gsub("_", "-", data$Phylogenetic_tree)
data$Phylogenetic_tree <- tools::toTitleCase(data$Phylogenetic_tree)  # Capitalize first letters

# Count the number of samples in each Phylogenetic_tree category
counts <- data %>% 
  count(Phylogenetic_tree) %>% 
  mutate(label = paste0(Phylogenetic_tree, " (n=", n, ")"))

# Create a named vector of colors for each category
colors <- c("Basal" = "orange", 
            "Drug-Type Feral" = "blue", 
            "Hemp-Type" = "green", 
            "Population-1" = "#9370DB",
            "Population-2" = "#008080")

# Create the plot with the world map, limit longitude and latitude, and color by Phylogenetic_tree
p <- ggplot() +
  # Plot the world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  # Add points for your latitude and longitude data, colored by Phylogenetic_tree
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Phylogenetic_tree), 
             size = 3, alpha = 0.7) +
  labs(title = "", 
       x = "Longitude", 
       y = "Latitude") + 
  # Apply custom colors and labels with counts in the legend
  scale_color_manual(
    name = "Group",  # Change legend title to "Group"
    values = colors,
    labels = setNames(counts$label, counts$Phylogenetic_tree)  # Use updated labels
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14)    # Increase legend text size
  ) +
  # Limit the longitude and latitude range to focus on the region you want
  coord_cartesian(xlim = c(-10, 150), ylim = c(10, 80))

# Display the plot
print(p)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/F1Amap_plot_updated.pdf",  
  plot = p,
  device = "pdf",
  width = 15,
  height = 8,
  units = "in"
)


################
#Figure 1B
################
library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)

# Load data
data2 <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy.csv')

# Selected variables and labels
selected_bio_vars <- c("bio_01", "bio_09", "bio_14", "bio_15")
selected_labels <- c(
  "BIO 1 - Annual Mean Temperature (PA=0.57)",
  "BIO 9 - Mean Temperature of the Driest Quarter (PA=0.77)",
  "BIO 14 - Precipitation of Driest Month (PA=0.71)",
  "BIO 15 - Precipitation Seasonality (Coefficient of Variation) (PA=0.71)"
)

# Capitalize the first letter of each Phylogenetic_tree group
data2$Phylogenetic_tree <- gsub("_", "-", data2$Phylogenetic_tree) # Replace underscores with dashes
data2$Phylogenetic_tree <- tools::toTitleCase(data2$Phylogenetic_tree) # Capitalize first letter of each word

# Calculate sample counts for each Phylogenetic_tree category
counts <- data2 %>%
  count(Phylogenetic_tree) %>%
  mutate(label = paste0(Phylogenetic_tree, " (n=", n, ")"))  # Format labels with counts

# Define colors for Phylogenetic_tree groups
colors <- c("Basal" = "orange", "Hemp-Type" = "green", "Drug-Type" = "red", 
            "Drug-Type Feral" = "blue", "Population-1" = "#9370DB", "Population-2" = "#008080")

# Reorder factor levels
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("Basal", "Hemp-Type", "Drug-Type Feral", "Drug-Type", "Population-1", "Population-2"))

# Create individual plots without legends
selected_plots <- list()
for (i in seq_along(selected_bio_vars)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = selected_bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = selected_labels[i], x = NULL, y = "GEAV") +  # Change y-axis label to GEAV
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      legend.position = "none",  # Hide legend for individual plots
      plot.title = element_text(size = 10),  # Reduce plot title size
      axis.text.y = element_text(size = 12),  # Increase y-axis tick label size
      axis.title.y = element_text(size = 14)  # Increase y-axis label size
    ) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  selected_plots[[i]] <- p
}

# Extract legend with updated labels showing counts, replacing _ with -, and capitalizing first letters
legend <- get_legend(
  ggplot(data2, aes(x = Phylogenetic_tree, y = bio_01, fill = Phylogenetic_tree)) +
    geom_boxplot() +
    scale_fill_manual(
      name = "Group",  # Change the legend title to "Group"
      values = colors,
      labels = setNames(counts$label, counts$Phylogenetic_tree)  # Use updated labels with counts
    ) +
    theme(legend.position = "right")
)

# Arrange plots in a single row with a single legend on the right
combined_plot <- plot_grid(
  plot_grid(plotlist = selected_plots, ncol = 4, align = 'v'),
  legend,
  ncol = 2,
  rel_widths = c(4, 0.5)
)

# Display the combined plot
print(combined_plot)


ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/F1B.pdf",  
  plot = combined_plot,                                                      
  device = "pdf",                                                
  width = 20,                                                   
  height = 6,                                                    
  units = "in"                                                   
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


# Group the data by Phylogenetic_tree and Climate_Class_Description to get the count
phylo_climate_summary <- data %>%
  group_by(Phylogenetic_tree, Climate_Class_Description) %>%
  summarise(Count = n(), .groups = "drop")

# Adjust x-axis labels: Capitalize first letter and replace _ with -
phylo_climate_summary <- phylo_climate_summary %>%
  mutate(Phylogenetic_tree = gsub("_", "-", tools::toTitleCase(Phylogenetic_tree)))

# Define climate colors

# Create the plot
p <- ggplot(phylo_climate_summary, aes(x = Phylogenetic_tree, y = Count, fill = Climate_Class_Description)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = Climate_Class_Description),
            position = position_fill(vjust = 0.5),
            size = 1.8, color = "white") +  # White text for visibility
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "",
    x = "",
    y = "Proportion of Samples",
    fill = "Climate Class"
  ) +
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

# Display the plot
print(p)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure2B.pdf",
  plot = p,
  device = "pdf",
  width = 16,
  height = 8,
  units = "in"
)


###########################
#Figure 3
###########################
#Figure 3A
#############
#groups with no other population level admixture
#############
library(ggplot2)

# Load data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups_5plus.csv", row.names = 1)

# Capitalize the first letter of each Phylogenetic_tree group and replace underscores with dashes
data$Phylogenetic_tree <- gsub("_", "-", data$Phylogenetic_tree)
data$Phylogenetic_tree <- tools::toTitleCase(data$Phylogenetic_tree)

# Define colors with updated group names
colors <- c(
  "Basal" = "orange",
  "Drug-Type Feral" = "blue",
  "Hemp-Type" = "green",
  "Population-1" = "#9370DB",
  "Population-2" = "#008080"
)

# Create the plot
plot <- ggplot(data, aes(x = pop, fill = Phylogenetic_tree)) +
  geom_bar() +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  theme_minimal() +
  labs(
    title = "",
    x = "Climate Niche",
    y = "Number of Samples",
    fill = "Group"  # Change legend title to "Group"
  ) +
  scale_fill_manual(
    name = "Group",  # Update legend title
    values = colors  # Apply updated colors
  ) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12),  # Adjust y-axis label size
    axis.text.x = element_text(size = 12),                        # Adjust x-axis label size
    axis.title.x = element_text(size = 14),                       # Adjust x-axis title size
    axis.title.y = element_text(size = 14),                       # Adjust y-axis title size
    plot.title = element_text(size = 16, face = "bold"),          # Adjust title size
    legend.title = element_text(size = 14),                       # Adjust legend title size
    legend.text = element_text(size = 12)                         # Adjust legend text size
  )

# Display the plot
print(plot)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure_3A.pdf",  # Specify the file path and name
  plot = plot,
  device = "pdf",
  width = 10,  # Set the width in inches
  height = 7,  # Set the height in inches
  units = "in"  # Specify the units (inches)
)


########################
#Figure 3B
########################
library(ggplot2)
library(reshape2)
library(dplyr)

# Define the custom color palette
climate_colors <- c(
  "BWh - Arid, desert, hot (n=9)" = "#482576FF", 
  "BWk - Arid, desert, cold (n=7)" = "#46307EFF",  
  "BSh - Arid, steppe, hot (n=6)" = "#443B84FF",  
  "BSk - Arid, steppe, cold (n=20)" = "#404688FF",  
  "Csa - Temperate, dry summer, hot summer (n=5)" = "#3C508BFF",  
  "Cwa - Temperate, dry winter, hot summer (n=6)" = "#2F6B8EFF",  
  "Dsa - Cold, dry summer, hot summer (n=5)" = "#1E9D89FF",  
  "Dfb - Cold, no dry season, warm summer (n=10)" = "#A2DA37FF"
)

# Read the data
data <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy.csv')

# List of sample IDs to keep
sample_ids <- c("SRR9089275.1", "SRR9089276.1", "SRR9089281.1", "SRR9089282.1", "SRR9089285.1", "SRR9089286.1", 
                "SRR9089230.1", "SRR9089255.1", "SRR9089256.1", "SRR9089257.1", "SRR9089258.1", "SRR9089263.1", 
                "SRR9089264.1", "SRR9089269.1", "SRR9089270.1", "SRR9089277.1", "SRR9089278.1", "SRR9089283.1", 
                "SRR9089284.1", "SRR9089289.1", "SRR9089290.1", "SRR9089296.1", "SRR9089301.1", "SRR9089308.1", 
                "SRR9089309.1", "SRR9089310.1", "SRR9089225.1", "SRR9089267.1", "SRR9089268.1", "SRR9089287.1", 
                "SRR9089288.1", "SRR9089299.1", "SRR9089303.1", "SRR9089304.1", "SRR9089317.1", "SRR14708203", 
                "SRR14708204", "SRR14708205", "SRR14708206", "SRR14708207", "SRR14708208", "SRR14708209", 
                "SRR14708199", "SRR14708198", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258", 
                "SRR14708213", "SRR14708217", "SRR14708230", "SRR14708231", "SRR14708232", "SRR14708244", 
                "SRR14708259", "SRR14708273", "SRR14708275", "SRR14708278", "SRR9089228.1", "SRR9089229.1", 
                "SRR9089291.1", "SRR9089305.1", "SRR9089306.1")

# Filter the data to keep only specified sample IDs
filtered_data <- data %>% filter(Sample_ID %in% sample_ids)

# Further filter down to specific climate classes
filtered_data <- filtered_data %>%
  mutate(Climate_Class_Description = recode(Climate_Class_Description,
                                            "BWh - Arid, desert, hot" = "BWh - Arid, desert, hot (n=9)",
                                            "BWk - Arid, desert, cold" = "BWk - Arid, desert, cold (n=7)",
                                            "BSh - Arid, steppe, hot" = "BSh - Arid, steppe, hot (n=6)",
                                            "BSk - Arid, steppe, cold" = "BSk - Arid, steppe, cold (n=20)",
                                            "Csa - Temperate, dry summer, hot summer" = "Csa - Temperate, dry summer, hot summer (n=5)",
                                            "Cwa - Temperate, dry winter, hot summer" = "Cwa - Temperate, dry winter, hot summer (n=6)",
                                            "Dsa - Cold, dry summer, hot summer" = "Dsa - Cold, dry summer, hot summer (n=5)",
                                            "Dfb - Cold, no dry season, warm summer" = "Dfb - Cold, no dry season, warm summer (n=10)"
  ))

# Reshape the data from wide to long format
data_long <- melt(filtered_data, id.vars = c("Sample_ID", "Climate_Class_Description"),
                  measure.vars = paste0("bio_", sprintf("%02d", 1:19)),
                  variable.name = "Bio_Trait", value.name = "GEBV")

# Separate the precipitation variables (bio_12 to bio_19)
precipitation_data <- subset(data_long, Bio_Trait %in% paste0("bio_", sprintf("%02d", 12:19)))

# Filter for specific precipitation traits: bio_12, bio_16, and bio_18
selected_precipitation_data <- precipitation_data %>%
  filter(Bio_Trait %in% c("bio_12", "bio_16", "bio_18"))

p_selected <- ggplot(selected_precipitation_data, aes(x = Bio_Trait, y = GEBV, fill = Climate_Class_Description)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),  # Increase x-axis label text size
    axis.text.y = element_text(size = 22),                         # Increase y-axis label text size
    axis.title.x = element_text(size = 15),                        # Increase x-axis title text size
    axis.title.y = element_text(size = 23),                        # Increase y-axis title text size
    plot.title = element_text(size = 14, face = "bold"),           # Increase plot title text size
    legend.title = element_text(size = 20),                        # Increase legend title text size
    legend.text = element_text(size = 20)                          # Increase legend text size
  ) +
  labs(
    title = "",  # Plot title
    x = "",      # x-axis title
    y = "GEAV"   # y-axis title
  ) +
  scale_x_discrete(labels = c(
    "bio_12" = "BIO 12 - Annual Precipitation (PA=0.56)",
    "bio_16" = "BIO 16 - Precipitation of Wettest Quarter (PA=0.73)",
    "bio_18" = "BIO 18 - Precipitation of Warmest Quarter (PA=0.77)"
  )) +  # Custom x-axis labels
  scale_fill_manual(values = climate_colors) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Climate Class Description"))

# Display the plot
print(p_selected)

ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Precipitation_GEBV_Boxplot.pdf",  
  plot = p_selected,                                                            
  device = "pdf",                                                               
  width = 30,                                                                   
  height = 12,                                                                   
  units = "in"                                                                  
)


########################
#FIGURE 3C
########################
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the filtered data
filtered_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_8X_climate_niches_to_5above.csv")

# Reshape the data into a long format
seasonal_long <- filtered_data %>%
  pivot_longer(
    cols = starts_with("wc2.1_30s"),       # Select columns that start with "wc2.1_30s" for monthly climate variables
    names_to = c("variable", "month"),     # Split into variable (e.g., "prec") and month (e.g., "01")
    names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",  # Regex pattern to capture variable and month
    values_to = "value"
  ) %>%
  mutate(
    month = as.numeric(month)              # Convert month to numeric for ordering
  )

# Load overlap data to get the count per climate class
overlap_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups_5plus.csv", row.names = 1)


table(is.na(overlap_data$pop))  # Check for missing values
unique(overlap_data$pop)       # Inspect unique values

# Perform the count
climate_counts <- overlap_data %>%
  count(pop)

# Rename the columns manually
colnames(climate_counts) <- c("climate_class", "count")

# Check the result
print(climate_counts)

# Extract counts from the overlap data using the "pop" column for climate class descriptions
climate_counts <- overlap_data %>%
  count(pop) %>%
  rename(climate_class = pop, count = n)


# Calculate median and quartiles by Climate_Class_Description and add climate class counts
summary_data <- seasonal_long %>%
  group_by(variable, month, Climate_Class_Description) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  ) %>%
  left_join(climate_counts, by = c("Climate_Class_Description" = "climate_class"))

# Update labels with counts
summary_data <- summary_data %>%
  mutate(label = paste0(Climate_Class_Description, " (n=", count, ")"))

# Filter data for only the precipitation variable
prec_data <- summary_data %>%
  filter(variable == "prec")

# Define custom colors for each Climate_Class_Description with counts
custom_colors <- c(
  "BSh - Arid, steppe, hot (n=6)" = "purple3",
  "BSk - Arid, steppe, cold (n=20)" = "purple",
  "BWh - Arid, desert, hot (n=9)" = "purple4",
  "BWk - Arid, desert, cold (n=7)" = "orange",
  "Cwa - Temperate, dry winter, hot summer (n=6)" = "blue",
  "Dfb - Cold, no dry season, warm summer (n=10)" = "green",
  "Dsa - Cold, dry summer, hot summer (n=5)" = "#1ABC9C"
)

# Plot seasonal precipitation data with median lines and IQR ribbons for each climate class
# Define the plot
precipitation_plot <- ggplot(prec_data, aes(x = month, y = median_value, color = label, fill = label)) +
  geom_line(size = 1) +  # Median line
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), alpha = 0.2) +  # Shaded area for IQR
  labs(
    title = "",
    x = "Month",
    y = "GEAV"
  ) + 
  scale_color_manual(values = custom_colors) +  # Custom colors for lines
  scale_fill_manual(values = custom_colors) +   # Matching fill color for ribbon
  scale_x_continuous(breaks = 1:12, labels = 1:12) +  # Ensure x-axis shows months 1 to 12
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),         # Increase legend text size
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank(),
    text = element_text(size = 14)  # Increase general text size for readability
  )

# Display the plot
print(precipitation_plot)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure_3C.pdf",  # Specify the file path and name
  plot = precipitation_plot,                                                      # Use the plot object
  device = "pdf",                                                                 # Specify the file format
  width = 15,                                                                     # Set the width in inches
  height = 8,                                                                     # Set the height in inches
  units = "in"                                                                    # Specify the units (inches)
)


#######
#Figure 4A
######
# Filter data for the Cwa climate class
cwa_data <- prec_data %>%
  filter(label == "Cwa - Temperate, dry winter, hot summer (n=6)")

# Plot seasonal precipitation data for Cwa climate class
cwa_plot <- ggplot(cwa_data, aes(x = month, y = median_value, color = label, fill = label)) +
  geom_line(size = 1) +  # Median line
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), alpha = 0.2) +  # Shaded area for IQR
  labs(
    title = "",
    x = "Month",
    y = "GEAV"
  ) + 
  scale_color_manual(values = custom_colors) +  # Custom color for line
  scale_fill_manual(values = custom_colors) +   # Matching fill color for ribbon
  scale_x_continuous(breaks = 1:12, labels = 1:12) +  # Ensure x-axis shows months 1 to 12
  theme_minimal() +
  theme(
    legend.position = "none",  # Hide legend as only one class is plotted
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(cwa_plot)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure_4A_Cwa.pdf",  # Specify the file path and name
  plot = cwa_plot,                                                  # Use the plot object
  device = "pdf",                                                   # Specify the file format
  width = 10,                                                       # Set the width in inches
  height = 6,                                                       # Set the height in inches
  units = "in"                                                      # Specify the units (inches)
)

###########################
#Figure 4
###########################
# conversion of vcf to genotype matrix
# filtering for no missingness 
# per snp and individual marker effect size calculation
# significantly changing SNPs between individuals in the two time points

############attempt 5
############### FROM VCF TO GENOTYPE MATRIX - rows=individuals  columns =SNPs, in 0/1/2 format
library(vcfR)
library(adegenet)
vcf <- read.vcfR("~/R/cannabis_GEAV/MERGE2/merged_Ren_Soorni_snps.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")


# Function to convert genotypes to numeric format
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, ifelse(x == "0/1" | x == "1/0", 1, ifelse(x == "1/1", 2, NA)))
}

# Apply this function to each column of the genotype matrix
geno_numeric <- apply(genotype_matrix, 2, convert_geno)


geno_numeric_t <- t(geno_numeric)
geno_df <- as.data.frame(geno_numeric_t)


rownames(geno_df) <- colnames(genotype_matrix) 

snp_names <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")
colnames(geno_df) <- snp_names

# Write the genotype matrix to a CSV file
write.csv(geno_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv") #geno type matrix

############################ 
#filter genotype matrix for no missingness
############################ 
library(vegan)
library(adegenet)
library(psych)
library(rrBLUP)
library(dplyr)
library(tidyr)

genotype_matrix <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv')  #rows = individuals and columns = SNPs, in 0/1/2 format

# Count the number of columns with no missing values
no_missing_count <- sum(colSums(is.na(genotype_matrix)) == 0)
cat("Number of SNPs with no missing values:", no_missing_count, "\n")

# Identify columns with no missing values
complete_markers <- genotype_matrix[, colSums(is.na(genotype_matrix)) == 0]
# Confirm the new dimensions and number of markers
dim(complete_markers)

# Save complete_markers as a CSV file
write.csv(complete_markers, "~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_genotype_matrix_no_missing_13670.csv", row.names = TRUE)

# Or save it as an RDS file to preserve data structure
#saveRDS(complete_markers, "~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_genotype_matrix_no_missing_13670.rds")

############################ 
#full try of per snps and individual marker effect size calculation - with no missingness
############################ 
library(vegan)
library(adegenet)
library(psych)
library(rrBLUP)
library(dplyr)
library(tidyr)

# Load genotype matrix with complete markers (no missing values)
genotype_matrix <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_genotype_matrix_no_missing_13670.csv', row.names = 1)

# Load phenotype data
#envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8X_core_n50.csv', header = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8x_prec_zero_replaced.csv', header = TRUE)

# Define the traits of interest
#trait_january <- "wc2.1_30s_prec_01"
trait_june <- "wc2.1_30s_prec_06"
trait_august <- "wc2.1_30s_prec_08"

# Reload the genotype_matrix with the first column as row names
genotype_matrix <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_genotype_matrix_no_missing_13670.csv', row.names = 1)
rownames(genotype_matrix) <- genotype_matrix$X  # Set the X column as row names
genotype_matrix <- genotype_matrix[ , -1]       # Remove the X column

# Define the core samples (n=50)
core_samples <- c("SRR14708197", "SRR14708199", "SRR14708202", "SRR14708203", "SRR14708204", 
                  "SRR14708205", "SRR14708206", "SRR14708207", "SRR14708208", "SRR14708209", 
                  "SRR14708210", "SRR14708212", "SRR14708213", "SRR14708229", "SRR14708230", 
                  "SRR14708231", "SRR14708234", "SRR14708235", "SRR14708236", "SRR14708243", 
                  "SRR14708246", "SRR14708248", "SRR14708251", "SRR14708275", "SRR14708276", 
                  "SRR9089225.1", "SRR9089228.1", "SRR9089230.1", "SRR9089236.1", "SRR9089253.1", 
                  "SRR9089257.1", "SRR9089259.1", "SRR9089261.1", "SRR9089262.1", "SRR9089263.1", 
                  "SRR9089266.1", "SRR9089270.1", "SRR9089275.1", "SRR9089276.1", "SRR9089282.1", 
                  "SRR9089285.1", "SRR9089295.1", "SRR9089297.1", "SRR9089299.1", "SRR9089300.1", 
                  "SRR9089301.1", "SRR9089303.1", "SRR9089304.1", "SRR9089308.1", "SRR9089310.1")

# Ensure there are no leading/trailing spaces in row names
rownames(genotype_matrix) <- trimws(rownames(genotype_matrix))

# Filter genotype_matrix to include only the core samples
genotype_core <- genotype_matrix[rownames(genotype_matrix) %in% core_samples, , drop = FALSE]

# Check dimensions to confirm filtering
cat("Dimensions of genotype_core after filtering:", dim(genotype_core), "\n")

envdat_core <- envdat[envdat$Sample_ID %in% core_samples, ]


# Reorder envdat_core to match the rownames of genotype_core
envdat_core <- envdat_core[match(rownames(genotype_core), envdat_core$Sample_ID), ]

# Check alignment again
if (!all(rownames(genotype_core) == envdat_core$Sample_ID)) {
  stop("Sample order mismatch still exists. Please check data.")
} else {
  cat("Sample order in envdat_core now matches genotype_core.\n")
}

genotype_core <- as.matrix(genotype_core)  # Ensure it's a matrix

# Run RR-BLUP on core individuals to get SNP effects for JUNE and August
solve_june <- mixed.solve(y = envdat_core[[trait_june]], Z = genotype_core)
u.hat_june <- solve_june$u  # SNP effects for JUNE

solve_aug <- mixed.solve(y = envdat_core[[trait_august]], Z = genotype_core)
u.hat_aug <- solve_aug$u  # SNP effects for August

# Apply SNP effects to all individuals in complete_markers
complete_markers <- as.matrix(genotype_matrix)
complete_markers <- as.matrix(complete_markers)  # Ensure matrix format

# Calculate SNP-level contributions for each individual for January and August
individual_marker_effects_june <- complete_markers * matrix(u.hat_june, nrow = nrow(complete_markers), ncol = ncol(complete_markers), byrow = TRUE)
individual_marker_effects_aug <- complete_markers * matrix(u.hat_aug, nrow = nrow(complete_markers), ncol = ncol(complete_markers), byrow = TRUE)

# Convert to data frames for easier interpretation
individual_marker_effects_june_df <- as.data.frame(individual_marker_effects_june)
colnames(individual_marker_effects_june_df) <- colnames(complete_markers)
individual_marker_effects_june_df$Individual <- rownames(complete_markers)
individual_marker_effects_june_df$Month <- "June"

individual_marker_effects_aug_df <- as.data.frame(individual_marker_effects_aug)
colnames(individual_marker_effects_aug_df) <- colnames(complete_markers)
individual_marker_effects_aug_df$Individual <- rownames(complete_markers)
individual_marker_effects_aug_df$Month <- "August"

# Combine January and August data into one data frame
individual_marker_effects_combined <- rbind(individual_marker_effects_june_df, individual_marker_effects_aug_df)

# Reshape to long format for per-SNP effect per individual
individual_marker_effects_long <- individual_marker_effects_combined %>%
  pivot_longer(
    cols = -c(Individual, Month),
    names_to = "SNP",
    values_to = "Effect"
  )

#per snp and individual marker effects
write.csv(individual_marker_effects_long, '~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_prec_per_SNP_june_aug.csv', row.names = FALSE) 

# Print a summary to check if both months are included
table(individual_marker_effects_combined$Month)


#################
#Figure 4B
#################
#signifiant marker effect sizes for dt-feral in cwa  - keeps only significants and plots the mean effect size over all the individuals in the plot (PURPLE)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the data
individual_marker_effects_long <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_prec_per_SNP_june_aug.csv')

# Define selected individuals and chromosome order
selected_individuals <- c("SRR14708198", "SRR14708199", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258")
chromosome_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", "NC_044373.1", "NC_044374.1", 
                      "NC_044377.1", "NC_044378.1", "NC_044379.1", "NC_044376.1", "NC_044370.1")

# Update chromosome labels to include both chromosome number and NC identifier
chromosome_mapping <- data.frame(
  original = c("NC_044371.1", "NC_044375.1", "NC_044372.1", "NC_044373.1", "NC_044374.1", 
               "NC_044377.1", "NC_044378.1", "NC_044379.1", "NC_044376.1", "NC_044370.1"),
  label = c("Chromosome 1 (NC_044371.1)", 
            "Chromosome 2 (NC_044375.1)", 
            "Chromosome 3 (NC_044372.1)", 
            "Chromosome 4 (NC_044373.1)", 
            "Chromosome 5 (NC_044374.1)", 
            "Chromosome 6 (NC_044377.1)", 
            "Chromosome 7 (NC_044378.1)", 
            "Chromosome 8 (NC_044379.1)", 
            "Chromosome 9 (NC_044376.1)", 
            "Chromosome 10 (X) (NC_044370.1)")
)

# Filter and process data
data_long <- individual_marker_effects_long %>%
  filter(Individual %in% selected_individuals) %>%
  mutate(
    Chromosome = factor(sub("_[0-9]+$", "", SNP), levels = chromosome_order),
    Position = as.numeric(sub(".*_", "", SNP))
  ) %>%
  pivot_wider(names_from = Month, values_from = Effect) %>%
  drop_na(June, August)  # Remove rows with missing values

# Filter for significant SNPs based on adjusted p-value
test_results <- data_long %>%
  group_by(SNP, Chromosome, Position) %>%
  summarize(
    t_stat = if (length(June) > 1 && length(August) > 1 && sd(June) != 0 && sd(August) != 0) {
      t.test(June, August, paired = TRUE)$statistic
    } else {
      NA
    },
    p_value = if (length(June) > 1 && length(August) > 1 && sd(June) != 0 && sd(August) != 0) {
      t.test(June, August, paired = TRUE)$p.value
    } else {
      NA
    },
    .groups = 'drop'
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "BH"),
    significant = adjusted_p_value < 0.01
  )

# Keep only significant SNPs
significant_snps_only <- data_long %>%
  filter(SNP %in% test_results$SNP[test_results$significant == TRUE]) %>%
  group_by(SNP, Chromosome, Position) %>%
  summarize(
    Mean_August = mean(August, na.rm = TRUE)
  ) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chromosome_mapping$original, 
                        labels = chromosome_mapping$label),
    Position = Position / 1e6  # Convert to Mb
  )

# Generate the plot
plot_significant_snps <- ggplot(significant_snps_only, aes(x = Position, y = Mean_August)) +
  geom_point(color = "purple3", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 2) +
  labs(
    title = "",
    x = "Genomic Position (Mb)",
    y = "Mean Effect Size"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

# Display the plot
print(plot_significant_snps)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure4B.pdf",
  plot = plot_significant_snps,
  device = "pdf",
  width = 12,
  height = 8,
  units = "in"
)

#################
#Figure 4C
#################
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the data
individual_marker_effects_long <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_prec_per_SNP_june_aug.csv')
selected_individuals <- c("SRR14708198", "SRR14708199", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258")                # drug-type feral in Cwa

# Define chromosome order for consistent plotting
chromosome_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", "NC_044373.1", "NC_044374.1", 
                      "NC_044377.1", "NC_044378.1", "NC_044379.1", "NC_044376.1", "NC_044370.1")

# Step 1: Filter data for the selected individuals
data_long <- individual_marker_effects_long %>%
  filter(Individual %in% selected_individuals) %>%
  mutate(
    Chromosome = factor(sub("_[0-9]+$", "", SNP), levels = chromosome_order),
    Position = as.numeric(sub(".*_", "", SNP))
  )

# Step 2: Reshape the data from long to wide format with separate columns for May and August
data_wide <- data_long %>%
  pivot_wider(names_from = Month, values_from = Effect) %>%
  drop_na(June, August)  # Remove rows with NA in either June or August


test_results <- data_wide %>%
  group_by(SNP, Chromosome, Position) %>%
  summarize(
    t_stat = if (all(!is.na(June), !is.na(August)) && length(June) > 1 && length(August) > 1 &&
                 sd(June) != 0 && sd(August) != 0) {
      t.test(June, August, paired = TRUE)$statistic
    } else {
      NA
    },
    p_value = if (all(!is.na(June), !is.na(August)) && length(June) > 1 && length(August) > 1 &&
                  sd(June) != 0 && sd(August) != 0) {
      t.test(June, August, paired = TRUE)$p.value
    } else {
      NA
    },
    .groups = 'drop'
  ) %>%
  filter(!is.na(p_value)) %>%  # Remove rows with NA p-values
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "BH"),  # Adjust for multiple testing
    significant = adjusted_p_value < 0.01  # Define significance level
  )


# Step 1: Filter only significant SNPs (788 SNPs)
significant_snps_only <- data_wide %>%
  filter(SNP %in% test_results$SNP[test_results$significant == TRUE]) %>%
  group_by(SNP, Chromosome, Position) %>%
  summarize(
    Mean_August = mean(August, na.rm = TRUE),
    Mean_June = mean(June, na.rm = TRUE),
    Mean_Difference = Mean_August - Mean_June
  )

##################################################################################################################
# ID SNPs in annotation 
library(GenomicRanges)
library(rtracklayer)

all_significant_snps <- significant_snps_only %>%
  mutate(
    Chromosome = sub("(_[0-9]+)$", "", SNP),  # Remove only the position part after the underscore
    Position = as.numeric(sub(".+_", "", SNP))  # Extract the position as numeric
  )

# Step 2: Convert all significant SNPs to GRanges
snp_gr <- GRanges(
  seqnames = all_significant_snps$Chromosome,  # Use only the chromosome identifier here
  ranges = IRanges(start = all_significant_snps$Position, end = all_significant_snps$Position),
  mcols = all_significant_snps
)

# Step 4: Import the GFF file as GRanges
annotation_gr <- import("~/R/cannabis_GEAV/Ren_Soorni/GEA/genomic_main_chromosomes.gff")

# Check if seqlevels match
cat("Chromosomes in SNP data (snp_gr):\n")
print(unique(seqlevels(snp_gr)))
cat("Chromosomes in annotation data (annotation_gr):\n")
print(unique(seqlevels(annotation_gr)))

# If they match, proceed with finding overlaps
if (all(unique(seqlevels(snp_gr)) %in% unique(seqlevels(annotation_gr)))) {
  # Step 5: Find overlaps between SNPs and annotation
  overlaps <- findOverlaps(snp_gr, annotation_gr)
  
  # Step 6: Extract overlapping SNPs and annotation features
  overlapping_snps <- snp_gr[queryHits(overlaps)]
  overlapping_features <- annotation_gr[subjectHits(overlaps)]
  
  # Define column names based on the available columns in the annotation
  snp_col_name <- if ("SNP" %in% names(mcols(overlapping_snps))) "SNP" else names(mcols(overlapping_snps))[1]
  feature_id_col <- if ("ID" %in% names(mcols(overlapping_features))) "ID" else names(mcols(overlapping_features))[1]
  feature_type_col <- if ("type" %in% names(mcols(overlapping_features))) "type" else names(mcols(overlapping_features))[2]
  
  # Combine the SNP and annotation information into a data frame
  snp_annotation_df <- data.frame(
    SNP_ID = mcols(overlapping_snps)[[snp_col_name]],  # Extract SNP ID
    SNP_Chromosome = as.character(seqnames(overlapping_snps)),
    SNP_Position = start(overlapping_snps),
    Feature_ID = mcols(overlapping_features)[[feature_id_col]],  # Extract feature ID
    Feature_Type = mcols(overlapping_features)[[feature_type_col]]  # Extract feature type
  )
  
  # View and save the annotated SNPs
  print(head(snp_annotation_df))
  write.csv(snp_annotation_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_High_Effect_june_aug.csv", row.names = FALSE)
} else {
  cat("Chromosome names still do not match. Please inspect manually.\n")
}

##################################################################################################################
# Filter for rows where Feature_Type is "gene"
gene_annotations <- snp_annotation_df %>% filter(Feature_Type == "gene")

# Count the unique genes
num_unique_genes <- gene_annotations %>% pull(Feature_ID) %>% unique() %>% length()

# Display the result
cat("Number of unique genes with overlapping high-effect SNPs:", num_unique_genes, "\n")
#Number of unique genes with overlapping high-effect SNPs: 566 

##################################################################################################################
#Pull LOCs in genes
# Load necessary libraries
library(dplyr)
library(stringr)

# Read in your CSV file
snp_annotation_df <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_High_Effect_june_aug.csv")

# Filter rows where Feature_Type is "gene" to focus only on genes
gene_rows <- snp_annotation_df %>% filter(Feature_Type == "gene")

# Extract and clean LOC IDs
gene_rows <- gene_rows %>%
  mutate(LOC_ID = str_remove(Feature_ID, "gene-")) %>%  # Remove "gene-" prefix
  distinct(LOC_ID)  # Keep only unique LOC IDs (xxxx)

# View the unique LOC IDs
print(gene_rows$LOC_ID)

# Optionally, save the unique LOC IDs to a new CSV file
write.csv(gene_rows$LOC_ID, "~/R/cannabis_GEAV/Ren_Soorni/GEA/Unique_LOC_IDs_june_august.csv", row.names = FALSE)

# Get unique gene IDs from the filtered gene annotations
#unique_gene_ids <- gene_annotations %>% pull(Feature_ID) %>% unique()

###################
#Get what LOCs associated with in Genes with description using rentrez
library(dplyr)
library(stringr)
library(rentrez)

# Read in your CSV file
snp_annotation_df <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_High_Effect_june_aug.csv")

# Extract unique gene IDs (LOC IDs) from the 'Feature_ID' column where 'Feature_Type' is "gene"
gene_ids <- snp_annotation_df %>%
  filter(Feature_Type == "gene") %>%
  mutate(LOC_ID = str_remove(Feature_ID, "gene-")) %>%
  pull(LOC_ID) %>%
  unique()

print(paste("Number of unique gene IDs:", length(gene_ids)))
if (length(gene_ids) == 0) {
  stop("No gene IDs found. Check the CSV file and the Feature_Type filtering.")
}

# Set email for Entrez
Entrez.email <- "mccormick.anna@hotmail.com"

# Define batch size for Entrez queries
batch_size <- 100

# Initialize a data frame to store results
results <- data.frame(GeneID=character(), Name=character(), Description=character(), Chromosome=character(),
                      OtherAliases=character(), GenomicInfo=character(), stringsAsFactors=FALSE)

# Loop through gene IDs in batches to query Entrez
for (i in seq(1, length(gene_ids), by=batch_size)) {
  # Select the batch of gene IDs
  batch <- gene_ids[i:min(i + batch_size - 1, length(gene_ids))]
  query_string <- paste(batch, collapse=" OR ")
  
  # Search for gene summaries
  search_results <- entrez_search(db="gene", term=query_string, retmax=100)
  
  print(paste("Batch", i, ":", length(search_results$ids), "results found"))
  
  if (length(search_results$ids) > 0) {
    # Fetch the summaries
    summaries <- entrez_summary(db="gene", id=search_results$ids)
    
    # Ensure summaries is a list of lists, even for a single summary
    if (!is.list(summaries)) {
      summaries <- list(summaries)
    }
    
    # Parse and store the results
    for (summary in summaries) {
      # Extract genomic info with additional checks
      if (is.list(summary$genomicinfo) && length(summary$genomicinfo) > 0 && is.data.frame(summary$genomicinfo[[1]])) {
        genomic_info <- paste(
          summary$genomicinfo[[1]]$chrloc %||% NA,
          summary$genomicinfo[[1]]$chraccver %||% NA,
          summary$genomicinfo[[1]]$chrstart %||% NA,
          summary$genomicinfo[[1]]$chrstop %||% NA,
          summary$genomicinfo[[1]]$exoncount %||% NA,
          sep="; "
        )
      } else {
        genomic_info <- NA
      }
      
      # Append data to results
      results <- rbind(results, data.frame(
        GeneID=summary$uid,
        Name=summary$name,
        Description=summary$description,
        Chromosome=summary$chromosome,
        OtherAliases=summary$otheraliases,
        GenomicInfo=genomic_info,
        stringsAsFactors=FALSE
      ))
    }
  }
  Sys.sleep(0.5)  # Rate limit
}

# Write the results to CSV
write.csv(results, "~/R/cannabis_GEAV/Ren_Soorni/GEA/gene_descriptions_dtferal_june_aug.csv", row.names=FALSE) #n=566 (pulls descriptions of cannabis LOC to gene - using rEntrez) Input for OrthoDB

#######################
#OrthoDB
#######################
#devtools::install_git("https://gitlab.com/ezlab/orthodb_r.git")       #need R 4.4.0 as a minimum to run locally
# on Koa HPC
#srun --mem=20G -t 0-04:00:00 --pty bash   
#module load lang/Anaconda3
#source activate new_R_env_4_4_0
#R
library(OrthoDB)
# Initialize API connection
api <- OrthoDB::OdbAPI$new("v12")

# Load your input CSV with all query GeneIDs
#gene_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/OrthoDB/gene_descriptions_dtferal.csv")
gene_data <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/OrthoDB/OrthoDB_may_aug/gene_descriptions_dtferal_may_aug.csv") #exact same 566 as above
gene_ids <- gene_data$GeneID  # Adjust if the column name is different

# Initialize a list to store results for each query
all_results <- list()

# Loop through each GeneID query
for (query_id in gene_ids) {
  cat("Processing GeneID:", query_id, "\n")
  
  # Search for orthologs for the current GeneID
  gs <- tryCatch({
    api$gene_search(as.character(query_id))
  }, error = function(e) {
    cat("Error occurred during gene search for query_id:", query_id, "\n")
    return(NULL)
  })
  
  # Validate response
  if (is.null(gs)) {
    cat("No valid response for query_id:", query_id, "\n")
    all_results[[as.character(query_id)]] <- NA
    next
  }
  
  # Get orthologs across model organisms
  orthologs <- gs$orthologs_in_model_organisms
  
  # Filter only for Nicotiana tabacum results
  tobacco_orthologs <- lapply(orthologs, function(ortholog_set) {
    if (!is.null(ortholog_set$organism$name) && ortholog_set$organism$name == "Nicotiana tabacum") {
      return(ortholog_set$df)  # Extract the data frame for tobacco orthologs
    } else {
      return(NULL)
    }
  })
  
  # Remove any NULLs and combine the tobacco ortholog results
  tobacco_orthologs <- Filter(Negate(is.null), tobacco_orthologs)
  tobacco_orthologs_df <- do.call(rbind, tobacco_orthologs)
  
  # Check if there are valid tobacco orthologs
  if (!is.null(tobacco_orthologs_df) && "genes" %in% names(tobacco_orthologs_df)) {
    # Extract gene IDs within the genes column
    extracted_gene_ids <- lapply(tobacco_orthologs_df$genes, function(gene_data) {
      if (!is.null(gene_data$gene_id$id)) {
        return(gene_data$gene_id$id)
      } else {
        return(NA)
      }
    })
    
    # Flatten and clean up the list
    all_results[[as.character(query_id)]] <- unlist(extracted_gene_ids)
  } else {
    cat("No valid Nicotiana tabacum orthologs found for query_id:", query_id, "\n")
    all_results[[as.character(query_id)]] <- NA
  }
  
  # Delay to avoid rate limiting
  Sys.sleep(1)
}


# Convert the list of results to a data frame with each GeneID as a separate column
max_length <- max(sapply(all_results, length))  # Find the maximum number of matches for alignment
results_df <- as.data.frame(lapply(all_results, function(x) {
  # Extend each list to max_length with NA to make them align as columns
  c(x, rep(NA, max_length - length(x)))
}))

# Save the results to a CSV file
#write.csv(results_df, "/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/OrthoDB/tobacco_gene_matches.csv", row.names = FALSE)
#write.csv(results_df, "/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/OrthoDB/OrthoDB_1814/tobacco_gene_matches_for1814_1736.csv", row.names = FALSE)
write.csv(results_df, "/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/OrthoDB/OrthoDB_may_aug/tobacco_gene_matches_may_aug.csv", row.names = FALSE)
#512 matches found in tobacco 

#############
#Take orthologs first hits and paste into DAVID
#Select the BP,CC,MF and KEGG boxes for output from DAVID 
#Manually copy and paste the functional annotation data with the p values etc.
#############
library(ggplot2)

data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/functional_annotation_david_may_aug.csv")

# Calculate -log10(P-Value) and add it as a new column
data$LogPValue <- -log10(data$P.Value)

# Reorder data by Category and -log10(P-Value)
data$Category <- factor(data$Category, levels = c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
data <- data[order(data$Category, -data$LogPValue), ]

# Convert Term to a factor to maintain the order in the plot
data$Term <- factor(data$Term, levels = unique(data$Term))

# Plot the reordered data with larger circles and updated legend size
functional_plot <- ggplot(data, aes(x = LogPValue, y = Term, color = Category)) +
  geom_point(aes(size = Count)) +
  scale_color_manual(
    values = c("GOTERM_BP_DIRECT" = "blue", 
               "GOTERM_CC_DIRECT" = "red", 
               "GOTERM_MF_DIRECT" = "orange",
               "KEGG_PATHWAY" = "green"),  # New color for KEGG_PATHWAY
    labels = c("GOTERM_BP_DIRECT" = "Biological Process", 
               "GOTERM_CC_DIRECT" = "Cellular Component", 
               "GOTERM_MF_DIRECT" = "Molecular Function",
               "KEGG_PATHWAY" = "KEGG Pathway")  # Label for KEGG_PATHWAY
  ) +
  scale_size_continuous(range = c(5, 15)) +  # Adjust range to make circles larger
  labs(
    title = "", 
    x = "-log10(P-Value)", 
    y = "Enrichment Term",
    color = "Category"
  ) +  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 14),   # Increase legend text size
    axis.text = element_text(size = 12),     # Adjust axis text size
    axis.title = element_text(size = 14),    # Adjust axis title size
    plot.title = element_text(size = 16, face = "bold")  # Adjust title size
  )

# Display the plot
print(functional_plot)

# Save the plot as a PDF
ggsave(
  filename = "~/R/cannabis_GEAV/Ren_Soorni/FIGS/Figure4C.pdf",  # Specify the file path and name
  plot = functional_plot,                                                      # 
  device = "pdf",                                                              # 
  width = 12,                                                                  # 
  height = 12,                                                                  # 
  units = "in"                                                                 # 
)


###################################################################################################################################################
#SUPPLEMENTAL FIGURES
###################################################################################################################################################

#########################
#Supplemental Figure 1
##########################
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


#########################
#Supplemental Figure 2
##########################
#Soorni alone
library(SNPRelate)
vcf.fn <- "~/R/cannabis_GEAV/MERGE2/Cannabis_sativa_PRJNA419020_Soorni_filtered.vcf.gz"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "soorni.gds", method="copy.num.of.ref")
snpgdsSummary("soorni.gds")
genofile <- snpgdsOpen("soorni.gds")
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

##########
#Figure S2A
##########
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
write.csv(tab, "~/R/cannabis_GEAV/MERGE2/soorni_alone_PC.csv")



pca <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_alone_PC.csv")
p <- ggplot(pca, aes(x = EV1, y = EV2, color = Population)) +
  geom_point(size = 3) +  # Create scatter plot points
  labs(title = "", 
       x = paste0("PC1 (3.70%)"), 
       y = paste0("PC2 (2.01%)")) +
  theme_minimal() +
  theme(legend.position = "right") +  # Add legend on the right
  scale_color_viridis_d(option = "viridis", name = "Climate Class")  # Use viridis color scale for better differentiation

p <- ggplot(pca, aes(x = EV1, y = EV2, color = Population)) +
  geom_point(size = 3) +  # Create scatter plot points
  labs(title = "", 
       x = paste0("PC1 (3.70%)"), 
       y = paste0("PC2 (2.01%)")) +
  theme_minimal() +
  theme(legend.position = "right") +  # Add legend on the right
  scale_color_viridis_d(option = "viridis", name = "Climate Class",
                        labels = gsub("_", " ", unique(pca$Population)))  # Remove underscores
p


# Save the plot as a PDF
ggsave(filename = "~/R/cannabis_GEAV/FIGS/Figure_S2A.pdf", plot = p, device = "pdf", width = 8, height = 6)


##########
#Figure S2B
##########
#MAP
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

data <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_alone_PC.csv")
iran_data <- data[data$Latitude >= 25 & data$Latitude <= 40 & data$Longitude >= 44 & data$Longitude <= 64, ]

world <- ne_countries(scale = "medium", returnclass = "sf")
iran <- world[world$name == "Iran", ]

p <- ggplot() +
  geom_sf(data = iran, fill = "lightgray", color = "black") +  # Base map of Iran
  geom_point(data = iran_data, aes(x = Longitude, y = Latitude, color = Population), size = 3, alpha = 0.8) +  # Plot points, colored by Population
  labs(title = "", x = "Longitude", y = "Latitude", color = "Population") +
  theme_minimal() +
  scale_color_manual(
    values = c("Population_1" = "purple", "Population_2" = "green", "Population_3" = "yellow"),
    labels = c("Population_1" = "Population 1", "Population_2" = "Population 2", "Population_3" = "Population 3")
  ) +  # Custom colors with labels
  coord_sf(xlim = c(44, 64), ylim = c(25, 40), expand = FALSE)  # Set coordinate limits to Iran

p
ggsave(filename = "~/R/cannabis_GEAV/FIGS/Figure_S2B.pdf", plot = p, device = "pdf", width = 8, height = 6)


##########
#Figure S2C
##########
#Type
p<- ggplot() +
  geom_sf(data = iran, fill = "lightgray", color = "black") +  # Base map of Iran
  geom_point(data = iran_data, aes(x = Longitude, y = Latitude, color = Type), size = 3, alpha = 0.8) +  # Plot points, colored by Population
  labs(title = "", x = "Longitude", y = "Latitude", color = "Type") +
  theme_minimal() +
  scale_color_manual(values = c("Type I" = "red", "Type II" = "orange", "Type III" = "green")) +  # Custom colors
  coord_sf(xlim = c(44, 64), ylim = c(25, 40), expand = FALSE)  # Set coordinate limits to Iran
p
ggsave(filename = "~/R/cannabis_GEAV/FIGs/Figure_S2C.pdf", plot = p, device = "pdf", width = 8, height = 6)


##########
#Figure S2D
##########
library(ggplot2)
#Climate niches in this dataset
data <- read.csv("~/R/cannabis_GEAV/MERGE2/soorni_alone_PC.csv")

data_summary <- data %>%
  group_by(Population, Climate_Class_Description) %>%
  summarise(count = n())

# Create the plot
ggplot(data_summary, aes(x = Population, y = count, fill = Climate_Class_Description)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot
  geom_text(aes(label = Climate_Class_Description), 
            position = position_fill(vjust = 0.5), 
            color = "white", size = 3) +  # Add climate class labels in white
  labs(title = "Distribution of Climate Classes by Population", 
       x = "Population", 
       y = "Proportion", 
       fill = "Climate Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_viridis_d(option = "viridis")  # Use a color palette





##########################
#Supplemental Figure 3
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

#Figure S3B
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

#Figure S3C
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
p<- ggplot(pca, aes(x = EV1, y = EV2, color = Type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, size = 1) +  # Add ellipses
  scale_color_manual(values = c("Type I" = "red",
                                "Type II" = "orange",
                                "Type III" = "green", 
                                "unknown" = "grey")) +  # Assign colors
  labs(title = "Secondary Chemistry", 
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
ggsave("~/R/cannabis_GEAV/FIGS/PCA_1D.pdf", plot = p, device = "pdf", width = 8, height = 6)


##########################
#Supplemental Figure 4
##########################
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


##########################
#Supplemental Figure 5
##########################
#56k tree
#Have a supplemental walk through document called "Hops_Cannabis_Tree_Steps" that has more exact details 
#convert vcf to phylip format and remove invariant sites (we had 18)
#match cannabis snp positions in hops for outgroup
#iqtree for model assessment for bes substitution model and included acertainment bias as this is a SNP set randomly spread across the genome
#after iqtree then add in names below with full metadata name and then import into FigTree
###############################################################################################
#TREE 12
###############################################################################################
library(ape)      # For tree manipulation
library(ggtree)   # For visualization
library(dplyr)    # For data manipulation
library(ggplot2)  # For additional visualization

# Load the tree file
tree <- read.tree("~/R/cannabis_GEAV/Ren_Soorni/Tree/Tree_12/filtered_variants.min4.phy.varsites.phy.contree")
class(tree)

# Load metadata
metadata <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/Tree/Tree_10/Cannabis_EGS_149_names.csv")
colnames(metadata) <- c("Current", "Full")

# Replace SRR IDs with full names
tree$tip.label <- metadata$Full[match(tree$tip.label, metadata$Current)]

# Assign groups, ensuring drug_type and hemp_type are separate
metadata <- metadata %>%
  mutate(
    Group = case_when(
      grepl("drug", Full, ignore.case = TRUE) ~ "drug_type",
      grepl("hemp", Full, ignore.case = TRUE) ~ "hemp_type",
      TRUE ~ sub(".*_", "", Full)  # Keep existing classifications
    )
  )

# Merge tree labels with metadata
tree_data <- data.frame(label = tree$tip.label) %>%
  left_join(metadata, by = c("label" = "Full"))

str(tree)

# Order tree from high to low (ladderize)
tree <- ladderize(tree, right = FALSE)

# Find the outgroup and root the tree
root_sample <- grep("Humulus_lupulus", tree$tip.label, value = TRUE)
tree <- root(tree, outgroup = root_sample, resolve.root = TRUE)

write.tree(tree, file = "~/R/cannabis_GEAV/Ren_Soorni/Tree/Tree_12/tree_12_processed.contree")


ggtree(tree) %<+% tree_data +
  geom_tiplab(size = 2) +  
  geom_tippoint(aes(color = Group), size = 2) +  
  scale_color_manual(values = c(
    "1" = "#9370DB",
    "2" = "#008080",
    "basal" = "orange",
    "feral" = "blue",
    "lupulus" = "black",
    "drug_type" = "red",
    "hemp_type" = "green"
  )) +
  theme_tree2() +  
  scale_x_continuous(
    trans = "log1p",  # Log transformation to expand small values
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 10, 15),  
    labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "1", "10", "15"),
    expand = c(0.2, 0.5)  # Adjust spacing for better readability
  ) +
  theme(legend.position = "right")



##########################
#Supplemental Figure 6
##########################
########
#Figure 6A - MAP
########
# Load necessary libraries
library(ggplot2)
library(maps)

# Load the world map data
world_map <- map_data("world")

# Load your data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_Lat_Long_n111_climate_class.csv")

# Create the plot with the world map, limit longitude and latitude, and color by Phylogenetic_tree
p <- ggplot() +
  # Plot the world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  # Add points for your latitude and longitude data, colored by Phylogenetic_tree
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Phylogenetic_tree), 
             size = 3, alpha = 0.7) +
  labs(title = "Geographical Distribution of Phylogenetic Categories", 
       x = "Longitude", 
       y = "Latitude") + 
  scale_color_manual(
    values = c(
      "basal" = "orange", 
      "drug-type feral" = "blue", 
      "hemp-type" = "green", 
      "Population_1" = "#9370DB",
      "Population_2" = "#008080"
    ),
    labels = c(
      "basal" = "Basal", 
      "drug-type feral" = "Drug-type feral", 
      "hemp-type" = "Hemp-type", 
      "Population_1" = "Population-1",
      "Population_2" = "Population-2"
    )
  ) +  # Customize colors and legend labels
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  # Limit the longitude and latitude range to focus on the region you want
  coord_cartesian(xlim = c(-10, 150), ylim = c(10, 80))
p

ggsave("~/R/cannabis_GEAV/FIGS/tree_ranges.pdf", plot = p, width = 12, height = 5)

#########################
#Supplemental Figure 6B
##########################
#Boxplots of raw data ranges for the groups

#Extract relevant environmental data from WorldClim (Bioclimatic 19X variables) with latitude and longitudes 
# For n50 for core greedy
library(raster)
worldclim_path <- "~/R/Data/climate_data/wc2.0_30s_bio/"
bio_files <- paste0(worldclim_path, "bio_", sprintf("%02d", 1:19), ".tif")
climate_layers <- stack(bio_files)

# Load your latitude and longitude dataset
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy.csv")
head(data)

# Prepare the coordinates (assuming Latitude and Longitude columns are present)
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "~/R/cannabis_GEAV/Ren_Soorni/outputs/GEBVs_for_n111_core_greedy_climate_data_raw.csv", row.names = FALSE)

########
#Boxplots 
########
library(dplyr)
library(tidyr)  # For pivot_longer
library(ggplot2)

# Read data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/GEBVs_for_n111_core_greedy_climate_data_raw.csv")

# List of WorldClim variables
worldclim_vars <- c("bio_01", "bio_02", "bio_03", "bio_04", "bio_05", "bio_06", "bio_07", 
                    "bio_08", "bio_09", "bio_10", "bio_11", "bio_12", "bio_13", 
                    "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")

# Reshape the data into long format
data_long <- data %>%
  pivot_longer(cols = all_of(worldclim_vars), 
               names_to = "WorldClim_variable", 
               values_to = "Value")

# Set the desired order for Phylogenetic_tree levels
data_long$Phylogenetic_tree <- factor(data_long$Phylogenetic_tree, 
                                      levels = c("basal", "hemp-type", "drug-type feral", "Population_1", "Population_2"))

# Updated facet labels with "BIO" and descriptive breakdown
worldclim_labels <- c(
  bio_01 = "BIO 01: Annual Mean Temperature (°C)",
  bio_02 = "BIO 02: Mean Diurnal Range",
  bio_03 = "BIO 03: Isothermality",
  bio_04 = "BIO 04: Temperature Seasonality (SD)",
  bio_05 = "BIO 05: Max Temp of Warmest Month (°C)",
  bio_06 = "BIO 06: Min Temp of Coldest Month (°C)",
  bio_07 = "BIO 07: Temperature Annual Range",
  bio_08 = "BIO 08: Mean Temp of Wettest Quarter (°C)",
  bio_09 = "BIO 09: Mean Temp of Driest Quarter (°C)",
  bio_10 = "BIO 10: Mean Temp of Warmest Quarter (°C)",
  bio_11 = "BIO 11: Mean Temp of Coldest Quarter (°C)",
  bio_12 = "BIO 12: Annual Precipitation (mm)",
  bio_13 = "BIO 13: Precipitation of Wettest Month (mm)",
  bio_14 = "BIO 14: Precipitation of Driest Month (mm)",
  bio_15 = "BIO 15: Precipitation Seasonality (CV)",
  bio_16 = "BIO 16: Precipitation of Wettest Quarter (mm)",
  bio_17 = "BIO 17: Precipitation of Driest Quarter (mm)",
  bio_18 = "BIO 18: Precipitation of Warmest Quarter (mm)",
  bio_19 = "BIO 19: Precipitation of Coldest Quarter (mm)"
)

# Custom color mapping for phylogenetic trees
color_mapping <- c("basal" = "orange", 
                   "drug-type feral" = "blue", 
                   "hemp-type" = "green", 
                   "Population_1" = "#9370DB", 
                   "Population_2" = "#008080")

# Updated labels for legend
legend_labels <- c(
  "basal" = "Basal", 
  "drug-type feral" = "Drug-type feral", 
  "hemp-type" = "Hemp-type", 
  "Population_1" = "Population-1",
  "Population_2" = "Population-2"
)

# Faceted boxplot with custom labels, colors, and updated legend
p <- ggplot(data_long, aes(x = Phylogenetic_tree, y = Value, fill = Phylogenetic_tree)) +
  geom_boxplot() +
  facet_wrap(~ WorldClim_variable, scales = "free_y", 
             labeller = as_labeller(worldclim_labels)) +
  scale_fill_manual(values = color_mapping, labels = legend_labels) +  # Apply the custom color mapping and legend labels
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

# Display the plot
p

ggsave("~/R/cannabis_GEAV/FIGS/worldclim_boxplots.pdf", plot = p, width = 20, height = 10)

#########################
#Supplemental Figure 7
##########################
########################
#Figure S7A
######################
#overall climate classes captured by datasets
#####################
library(ggplot2)
library(dplyr)

# Load the data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy_trimmed_for_climate.csv", header = TRUE)

# Filter for datasets
filtered_data <- data %>% 
  filter(Dataset %in% c("Ren et al 2021", "Soorni et al 2017"))

# Create the plot
p<- ggplot(filtered_data, aes(x = Climate_Class_Description)) +
  geom_bar(position = "dodge", fill = "black") +  # Use black fill for bars
  geom_hline(yintercept = 5, color = "red", linetype = "dashed", size = 1) +  # Add a red dashed line at y = 5
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_minimal() +
  labs(title = "Climate Class Distribution for Ren et al 2021 and Soorni",
       x = "Climate Class",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # Adjust text angle and size for readability
        strip.text = element_text(size = 12))  # Adjust facet title size

ggsave("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Climate_niche_covered_by_datasets.pdf", plot = p, width = 10, height = 8)

########################
#Figure S7B
########################
#niches captured by core
########################
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy_trimmed_for_climate.csv", header = TRUE)
core_data <- data %>% filter(Core_greedy == TRUE)

# Summarize the counts for each unique climate class by counting rows and assign a numeric label
climate_summary <- core_data %>%
  group_by(Climate_Class_Description) %>%
  summarise(Total_Count = n()) %>%
  mutate(Label = row_number())  # Add a numeric label

# Plot the climate class counts with numeric labels
p<- ggplot(climate_summary, aes(x = reorder(as.factor(Label), -Total_Count), y = Total_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  scale_x_discrete(labels = paste(climate_summary$Label, climate_summary$Climate_Class_Description, sep = " - ")) +
  scale_y_continuous(breaks = seq(0, max(climate_summary$Total_Count), by = 1)) +  # Set y-axis increments to 1
  coord_flip() +
  theme_minimal() +
  labs(title = "Number of Samples in Each Climate Class Captured by Core",
       x = "Climate Class",
       y = "Number of Samples") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size as needed
p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Climate_niche_covered_by_core.pdf", plot = p, width = 10, height = 8)

########################
#Figure S7C
########################
#Boxplot on side for how groupings cover cliamte niches above 5
########################
# Load the data
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups.csv", row.names = 1)

# Create the plot with larger text
p <- ggplot(data, aes(x = pop, fill = Phylogenetic_tree)) +
  geom_bar() +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  theme_minimal() +
  labs(title = "Samples in Each Climate Niche by Phylogenetic Tree Status",
       x = "Climate Niche (pop)",
       y = "Number of Samples",
       fill = "Phylogenetic Tree Status") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),        # Increase title size
    axis.title.x = element_text(size = 14),                     # Increase x-axis title size
    axis.title.y = element_text(size = 14),                     # Increase y-axis title size
    axis.text.x = element_text(size = 12),                      # Increase x-axis text size
    axis.text.y = element_text(size = 12),                      # Increase y-axis text size
    legend.title = element_text(size = 14),                     # Increase legend title size
    legend.text = element_text(size = 12)                       # Increase legend text size
  ) +
  # Manually define the colors for each Phylogenetic_tree status
  scale_fill_manual(values = c(
    "basal" = "orange",
    "drug-type feral" = "blue",
    "hemp-type" = "green",
    "Population_1" = "#9370DB",
    "Population_2" = "#008080"
  ))

# Display and save the plot
p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Climate_niche_covered_by_dataset_partitions.pdf", plot = p, width = 10, height = 8)


########################
#Figure S7D
########################
#Boxplot on side for how groupings cover cliamte niches above 5 by Type I/II/III
########################
library(ggplot2)
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups.csv", row.names = 1)

# Create the plot with larger text
p<- ggplot(data, aes(x = pop, fill = Type)) +
  geom_bar() +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  theme_minimal() +
  labs(title = "Samples in Each Climate Niche by Type",
       x = "Climate Niche (pop)",
       y = "Number of Samples",
       fill = "Type") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),         # Increase title size
    axis.title.x = element_text(size = 14),                      # Increase x-axis title size
    axis.title.y = element_text(size = 14),                      # Increase y-axis title size
    axis.text.x = element_text(size = 12),                       # Increase x-axis text size
    axis.text.y = element_text(size = 12),                       # Increase y-axis text size
    legend.title = element_text(size = 14),                      # Increase legend title size
    legend.text = element_text(size = 12)                        # Increase legend text size
  ) +
  # Manually define the colors for each 'Type'
  scale_fill_manual(values = c(
    "unknown" = "grey",
    "Type I" = "red",
    "Type II" = "orange",
    "Type III" = "green"
  ))
ggsave("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Climate_niche_covered_by_dataset_partitions_chemistry.pdf", plot = p, width = 10, height = 8)


#########################
#Supplemental Figure 8
##########################
#Cross validation (WorldClim 19) for Ren dataset
#ON HPC
library(rrBLUP) #rrBLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) #hibayes package for Bayesian models
library(dplyr)

#load functions 
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



#########################
#Supplemental Figure 9
##########################
#Cross validation (WorldClim 19) for Ren + Soorni dataset
#HPC
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/xval_kfold_functions.R")

# Read genotype data in rrBLUP format
data <- read.delim("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(data)[1] <- "rs." #needed to rename the sample column here
#colnames(data) <- gsub("^X", "", colnames(data)) #needed to remove the extra X added in 


# Read environmental data from CSV
envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/ren_soorni_n111_extracted_climate_data.csv', head = TRUE)  # full environmental dataset

# Read training dataset from CSV (only core environmental data)
trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/Ren_Soorni_n50_greedy_climate_data.csv", head = TRUE)

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

#####
#PLOTS
#####
library(ggplot2)
setwd("~/R/cannabis_GEAV/Ren_Soorni/cross_validation_19/")

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
ggsave("~/R/cannabis_GEAV/Ren_Soorni/cross_validation_19/cross_val_ren_soorni.pdf", plot = p, device = "pdf", width = 18, height = 10)
rrblup_results <- all_bio[all_bio$model == "rrBLUP", ]

# Export the results to a CSV file
write.csv(rrblup_results, "~/R/cannabis_GEAV/Ren_Soorni/cross_validation_19/rrblup_prediction_accuracy_ren.csv", row.names = FALSE)


#########################
#Supplemental Figure 10
##########################
#Predictions for (WorldClim 19) for Ren (phylogenetic groupings)

#REN ALONE
#(Temp and ppt)

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


#########################
#Supplemental Figure 11
##########################
#Predictions for (WorldClim 19) for Ren (cultivated status groupings)
#REN ALONE

#(Temp and ppt)
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


#########################
#Supplemental Figure 12
##########################
#Predictions for (WorldClim 19) for Ren + Soorni (PCA groupings)

###################
# Step 5: Genomic Selection with RR-BLUP
###################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(gd1)[1] <- "rs."

#envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n50_extracted_climate_data.csv', head = TRUE)
#envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_n50_even_extracted_climate_data.csv', head = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_n50_greedy_climate_data.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

#row.names(g.in) <- gsub("^X", "", row.names(g.in)) #accidental X added
#print(g.in [1:10, 1:10])  # Visualize 

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
pred <- setdiff(row.names(g.in), train)  # This gives # individuals not in the training set
g.pred <- g.in[pred, ]

##################################################
# Calculate the overall percentage of missing values in the entire matrix
total_missing <- sum(is.na(g.train))
total_values <- prod(dim(g.train))
percentage_missing <- (total_missing / total_values) * 100

cat("Overall percentage of missing data:", percentage_missing, "%\n")
#4.201%
##################################################

common_snps <- intersect(colnames(g.train), colnames(g.pred))

# Subset both g.train and g.pred to keep only the common SNPs
g.train <- g.train[, common_snps]
g.pred <- g.pred[, common_snps]

# Impute missing values in g.train and g.pred
# Impute missing values with the column mean
g.train[is.na(g.train)] <- apply(g.train, 2, function(x) mean(x, na.rm = TRUE)) #warning is ok
g.pred[is.na(g.pred)] <- apply(g.pred, 2, function(x) mean(x, na.rm = TRUE))#warning is ok

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
#write.csv(gebv_df, '~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI.csv')
#write.csv(gebv_df, '~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_even.csv')
write.csv(gebv_df, '~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy.csv')

###################
#Visualise results
###################
data2 <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy.csv')
library(gridExtra)
library(ggplot2)

# Reorder the factor levels for Phylogenetic_tree in the desired order
data2$Phylogenetic_tree <- factor(data2$Phylogenetic_tree, levels = c("basal", "hemp-type", "drug-type feral", "drug-type", "Population_1", "Population_2"))

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
bio_vars <- sprintf("bio_%02d", 1:19)
# Loop through each bio variable to generate the plots for temperature variables
for (i in 1:length(temperature_associations)) {
  p <- ggplot(data2, aes_string(x = "Phylogenetic_tree", y = bio_vars[i], fill = "Phylogenetic_tree")) +
    geom_boxplot(color = "black") +
    labs(title = temperature_associations[i], x = NULL, y = bio_vars[i]) +
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue", "Population_1" = "#9370DB","Population_2" = "#008080")) +
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
    scale_fill_manual(values = c("basal" = "orange", "hemp-type" = "green", "drug-type" = "red", "drug-type feral" = "blue", "Population_1" = "#9370DB","Population_2" = "#008080")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 4)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5, linetype = "longdash")
  
  # Add to precipitation plots
  precipitation_plots[[i]] <- p
}

# Save combined temperature plots into one PDF
pdf("~/R/cannabis_GEAV/Ren_Soorni/outputs/GEAV_combined_Soorni_REN_temperature_plots_core_greedy.pdf", width = 8, height = 10)
grid.arrange(grobs = temperature_plots, ncol = 3)  # Arrange the temperature plots in a grid
dev.off()

# Save combined precipitation plots into another PDF
pdf("~/R/cannabis_GEAV/Ren_Soorni/outputs/GEAV_combined_Soorni_REN_precipitation_plots_core_greedy.pdf", width = 8, height = 10)
grid.arrange(grobs = precipitation_plots, ncol = 3)  # Arrange the precipitation plots in a grid
dev.off()

#########################
#Supplemental Figure 13
##########################
#Latitude EGS

###################
# Step 5: Genomic Selection with RR-BLUP
###################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(gd1)[1] <- "rs."

#envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n50_extracted_climate_data.csv', head = TRUE)
#envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_n50_even_extracted_climate_data.csv', head = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_n50_greedy_climate_data.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

#row.names(g.in) <- gsub("^X", "", row.names(g.in)) #accidental X added
#print(g.in [1:10, 1:10])  # Visualize 

# Set row names for environmental data using 'Sample_ID'
row.names(envdat) <- envdat$Sample_ID

# Subset environmental data for analysis (start from the 5th column)
#y.in.rr <- envdat[, 5:ncol(envdat)]  
y.in.rr <- envdat %>% select(Latitude, bio_01, bio_02)
minicore_entries <- which(envdat$Core == TRUE)  # Subset to core lines where 'Core' is TRUE
y.in.rr <- y.in.rr[minicore_entries, ]  
y.in.mat <- as.matrix(y.in.rr)
print(y.in.mat [1:3, 1:3])  # Visualize 

dim(y.in.mat)

# Training set
train <- row.names(y.in.mat)  # Names of training lines
g.train <- g.in[train, ]  # Subset rows (individuals) directly

print(g.train[1:10, 1:10])  # Visualize part of g.train for debugging
dim(g.train)

# Subset g.pred to include only individuals in `pred`
pred <- setdiff(row.names(g.in), train)  # This gives # individuals not in the training set
g.pred <- g.in[pred, ]


##################################################
# Calculate the overall percentage of missing values in the entire matrix
total_missing <- sum(is.na(g.train))
total_values <- prod(dim(g.train))
percentage_missing <- (total_missing / total_values) * 100

cat("Overall percentage of missing data:", percentage_missing, "%\n")
#4.201%
##################################################

common_snps <- intersect(colnames(g.train), colnames(g.pred))

# Subset both g.train and g.pred to keep only the common SNPs
g.train <- g.train[, common_snps]
g.pred <- g.pred[, common_snps]

# Impute missing values in g.train and g.pred
# Impute missing values with the column mean
g.train[is.na(g.train)] <- apply(g.train, 2, function(x) mean(x, na.rm = TRUE)) #warning is ok
g.pred[is.na(g.pred)] <- apply(g.pred, 2, function(x) mean(x, na.rm = TRUE))#warning is ok

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
write.csv(gebv_df, '~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy_LAT.csv')

###################
#PLOT
###################
library(ggplot2)
library(dplyr)

gebv_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy_LAT.csv")

# Ensure Latitude_GEBV and Latitude columns are numeric
gebv_data$Latitude_GEBV <- as.numeric(gebv_data$Latitude_GEBV)
gebv_data$Latitude <- as.numeric(gebv_data$Latitude)

# Sort data by Latitude_GEBV
gebv_data <- gebv_data %>% arrange(Latitude_GEBV)

# Add columns to indicate matching latitude and core status
gebv_data <- gebv_data %>%
  mutate(
    Matching_Latitude = ifelse(is.na(Latitude), "No Match", "Match"),
    Core_Status = ifelse(Core_greedy == TRUE & !is.na(Latitude), "Core", "Non-Core")
  )

# Define colors for Phylogenetic_tree based on the standard colors
phylogenetic_colors <- c(
  "basal" = "orange",
  "hemp-type" = "green",
  "drug-type" = "red",
  "drug-type feral" = "blue",
  "Population_1" = "#9370DB",
  "Population_2" = "#008080"
)

# Create a named vector for label colors based on the Phylogenetic_tree column
label_colors <- phylogenetic_colors[gebv_data$Phylogenetic_tree]
names(label_colors) <- gebv_data$Name

ggplot(gebv_data, aes(x = reorder(Name, Latitude_GEBV), y = Latitude_GEBV)) +
  # Plot GEBVs for latitude with color based on matching
  geom_point(aes(color = Matching_Latitude), size = 8) +
  # Overlay actual latitude values with core color distinction
  geom_point(aes(y = Latitude, color = Core_Status), shape = 1, size = 8) +
  coord_flip() +  # Flip coordinates for a horizontal layout
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  # Color scale for Phylogenetic_tree
  scale_fill_manual(
    name = "Phylogenetic Group",
    values = phylogenetic_colors
  ) +
  # Customize color for Matching_Latitude and Core_Status
  scale_color_manual(
    name = "Legend",
    values = c("Match" = "black", "No Match" = "red", "Core" = "red", "Non-Core" = "blue"),
    labels = c("GEBV Match", "GEBV No Match", "Actual Latitude (Core)", "Actual Latitude (Non-Core)")
  ) +
  scale_y_continuous(
    breaks = c(seq(floor(min(gebv_data$Latitude_GEBV)), ceiling(max(gebv_data$Latitude_GEBV)), by = 10), 20, 30, 40, 50, 60),
    minor_breaks = NULL
  ) +
  geom_hline(yintercept = seq(20, 60, by = 10), linetype = "dotted", color = "grey50") +
  theme(
    axis.text.y = element_text(color = label_colors, size = 20),  # Increased y-axis text size
    axis.text.x = element_text(angle = 270, hjust = 1, size = 30),  # Increased x-axis text size
    legend.position = "bottom"
  )

ggsave("~/R/cannabis_GEAV/FIGS/Latitude_GEBV_dot_plot_with_core_overlay.pdf", width = 20, height = 40)


############### THE SUGGESTION FOR PER INDIVIDUAL EFFECTS FOR LATITUDE
# Load required libraries
library(rrBLUP)
library(dplyr)
library(tidyr)

# Load the genotype matrix
genotype_matrix <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_genotype_matrix_no_missing_13670.csv', row.names = 1)

# Ensure row names are correctly set and no extra columns
rownames(genotype_matrix) <- genotype_matrix$X  # Set the X column as row names
genotype_matrix <- genotype_matrix[ , -1]       # Remove the X column

# Load environmental data
envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_n50_greedy_climate_data.csv', head = TRUE)

# Set row names for environmental data
row.names(envdat) <- envdat$Sample_ID

# Define the core samples (n=50)
core_samples <- c("SRR14708197", "SRR14708199", "SRR14708202", "SRR14708203", "SRR14708204", 
                  "SRR14708205", "SRR14708206", "SRR14708207", "SRR14708208", "SRR14708209", 
                  "SRR14708210", "SRR14708212", "SRR14708213", "SRR14708229", "SRR14708230", 
                  "SRR14708231", "SRR14708234", "SRR14708235", "SRR14708236", "SRR14708243", 
                  "SRR14708246", "SRR14708248", "SRR14708251", "SRR14708275", "SRR14708276", 
                  "SRR9089225.1", "SRR9089228.1", "SRR9089230.1", "SRR9089236.1", "SRR9089253.1", 
                  "SRR9089257.1", "SRR9089259.1", "SRR9089261.1", "SRR9089262.1", "SRR9089263.1", 
                  "SRR9089266.1", "SRR9089270.1", "SRR9089275.1", "SRR9089276.1", "SRR9089282.1", 
                  "SRR9089285.1", "SRR9089295.1", "SRR9089297.1", "SRR9089299.1", "SRR9089300.1", 
                  "SRR9089301.1", "SRR9089303.1", "SRR9089304.1", "SRR9089308.1", "SRR9089310.1")

# Filter genotype matrix to include only core samples
genotype_core <- genotype_matrix[rownames(genotype_matrix) %in% core_samples, , drop = FALSE]

# Filter environmental data for core samples and ensure correct ordering
envdat_core <- envdat[envdat$Sample_ID %in% core_samples, ]
envdat_core <- envdat_core[match(rownames(genotype_core), envdat_core$Sample_ID), ]

# Check alignment
if (!all(rownames(genotype_core) == envdat_core$Sample_ID)) {
  stop("Sample order mismatch still exists. Please check data.")
} else {
  cat("Sample order in envdat_core now matches genotype_core.\n")
}

# Convert genotype_core to matrix
genotype_core <- as.matrix(genotype_core)

# Run RR-BLUP on the core samples for latitude
solve_lat <- mixed.solve(y = envdat_core$Latitude, Z = genotype_core)
u.hat_lat <- solve_lat$u  # SNP effects for Latitude

# Apply SNP effects to all individuals in complete genotype matrix
complete_markers <- as.matrix(genotype_matrix)

# Calculate SNP-level contributions for each individual for Latitude
individual_marker_effects_lat <- complete_markers * matrix(u.hat_lat, nrow = nrow(complete_markers), ncol = ncol(complete_markers), byrow = TRUE)

# Convert to data frame for easier interpretation
individual_marker_effects_lat_df <- as.data.frame(individual_marker_effects_lat)
colnames(individual_marker_effects_lat_df) <- colnames(complete_markers)
individual_marker_effects_lat_df$Individual <- rownames(complete_markers)
individual_marker_effects_lat_df$Trait <- "Latitude"

# Reshape to long format for per-SNP effect per individual
individual_marker_effects_lat_long <- individual_marker_effects_lat_df %>%
  pivot_longer(
    cols = -c(Individual, Trait),
    names_to = "SNP",
    values_to = "Effect"
  )

# Save the per-SNP, per-individual effect data for Latitude to CSV
write.csv(individual_marker_effects_lat_long, '~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_latitude_per_SNP_NEWNEW.csv', row.names = FALSE)

# Print a summary to check the output
head(individual_marker_effects_lat_long)



#########################
#Supplemental Figure 14
##########################
#RDA
###############
library(vcfR)
library(adegenet)
vcf <- read.vcfR("~/R/cannabis_GEAV/MERGE2/merged_Ren_Soorni_snps.vcf.gz")
genotype_matrix <- extract.gt(vcf, element = "GT")


# Function to convert genotypes to numeric format
convert_geno <- function(x) {
  ifelse(x == "0/0", 0, ifelse(x == "0/1" | x == "1/0", 1, ifelse(x == "1/1", 2, NA)))
}

# Apply this function to each column of the genotype matrix
geno_numeric <- apply(genotype_matrix, 2, convert_geno)


geno_numeric_t <- t(geno_numeric)
geno_df <- as.data.frame(geno_numeric_t)


rownames(geno_df) <- colnames(genotype_matrix) 

snp_names <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")
colnames(geno_df) <- snp_names

# Write the genotype matrix to a CSV file
write.csv(geno_df, "~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv")

############################# WORKS
library(vegan)
library(adegenet)
library(psych)

# Load your genotype and environment matrices
genotype_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv", row.names = 1)
env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", row.names = 1)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

# Find common samples (individuals) between genotype and environmental data
common_samples <- intersect(rownames(genotype_data), rownames(env_data_clean))

# Subset both genotype and environmental data to keep only the common samples
genotype_data_clean <- genotype_data[common_samples, ]
env_data_clean <- env_data_clean[common_samples, ]

# Check if the number of rows now match
cat("Number of rows in genotype data:", nrow(genotype_data_clean), "\n")
cat("Number of rows in environmental data:", nrow(env_data_clean), "\n")

# Check for missing or non-finite values in genotype data
cat("Checking for remaining NA values in genotype data...\n")
missing_genotypes <- sum(is.na(genotype_data_clean))
cat("Missing values in genotype data: ", missing_genotypes, "\n")

cat("Checking for remaining non-finite values in genotype data...\n")
non_finite_genotypes <- sum(!is.finite(as.matrix(genotype_data_clean)))
cat("Non-finite values in genotype data: ", non_finite_genotypes, "\n")

# Impute missing values in genotype data by filling with the most common value (0, 1, or 2)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with most common value
  return(x)
}

# Apply the imputation function to each SNP column in the genotype data
cat("Imputing missing values in genotype data...\n")
genotype_data_clean <- apply(genotype_data_clean, 2, impute_geno)

# Check the dimensions after imputation
cat("Final genotype data dimensions after imputing missing values: ", dim(genotype_data_clean), "\n")

# Ensure the common samples still match after cleaning
common_samples <- intersect(rownames(genotype_data_clean), rownames(env_data_clean))
genotype_data_clean <- genotype_data_clean[common_samples, ]
env_data_clean <- env_data_clean[common_samples, ]

# Perform RDA with the cleaned and imputed/filtered data
cat("Performing RDA...\n")
rda_result <- rda(genotype_data_clean ~ ., data = env_data_clean)

# View the summary of RDA results
cat("Summary of RDA results:\n")
print(summary(rda_result)) #captures the variance explained by each axis

# Perform a permutation test for significance
cat("Performing permutation test...\n")
rda_test <- anova(rda_result, permutations = 999) #slow step get a cup of tea
print(rda_test)

# Plot the RDA results
cat("Plotting RDA results...\n")
plot(rda_result, scaling = 3)

pdf("RDA_plot_customized.pdf", width = 18, height = 14)
plot(rda_result, scaling = 3, main = "RDA for Latitude, Longitude and 19 WorldClim variables", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
dev.off()

library(vegan)
RsquareAdj(rda_result)

############################################
#Effect of population Structure in the RDA
############################################
genotype_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv", row.names = 1)
population_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy_trimmed_for_climate.csv")
env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", row.names = 1)
env_data_clean <- na.omit(env_data)

# Subset both genotype and environmental data to keep only the common samples
common_samples <- intersect(rownames(genotype_data), rownames(env_data_clean))
env_data_clean <- env_data_clean[!is.na(env_data_clean$Phylogenetic_tree), ]

genotype_data_clean <- genotype_data[common_samples, ]
env_data_clean <- env_data_clean[common_samples, ]

env_data$Sample_ID <- rownames(env_data) 
env_data_clean <- merge(env_data, population_data, by = "Sample_ID", all.x = TRUE)

if ("Phylogenetic_tree" %in% colnames(env_data_clean)) {
  cat("Phylogenetic_tree column found in env_data_clean.\n")
} else {
  cat("Phylogenetic_tree column not found. Check population_data and merge key.\n")
}

structure_var <- as.factor(env_data_clean$Phylogenetic_tree)
unique(structure_var)


# Run the RDA model with structure as a predictor
rda_structure <- rda(genotype_data_clean ~ structure_var, data = env_data_clean)
summary(rda_structure)

# Permutation test for significance of the RDA model
rda_structure_test <- anova(rda_structure, permutations = 999)
print(rda_structure_test)

# Calculate adjusted R-squared
RsquareAdj(rda_structure)

###WORKING FOR STRUCTURE
# Load data
genotype_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv", row.names = 1)
population_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy_trimmed_for_climate.csv")
env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", row.names = 1)

# Add Sample_ID as a column in env_data
env_data$Sample_ID <- rownames(env_data)

# Merge environmental data with population data
env_data_clean <- merge(env_data, population_data, by = "Sample_ID", all.x = TRUE)

# Remove rows with NA in 'Phylogenetic_tree' (structure variable)
env_data_clean <- env_data_clean[!is.na(env_data_clean$Phylogenetic_tree), ]

# Ensure common samples are used in both datasets
common_samples <- intersect(rownames(genotype_data), env_data_clean$Sample_ID)
genotype_data_clean <- genotype_data[common_samples, ]
env_data_clean <- env_data_clean[env_data_clean$Sample_ID %in% common_samples, ]


# Ensure structure_var is properly defined as a factor and check for NA
structure_var <- as.factor(env_data_clean$Phylogenetic_tree)
if (any(is.na(structure_var))) {
  stop("NA values remain in structure_var after filtering. Please check the input data.")
}

# Check for and handle NA values in genotype_data_clean
if (any(is.na(genotype_data_clean))) {
  cat("NA values detected in genotype_data_clean; imputing missing values...\n")
  # Impute NA values in each column with the column mean
  for (col in 1:ncol(genotype_data_clean)) {
    genotype_data_clean[, col][is.na(genotype_data_clean[, col])] <- mean(genotype_data_clean[, col], na.rm = TRUE)
  }
}

# Run the RDA model with structure as a predictor
rda_structure <- rda(genotype_data_clean ~ structure_var, data = env_data_clean)
summary(rda_structure)

# Perform a permutation test
rda_structure_test <- anova(rda_structure, permutations = 999)
print(rda_structure_test)

# Calculate adjusted R-squared
RsquareAdj(rda_structure)


############
#separating out variables
############
library(vegan)
library(adegenet)
library(psych)

# Load your genotype and environment matrices
genotype_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/genotype_matrix.csv", row.names = 1)
env_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/GEA/GEA_environment_input.csv", row.names = 1)

# Remove rows with NA values from environmental data
env_data_clean <- na.omit(env_data)

# Find common samples (individuals) between genotype and environmental data
common_samples <- intersect(rownames(genotype_data), rownames(env_data_clean))

# Subset both genotype and environmental data to keep only the common samples
genotype_data_clean <- genotype_data[common_samples, ]
env_data_clean <- env_data_clean[common_samples, ]

# Impute missing values in genotype data by filling with the most common value (0, 1, or 2)
impute_geno <- function(x) {
  x[is.na(x)] <- as.numeric(names(which.max(table(x, useNA = "no"))))  # Replace NA with most common value
  return(x)
}

# Apply the imputation function to each SNP column in the genotype data
genotype_data_clean <- apply(genotype_data_clean, 2, impute_geno)

# Subset environmental data into temperature, precipitation, and location variables
temperature_vars <- env_data_clean[, c("bio_01","bio_02","bio_03","bio_04","bio_05","bio_06","bio_07", "bio_08", "bio_09", "bio_10", "bio_11")]  #temperature
precipitation_vars <- env_data_clean[, c("bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")]  #precipitation
location_vars <- env_data_clean[, c("Latitude", "Longitude")] #geographic
location_vars_2 <- data.frame(Latitude = env_data_clean$Latitude) #latitude alone

# Perform RDA with Latitude as the explanatory variable


# Perform RDA for temperature variables
cat("Performing RDA for temperature variables...\n")
rda_temp <- rda(genotype_data_clean ~ ., data = temperature_vars)
print(summary(rda_temp)) # View the summary of RDA results for temperature
rda_temp_test <- anova(rda_temp, permutations = 999) # Permutation test - slow step
print(rda_temp_test)
#plot(rda_temp, scaling = 3, main = "RDA for Temperature Variables (Variance explained 15.02%)")
text(rda_temp, display = "bp", cex = 2, col = "blue")
plot(rda_temp, scaling = 3, main = "RDA for Temperature Variables (Variance explained 15.02%)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

pdf("RDA_plot_temp.pdf", width = 18, height = 14)
plot(rda_temp, scaling = 3, main = "RDA for Temperature Variables (BIO1-11)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
dev.off()

library(vegan)
RsquareAdj(rda_temp)


# Perform RDA for precipitation variables
cat("Performing RDA for precipitation variables...\n")
rda_precip <- rda(genotype_data_clean ~ ., data = precipitation_vars)
print(summary(rda_precip)) # View the summary of RDA results for precipitation
rda_precip_test <- anova(rda_precip, permutations = 999) # Permutation test- slow step
print(rda_precip_test)
#text(rda_precip, display = "bp", cex = 2, col = "blue")
#plot(rda_precip, scaling = 3, main = "RDA for Precipitation Variables (Variance explained 12.14%)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

pdf("RDA_plot_ppt.pdf", width = 18, height = 14)
plot(rda_precip, scaling = 3, main = "RDA for Precipitation Variables (BIO12-19)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
dev.off()


library(vegan)
RsquareAdj(rda_precip)

# Perform RDA for location variables (latitude and longitude)
cat("Performing RDA for location variables...\n")
rda_location <- rda(genotype_data_clean ~ ., data = location_vars)
print(summary(rda_location)) # View the summary of RDA results for location
rda_location_test <- anova(rda_location, permutations = 999) # Permutation test- slow step
print(rda_location_test)
text(rda_location, display = "bp", cex = 2, col = "blue")
plot(rda_location, scaling = 3, main = "RDA for Latitude and Longitude (Variance explained 6.35%)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

library(vegan)
RsquareAdj(rda_location)

# Perform RDA for location variables (latitude)
cat("Performing RDA for location variables...\n")
rda_location_2 <- rda(genotype_data_clean ~ Latitude, data = location_vars_2)
print(summary(rda_location_2))
rda_location_test_2 <- anova(rda_location_2, permutations = 999) # Permutation test- slow step
print(rda_location_test_2)
#text(rda_location_2, display = "bp", cex = 2, col = "blue")
#plot(rda_location_2, scaling = 3, main = "RDA for Latitude (Variance explained 4.84%)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)

pdf("RDA_plot_lat.pdf", width = 18, height = 14)
plot(rda_location_2, scaling = 3, main = "RDA for Latitude", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
dev.off()


library(vegan)
RsquareAdj(rda_location_2)



#########################
#Supplemental Figure 15
##########################
#Figure S15A
########################
#Boxplot on side for how groupings cover cliamte niches above 5
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups.csv", row.names = 1)

p <- ggplot(data, aes(x = pop, fill = Phylogenetic_tree)) +
  geom_bar() +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  theme_minimal() +
  labs(title = "",
       x = "Climate Niche",
       y = "Number of Samples",
       fill = "Group") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),        # Increase title size
    axis.title.x = element_text(size = 14),                     # Increase x-axis title size
    axis.title.y = element_text(size = 14),                     # Increase y-axis title size
    axis.text.x = element_text(size = 12),                      # Increase x-axis text size
    axis.text.y = element_text(size = 12),                      # Increase y-axis text size
    legend.title = element_text(size = 14),                     # Increase legend title size
    legend.text = element_text(size = 12)                       # Increase legend text size
  ) +
  scale_fill_manual(values = c(
    "basal" = "orange",
    "drug-type feral" = "blue",
    "hemp-type" = "green",
    "Population_1" = "#9370DB",
    "Population_2" = "#008080"
  ), labels = c("Basal", "Drug-type feral", "Hemp-type", "Population-1", "Population-2"))

# Display and save the plot
p
ggsave("~/R/cannabis_GEAV/FIGS/15A.pdf", plot = p, width = 10, height = 8)

#############
#Figure S15B
#############
#groups with no other population level admixture
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/Climates_overlap_with_groups_5plus.csv", row.names = 1)

p <- ggplot(data, aes(x = pop, fill = Phylogenetic_tree)) +
  geom_bar() +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  theme_minimal() +
  labs(title = "",
       x = "Climate Niche",
       y = "Number of Samples",
       fill = "Group") +
  scale_fill_manual(values = c(
    "basal" = "orange",
    "drug-type feral" = "blue",
    "hemp-type" = "green",
    "Population_1" = "#9370DB",
    "Population_2" = "#008080"
  ), labels = c("Basal", "Drug-type feral", "Hemp-type", "Population-1", "Population-2")) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12),  # Adjust y-axis label size
    axis.text.x = element_text(size = 12),                        # Adjust x-axis label size
    axis.title.x = element_text(size = 14),                       # Adjust x-axis title size
    axis.title.y = element_text(size = 14),                       # Adjust y-axis title size
    plot.title = element_text(size = 16, face = "bold"),          # Adjust title size
    legend.title = element_text(size = 14),                       # Adjust legend title size
    legend.text = element_text(size = 12)                         # Adjust legend text size
  )

# Save the plot
p
ggsave("~/R/cannabis_GEAV/FIGS/S15B.pdf", plot = p, width = 10, height = 8)

#############
#Figure S15C
#############
#pi
#use vcftools on merged vcf files using a population .txt file 
#bcftools view -S BSh_population_file.txt -Oz -o BSh_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S BWh_population_file.txt -Oz -o BWh_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S BWk_population_file.txt -Oz -o BWk_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S BSk_population_file.txt -Oz -o BSk_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S DSa_population_file.txt -Oz -o DSa_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S Cwa_population_file.txt -Oz -o Cwa_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz
#bcftools view -S Dfb_population_file.txt -Oz -o Dfb_population.vcf.gz merged_Ren_Soorni_snps.vcf.gz

#getting separate pi for each new vcf file
#vcftools --gzvcf BSh_population.vcf.gz --site-pi --out BSh_diversity
#vcftools --gzvcf BWh_population.vcf.gz --site-pi --out BWh_diversity
#vcftools --gzvcf BWk_population.vcf.gz --site-pi --out BWk_diversity
#vcftools --gzvcf BSk_population.vcf.gz --site-pi --out BSk_diversity
#vcftools --gzvcf DSa_population.vcf.gz --site-pi --out DSa_diversity
#vcftools --gzvcf Cwa_population.vcf.gz --site-pi --out Cwa_diversity
#vcftools --gzvcf Dfb_population.vcf.gz --site-pi --out Dfb_diversity

library(ggplot2)
library(dplyr)

# Load the populations
bsh_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/BSh_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "BSh: Population 1, n=6") # Replace XX with the actual sample size
bwh_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/BWh_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "BWh: Population 1, n=9")
bwk_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/BWk_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "BWk: Basal, n=7")
bsk_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/BSk_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "BSk: Population 1, n=20")
dsa_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Dsa_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Dsa: Population 2, n=5")
cwa_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Cwa_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Cwa: Drug-type-feral, n=6")
dfb_pi <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Dfb_diversity.sites.pi", header = TRUE) %>% 
  mutate(Population = "Dfb: Hemp, n=10") # Include n=10 as an example

# Combine all population data into a single dataframe
pi_data <- bind_rows(bsh_pi, bwh_pi, bwk_pi, bsk_pi, dsa_pi, cwa_pi, dfb_pi)

# Define chromosome mapping with both number and NC identifier
chromosome_mapping <- c(
  "NC_044371.1" = "1 - NC_044371.1",
  "NC_044375.1" = "2 - NC_044375.1",
  "NC_044372.1" = "3 - NC_044372.1",
  "NC_044373.1" = "4 - NC_044373.1",
  "NC_044374.1" = "5 - NC_044374.1",
  "NC_044377.1" = "6 - NC_044377.1",
  "NC_044378.1" = "7 - NC_044378.1",
  "NC_044379.1" = "8 - NC_044379.1",
  "NC_044376.1" = "9 - NC_044376.1",
  "NC_044370.1" = "10 - NC_044370.1"
)

# Add Chromosome column with labels and ordered levels
pi_data <- pi_data %>%
  mutate(Chromosome = factor(chromosome_mapping[CHROM], 
                             levels = c("1 - NC_044371.1", "2 - NC_044375.1", "3 - NC_044372.1", 
                                        "4 - NC_044373.1", "5 - NC_044374.1", "6 - NC_044377.1", 
                                        "7 - NC_044378.1", "8 - NC_044379.1", "9 - NC_044376.1", 
                                        "10 - NC_044370.1")))

# Plot by chromosome with combined labels
# Define custom colors for each population to improve visibility
custom_colors <- c("BSh: Population 1, n=6" = "purple",
                   "BWh: Population 1, n=9" = "purple4",
                   "BWk: Basal, n=7" = "orange",
                   "BSk: Population 1, n=20" = "#e78ac3",
                   "Dsa: Population 2, n=5" = "#66c2a5",
                   "Cwa: Drug-type-feral, n=6" = "blue",  # Use a darker yellow for better visibility
                   "Dfb: Hemp, n=10" = "green")

p <- ggplot(pi_data, aes(x = POS, y = PI, color = Population)) +
  geom_smooth(alpha = 0.5, size = 0.8, se = FALSE) +
  labs(title = "",
       x = "Genomic Position (Mb)",
       y = "Nucleotide Diversity (pi)") +
  facet_wrap(~ Chromosome, ncol = 2, nrow = 5, scales = "free_x") +  # Set two columns and five rows
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(0, max(pi_data$POS), by = 1e7),  # Set tick marks every 10 Mb
                     labels = function(x) paste0(x / 1e6, " Mb")) +  # Convert to Mb for readability
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),  # Increase font size for x-axis text
    axis.text.y = element_text(size = 20),  # Increase font size for y-axis text
    axis.title.x = element_text(size = 20, face = "plain"),  # Increase font size for x-axis title
    axis.title.y = element_text(size = 20, face = "plain"),  # Increase font size for y-axis title
    plot.title = element_text(size = 20, face = "plain", hjust = 0.5),  # Increase font size for plot title
    strip.text = element_text(size = 20, face = "plain"),  # Increase font size for facet labels
    legend.text = element_text(size = 20),  # Increase font size for legend text
    legend.title = element_text(size = 20, face = "plain")  # Increase font size for legend title
  ) +
  theme(panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank())  # Remove minor grid lines

p

ggsave("~/R/cannabis_GEAV/FIGS/Climate_niche_pi_2.pdf", plot = p, width = 20, height = 22)

#########################
#Supplemental Figure 16
##########################
#LD plots

################
#Inversion LD
################

#Prior steps locally using bcftools
#bcftools index Dsa_population.vcf.gz
#bcftools index Dfb_population.vcf.gz

#bcftools view -r NC_044373.1:70000000-90000000 -o Dsa_population_chr4_region.vcf Dsa_population.vcf.gz
#bcftools view -r NC_044373.1:70000000-90000000 -o Dfb_population_chr4_region.vcf Dfb_population.vcf.gz

#./plink --vcf Dsa_population_chr4_region.vcf --make-bed --out Dsa_chr4_region --allow-extra-chr
#./plink --vcf Dfb_population_chr4_region.vcf --make-bed --out Dfb_chr4_region --allow-extra-chr

#calculate LD
#./plink --bfile Dsa_chr4_region --r2 --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0.8 --out Dsa_chr4_ld --allow-extra-chr
#./plink --bfile Dfb_chr4_region --r2 --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0.8 --out Dfb_chr4_ld --allow-extra-chr

library(ggplot2)
dsa_ld <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Dsa_chr4_ld.ld", header = TRUE)
dfb_ld <- read.table("~/R/cannabis_GEAV/Ren_Soorni/pi_vcftools/Dfb_chr4_ld.ld", header = TRUE)

# Plot for DSa Population 2
p <- ggplot(dsa_ld, aes(x = BP_A / 1e6, y = BP_B / 1e6, color = R2)) +  # Convert positions to Mb
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(min(dsa_ld$BP_A) / 1e6, max(dsa_ld$BP_A) / 1e6),  # Adjust x-axis limits
                     labels = scales::number_format(suffix = " Mb")) +            # Add "Mb" suffix
  scale_y_continuous(limits = c(min(dsa_ld$BP_B) / 1e6, max(dsa_ld$BP_B) / 1e6),  # Adjust y-axis limits
                     labels = scales::number_format(suffix = " Mb")) +            # Add "Mb" suffix
  labs(title = "Linkage Disequilibrium in Chr 4 (DSa Population 2)",
       x = "Position A (Mb)",
       y = "Position B (Mb)",
       color = "R2") +
  theme_minimal()
ggsave("~/R/cannabis_GEAV/FIGS/Dsa_chr4.pdf", plot = p, width = 10, height = 10)

# Plot for Dfb Hemp Population
p2 <- ggplot(dfb_ld, aes(x = BP_A / 1e6, y = BP_B / 1e6, color = R2)) +  # Convert positions to Mb
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(min(dfb_ld$BP_A) / 1e6, max(dfb_ld$BP_A) / 1e6),  # Adjust x-axis limits
                     labels = scales::number_format(suffix = " Mb")) +            # Add "Mb" suffix
  scale_y_continuous(limits = c(min(dfb_ld$BP_B) / 1e6, max(dfb_ld$BP_B) / 1e6),  # Adjust y-axis limits
                     labels = scales::number_format(suffix = " Mb")) +            # Add "Mb" suffix
  labs(title = "Linkage Disequilibrium in Chr 4 (Dfb Hemp Population)",
       x = "Position A (Mb)",
       y = "Position B (Mb)",
       color = "R2") +
  theme_minimal()
ggsave("~/R/cannabis_GEAV/FIGS/Dfb_chr4.pdf", plot = p2, width = 10, height = 10)

########################
#Supplemental Figure 17
########################
#trimming for clean samples per environmental niche 
#wrapping GEBVs per climate niche for subset to see how GEBV driven by climate
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Define the custom color palette
climate_colors <- c(
  "BWh - Arid, desert, hot" = "#482576FF", 
  "BWk - Arid, desert, cold" = "#46307EFF",  
  "BSh - Arid, steppe, hot" = "#443B84FF",  
  "BSk - Arid, steppe, cold" = "#404688FF",  
  "Csa - Temperate, dry summer, hot summer" = "#3C508BFF",  
  "Cwa - Temperate, dry winter, hot summer" = "#2F6B8EFF",  
  "Dsa - Cold, dry summer, hot summer" = "#1E9D89FF",  
  "Dfb - Cold, no dry season, warm summer" = "#A2DA37FF"
)

# Read the data
data <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/inputs/GEBVs_for_n111_core_greedy.csv')

# List of sample IDs to keep
sample_ids <- c("SRR9089275.1", "SRR9089276.1", "SRR9089281.1", "SRR9089282.1", "SRR9089285.1", "SRR9089286.1", 
                "SRR9089230.1", "SRR9089255.1", "SRR9089256.1", "SRR9089257.1", "SRR9089258.1", "SRR9089263.1", 
                "SRR9089264.1", "SRR9089269.1", "SRR9089270.1", "SRR9089277.1", "SRR9089278.1", "SRR9089283.1", 
                "SRR9089284.1", "SRR9089289.1", "SRR9089290.1", "SRR9089296.1", "SRR9089301.1", "SRR9089308.1", 
                "SRR9089309.1", "SRR9089310.1", "SRR9089225.1", "SRR9089267.1", "SRR9089268.1", "SRR9089287.1", 
                "SRR9089288.1", "SRR9089299.1", "SRR9089303.1", "SRR9089304.1", "SRR9089317.1", "SRR14708203", 
                "SRR14708204", "SRR14708205", "SRR14708206", "SRR14708207", "SRR14708208", "SRR14708209", 
                "SRR14708199", "SRR14708198", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258", 
                "SRR14708213", "SRR14708217", "SRR14708230", "SRR14708231", "SRR14708232", "SRR14708244", 
                "SRR14708259", "SRR14708273", "SRR14708275", "SRR14708278", "SRR9089228.1", "SRR9089229.1", 
                "SRR9089291.1", "SRR9089305.1", "SRR9089306.1")

# Filter the data to keep only specified sample IDs
filtered_data <- data %>% filter(Sample_ID %in% sample_ids)

# Further filter down to specific climate classes
filtered_data <- filtered_data %>%
  filter(Climate_Class_Description %in% c("BSk - Arid, steppe, cold", 
                                          "BWk - Arid, desert, cold", 
                                          "Cwa - Temperate, dry winter, hot summer", 
                                          "Dfb - Cold, no dry season, warm summer", 
                                          "BSh - Arid, steppe, hot", 
                                          "BWh - Arid, desert, hot", 
                                          "Csa - Temperate, dry summer, hot summer", 
                                          "Dsa - Cold, dry summer, hot summer"))

# Reshape the data from wide to long format
data_long <- melt(filtered_data, id.vars = c("Sample_ID", "Climate_Class_Description"),
                  measure.vars = paste0("bio_", sprintf("%02d", 1:19)),
                  variable.name = "Bio_Trait", value.name = "GEBV")

# Separate the temperature variables (bio_01 to bio_11)
temperature_data <- subset(data_long, Bio_Trait %in% paste0("bio_", sprintf("%02d", 1:11)))

# Separate the precipitation variables (bio_12 to bio_19)
precipitation_data <- subset(data_long, Bio_Trait %in% paste0("bio_", sprintf("%02d", 12:19)))

# Plot for Temperature Variables
# Plot for Temperature Variables with increased text size
p1 <- ggplot(temperature_data, aes(x = Bio_Trait, y = GEBV, fill = Climate_Class_Description)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),                        # Increase y-axis text size
    axis.title.x = element_text(size = 18),                       # Increase x-axis title size
    axis.title.y = element_text(size = 18),                       # Increase y-axis title size
    plot.title = element_text(size = 18, face = "bold"),          # Increase plot title size
    legend.text = element_text(size = 18),                        # Increase legend text size
    legend.title = element_text(size = 18)                        # Increase legend title size
  ) +
  labs(
    title = "GEBVs for Temperature Traits (BIO1 to BIO11)",
    x = "Temperature Traits",
    y = "GEBV"
  ) +
  scale_fill_manual(values = climate_colors) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Climate Class Description"))
p1

ggsave("~/R/cannabis_GEAV/FIGS/Climate_wrap_gebv_bio1_to_bio11.pdf", plot = p1, width = 20, height = 8)


p2 <- ggplot(precipitation_data, aes(x = Bio_Trait, y = GEBV, fill = Climate_Class_Description)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),  # Increase x-axis text size
    axis.text.y = element_text(size = 18),                        # Increase y-axis text size
    axis.title.x = element_text(size = 18),                       # Increase x-axis title size
    axis.title.y = element_text(size = 18),                       # Increase y-axis title size
    plot.title = element_text(size = 18, face = "bold"),          # Increase plot title size
    legend.text = element_text(size = 18),                        # Increase legend text size
    legend.title = element_text(size = 18)                        # Increase legend title size
  ) +
  labs(
    title = "GEBVs for Precipitation Traits (BIO12 to BIO19)",
    x = "Precipitation Traits",
    y = "GEBV"
  ) +
  scale_fill_manual(values = climate_colors) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Climate Class Description"))
p2

ggsave("~/R/cannabis_GEAV/FIGS/Climate_wrap_gebv_bio11_to_bio19.pdf", plot = p2, width = 20, height = 8)



#########################
#Supplemental Figure 18
##########################
#Cross validation (WorldClim 8X) for Ren dataset
#on HPC

#running in batches of 12
#####
#PLOTS 
#####
#PREC
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/prec/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_prec_01', 'wc2.1_30s_prec_02', 'wc2.1_30s_prec_03', 
                         'wc2.1_30s_prec_04', 'wc2.1_30s_prec_05', 'wc2.1_30s_prec_06', 
                         'wc2.1_30s_prec_07', 'wc2.1_30s_prec_08', 'wc2.1_30s_prec_09', 
                         'wc2.1_30s_prec_10', 'wc2.1_30s_prec_11', 'wc2.1_30s_prec_12')

all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_prec",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/prec/prec_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")

##########################
#SRAD
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/srad/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_srad_01', 'wc2.1_30s_srad_02', 'wc2.1_30s_srad_03', 
                         'wc2.1_30s_srad_04', 'wc2.1_30s_srad_05', 'wc2.1_30s_srad_06', 
                         'wc2.1_30s_srad_07', 'wc2.1_30s_srad_08', 'wc2.1_30s_srad_09', 
                         'wc2.1_30s_srad_10', 'wc2.1_30s_srad_11', 'wc2.1_30s_srad_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_srad",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/srad/srad_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")

##########################
# TAVG
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tavg/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_tavg_01', 'wc2.1_30s_tavg_02', 'wc2.1_30s_tavg_03', 
                         'wc2.1_30s_tavg_04', 'wc2.1_30s_tavg_05', 'wc2.1_30s_tavg_06', 
                         'wc2.1_30s_tavg_07', 'wc2.1_30s_tavg_08', 'wc2.1_30s_tavg_09', 
                         'wc2.1_30s_tavg_10', 'wc2.1_30s_tavg_11', 'wc2.1_30s_tavg_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_tavg",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tavg/tavg_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")


##########################
# TMIN
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tmin/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_tmin_01', 'wc2.1_30s_tmin_02', 'wc2.1_30s_tmin_03', 
                         'wc2.1_30s_tmin_04', 'wc2.1_30s_tmin_05', 'wc2.1_30s_tmin_06', 
                         'wc2.1_30s_tmin_07', 'wc2.1_30s_tmin_08', 'wc2.1_30s_tmin_09', 
                         'wc2.1_30s_tmin_10', 'wc2.1_30s_tmin_11', 'wc2.1_30s_tmin_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_tmin",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tmin/tmin_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")

##########################
# TMAX
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tmax/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_tmax_01', 'wc2.1_30s_tmax_02', 'wc2.1_30s_tmax_03', 
                         'wc2.1_30s_tmax_04', 'wc2.1_30s_tmax_05', 'wc2.1_30s_tmax_06', 
                         'wc2.1_30s_tmax_07', 'wc2.1_30s_tmax_08', 'wc2.1_30s_tmax_09', 
                         'wc2.1_30s_tmax_10', 'wc2.1_30s_tmax_11', 'wc2.1_30s_tmax_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_tmax",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/tmax/tmax_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")

##########################
# VAPR
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/vapr/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_vapr_01', 'wc2.1_30s_vapr_02', 'wc2.1_30s_vapr_03', 
                         'wc2.1_30s_vapr_04', 'wc2.1_30s_vapr_05', 'wc2.1_30s_vapr_06', 
                         'wc2.1_30s_vapr_07', 'wc2.1_30s_vapr_08', 'wc2.1_30s_vapr_09', 
                         'wc2.1_30s_vapr_10', 'wc2.1_30s_vapr_11', 'wc2.1_30s_vapr_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_vapr",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/vapr/vapr_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")

##########################
# WIND
library(ggplot2)
setwd("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/wind/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

#Remove any NA values from the model results
model_list1 <- lapply(model_list, na.omit)

#Combine all models into a single dataframe
all_models <- do.call("rbind", model_list1)

#Convert standard deviation values to numeric
all_models$r.sd <- as.numeric(all_models$r.sd)

#Ensure cross-validation type is a factor with specified levels
all_models$xval <- factor(all_models$xval, levels = c("Ten-Fold"))

monthly_prec_traits <- c('wc2.1_30s_wind_01', 'wc2.1_30s_wind_02', 'wc2.1_30s_wind_03', 
                         'wc2.1_30s_wind_04', 'wc2.1_30s_wind_05', 'wc2.1_30s_wind_06', 
                         'wc2.1_30s_wind_07', 'wc2.1_30s_wind_08', 'wc2.1_30s_wind_09', 
                         'wc2.1_30s_wind_10', 'wc2.1_30s_wind_11', 'wc2.1_30s_wind_12')


all_bio <- all_models[all_models$trait %in% monthly_prec_traits,]

# Adjust the trait levels for the order
all_bio$trait <- factor(all_bio$trait, levels = monthly_prec_traits)

# Create a ggplot layout similar to the single plot
p <- ggplot(all_bio, aes(x = as.numeric(factor(trait)), y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Label months as Jan, Feb, etc.
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(title = "Cross-Validation Results for wc2.1_30s_wind",
       x = "Month", y = "Prediction Accuracy (r.mean)")

p
# Save the plot as a PDF
ggsave("~/R/cannabis_GEAV/cross_validation/cross_validation_8X/wind/wind_combined8x.pdf", plot = p, width = 8, height = 3, units = "in")


#########################
#Supplemental Figure 19
##########################
#Predictions for (WorldClim 8X) for Ren (phylogenetic groupings)

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


######################################################
#Plots for 8X traits
######################################################
# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

####################
# Median Plot Line with Interquartile Range
####################

# Read the data
data <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_8X_2.csv')

# Reshape the data into a long format
data_long <- pivot_longer(data, 
                          cols = starts_with("wc2.1_30s"), 
                          names_to = c("variable", "month"), 
                          names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)", 
                          values_to = "value")

# Convert month to numeric so that ggplot orders the months correctly
data_long$month <- as.numeric(data_long$month)

# Ensure `Phylogenetic_tree` is a factor with correct levels before summarization
data_long$Phylogenetic_tree <- factor(data_long$Phylogenetic_tree, 
                                      levels = c("basal", "hemp-type", "drug-type", "drug-type feral"))

# Check for unexpected levels (for debugging)
print(unique(data_long$Phylogenetic_tree))

# Calculate the median, lower (25th) and upper (75th) quartiles for each group
summary_data <- data_long %>%
  group_by(variable, month, Phylogenetic_tree) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Remove the elevation variable
summary_data <- summary_data %>% filter(variable != "elevation")

# Specify the desired order of variables
variable_order <- c("tavg", "tmin", "tmax", "prec", "vapr", "wind", "srad")
facet_labels <- c(
  "tavg" = "Average temperature (tavg °C)",
  "tmin" = "Minimum temperature (tmin °C)",
  "tmax" = "Maximum temperature (tmax °C)",
  "prec" = "Precipitation (mm)",
  "vapr" = "Water vapor pressure (vapr kPa)",
  "wind" = "Wind speed (wind m/s)",
  "srad" = "Solar radiation (srad kJ/m²/day)"
)

# Reorder variables
summary_data$variable <- factor(summary_data$variable, levels = variable_order)

# Define color mapping explicitly
phylo_colors <- c("basal" = "orange", 
                  "hemp-type" = "green", 
                  "drug-type" = "red", 
                  "drug-type feral" = "blue")

# Define labels explicitly
phylo_labels <- c("Basal", "Hemp-type", "Drug-type", "Drug-type feral")

# Plot with lines for the median and ribbons for the quartiles
p <- ggplot(summary_data, aes(x = month, y = median_value, color = Phylogenetic_tree, fill = Phylogenetic_tree)) +
  geom_line(size = 1) +  # Median line
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), alpha = 0.2) +  # Shaded area for IQR
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = facet_labels)) +  # Custom facet labels
  labs(title = "", x = "Month", y = "GEAV") +  
  scale_color_manual(name = "Phylogenetic Tree", values = phylo_colors, labels = phylo_labels) +
  scale_fill_manual(name = "Phylogenetic Tree", values = phylo_colors, labels = phylo_labels) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 24),  # Adjust the size of the y-axis numbers
    axis.text.x = element_text(size = 24, angle = 0, vjust = 0.5, hjust = 0.5),  # Adjust x-axis numbers
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    strip.text = element_text(size = 24, face = "plain", hjust = 0),
    panel.spacing = unit(5, "lines")  # Increase space between plots
  ) +
  scale_x_continuous(breaks = seq(1, 12, 1))  # Show months as integers

# Print the plot
p

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S19_fixed.pdf", plot = p, width = 30, height = 18)


##################
#Elevation plot
##################
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the data
data <- read.csv('~/R/cannabis_GEAV/Outputs/rrblup_GEBV_values_cannabis_8X_2.csv')

# Create the boxplot
p <- ggplot(data, aes(x = Phylogenetic_tree, y = wc2.1_30s_elev, fill = Phylogenetic_tree)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +  # Add boxplot
  scale_fill_manual(
    name = "Phylogenetic Tree",
    values = c(
      "basal" = "orange",
      "hemp-type" = "green",
      "drug-type" = "red",
      "drug-type feral" = "blue"
    ),
    labels = c("Basal", "Hemp-type", "Drug-type", "Drug-type feral")  # Capitalized labels
  ) +
  labs(
    title = "Elevation (m)",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "plain", hjust = 0),  # Increase title size
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_text(size = 20),  # Increase y-axis text size
    axis.title.x = element_text(size = 20, face = "bold"),  # Style x-axis title
    axis.title.y = element_text(size = 20, face = "bold"),  # Style y-axis title
    legend.title = element_text(size = 10, face = "bold"),  # Style legend title
    legend.text = element_text(size = 20)  # Style legend text
  )

# Display the plot
print(p)

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S19_elev.pdf", plot = p, width = 12, height = 6)


#########################
#Supplemental Figure 20
##########################
#Cross validation (WorldClim 8X) for Ren + Soorni dataset
#HPC or Local O/N
library(rrBLUP)  # RR-BLUP package for Ridge Regression Best Linear Unbiased Prediction
library(hibayes) # hibayes package for Bayesian models
library(dplyr)

# Load custom functions for k-fold cross-validation
#source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/xval_kfold_functions.R")
#source("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_8x/xval_kfold_functions.R")
source("~/R/cannabis_GEAV/functions/xval_kfold_functions.R")


# Read genotype data in rrBLUP format
#data <- read.delim("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
data <- read.delim("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(data)[1] <- "rs." #needed to rename the sample column here
#colnames(data) <- gsub("^X", "", colnames(data)) #needed to remove the extra X added in 


# Read environmental data from CSV
#envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/ren_soorni_n111_extracted_climate_data.csv', head = TRUE)  # full environmental dataset
#envdat <- read.csv('/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_8x/ren_soorni_n111_extracted_climate_data_8X_ALL.csv', head = TRUE)
envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8X_ALL.csv', head = TRUE)

# Read training dataset from CSV (only core environmental data)
#trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_greedy/Ren_Soorni_n50_greedy_climate_data.csv", head = TRUE)
#trainingset <- read.csv("/home/ahmccorm/kantar_koastore/anna/Cannabis_GEAV/Ren_SOORNI/cross_validation_8x/ren_soorni_n111_extracted_climate_data_8X_core_n50.csv", head = TRUE)
trainingset <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8X_core_n50.csv", head = TRUE)

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
y.in <- envdat[, c(1, 4:ncol(envdat))]

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
saveRDS(xval_k10_rrblup, "~/R/cannabis_GEAV/Ren_Soorni/outputs/8X_xval_rrblup_kfold_10.RData")

#####
#PLOTS
#####
library(ggplot2)
setwd("~/R/cannabis_GEAV/Ren_Soorni/outputs/")

# RR-BLUP - K-Fold Cross-Validation
#Load cross-validation results for rrBLUP
rrblup_kfold10 <- readRDS("8X_xval_rrblup_kfold_10.RData")
rrblup_kfold10 <- rrblup_kfold10$xval.result
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)


# Organize all dataframes for merging
## Rename model names
rrblup_kfold10$model <- "rrBLUP"

## Input xval type
rrblup_kfold10$xval <- "Ten-Fold"

#Combine all model results into a single list
#model_list <- list(rrblup_kfold10, gauss_kfold_10, EXP_kfold_10,bayescpi_kfold_10)
model_list <- list(rrblup_kfold10)

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

###############
rrblup_results <- all_bio[all_bio$model == "rrBLUP", ]

# Export the results to a CSV file
write.csv(rrblup_results, "~/R/cannabis_GEAV/Ren_Soorni/cross_validation_19/8X_rrblup_prediction_accuracy_ren.csv", row.names = FALSE)
###############
library(ggplot2)
# Loop through each variable to create a combined plot for its 12 monthly values
for (variable_name in unique(all_models$variable)) {
  # Filter data for the current variable
  variable_data <- subset(all_models, variable == variable_name)
  
  # Create plot for the current variable across 12 months
  p <- ggplot(variable_data, aes(x = month, y = r.mean, color = model)) +
    theme_bw() +
    geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
    geom_point(size = 3) +
    ggtitle(paste("Cross-Validation Results for", variable_name)) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Labels months as Jan, Feb, etc.
    geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 15)) +
    ylim(0, 1) +
    labs(x = "Month", y = "Prediction Accuracy (r.mean)")
  
  # Save each plot to a file
  ggsave(filename = paste0("xval_", variable_name, "_monthly.pdf"), plot = p, width = 8, height = 3)
}

################
#Subset the data for elevation only
elevation_data <- subset(all_models, variable == "wc2.1_30s_elev.1")

# Plot for elevation with the same aesthetics
p_elev <- ggplot(elevation_data, aes(x = model, y = r.mean, color = model)) +
  theme_bw() +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), width = 0.3) +
  geom_point(size = 3) +
  ggtitle("Cross-Validation Results for Elevation") +
  geom_hline(yintercept = 0.5, color = "red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10),
        axis.text.y = element_text(size = 15)) +
  ylim(0, 1) +
  labs(x = "", y = "Prediction Accuracy (r.mean)")
p_elev
# Save the plot to a file
ggsave(filename = "xval_elevation_monthly.pdf", plot = p_elev, width = 8, height = 3)


#########################
#Supplemental Figure 21
##########################
#Predictions for (WorldClim 8X) for Ren + soorni datasets (PCA groupings)

################################
#GS on 8X - Ren + Soorni
################################
library(raster)
# Define the worldclim path
worldclim_path <- "~/R/Data/climate_data_8X/"

# Define the file patterns based on the types of climate variables
prec_files <- paste0(worldclim_path, "wc2.1_30s_prec_", sprintf("%02d", 1:12), ".tif")  # Precipitation
srad_files <- paste0(worldclim_path, "wc2.1_30s_srad_", sprintf("%02d", 1:12), ".tif")  # Solar radiation
tavg_files <- paste0(worldclim_path, "wc2.1_30s_tavg_", sprintf("%02d", 1:12), ".tif")  # Mean temperature
tmax_files <- paste0(worldclim_path, "wc2.1_30s_tmax_", sprintf("%02d", 1:12), ".tif")  # Maximum temperature
tmin_files <- paste0(worldclim_path, "wc2.1_30s_tmin_", sprintf("%02d", 1:12), ".tif")  # Minimum temperature
vapr_files <- paste0(worldclim_path, "wc2.1_30s_vapr_", sprintf("%02d", 1:12), ".tif")  # Water vapor pressure
wind_files <- paste0(worldclim_path, "wc2.1_30s_wind_", sprintf("%02d", 1:12), ".tif")  # Wind speed
elev_file <- paste0(worldclim_path, "wc2.1_30s_elev.tif")  # Elevation

# Load one of the precipitation raster files to use its CRS
sample_raster <- raster(prec_files[1])  # Load a raster to get the CRS

# Load the elevation layer
elev_layer <- raster(elev_file)

# Ensure all layers have the same CRS as the sample raster
crs(elev_layer) <- crs(sample_raster)  # Set CRS for elevation layer

# Now you can load the rest of the raster files and ensure they are consistent
prec_layers <- stack(prec_files)
srad_layers <- stack(srad_files)
tavg_layers <- stack(tavg_files)
tmax_layers <- stack(tmax_files)
tmin_layers <- stack(tmin_files)
vapr_layers <- stack(vapr_files)
wind_layers <- stack(wind_files)

# Elevation layer: replicate it 12 times to match other layers
elev_layers <- stack(replicate(12, elev_layer))

# Combine all layers into a single stack
climate_layers <- stack(prec_layers, srad_layers, tavg_layers, tmax_layers, tmin_layers, vapr_layers, wind_layers, elev_layers)

# Ensure the CRS for all layers is consistent
crs(climate_layers) <- crs(sample_raster)

# Continue with extracting the climate data for your coordinates
#data <- read.csv("~/R/cannabis_GEAV/DataSetMerged/inputs/ren_lw_lat_long.csv")
data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/inputs/ren_soorni_n111_for_climate.csv")
coords <- data.frame(lon = data$Longitude, lat = data$Latitude)

# Extract climate data for these coordinates
climate_values <- raster::extract(climate_layers, coords)

# Combine the extracted climate data with your original dataset
result <- cbind(data, climate_values)

# Save the result to a CSV file
write.csv(result, "~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8X_ALL.csv", row.names = FALSE)

######################################################
# Step 5: Genomic Selection with RR-BLUP
######################################################
library(rrBLUP)
library(dplyr)

# Load genotype and environmental data
gd1 <- read.delim("~/R/cannabis_GEAV/Ren_Soorni/outputs/Ren_Soorni_rrBLUP_format.txt", sep = "\t", header = TRUE)
colnames(gd1)[1] <- "rs."

envdat <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/ren_soorni_n111_extracted_climate_data_8X_core_n50.csv', head = TRUE)

# Set row names for genotype data (use 'rs.' column which contains SNP IDs)
row.names(gd1) <- gd1$rs.

# Remove the non-genotype columns ('rs.', 'allele', 'chrom', 'pos') for analysis
gd3 <- gd1[, -c(1)]  # SNP columns only
g.in <- as.matrix(gd3)  # Convert to matrix
g.in <- t(g.in)  # Transpose g.in to align with g.train structure (SNPs as columns, individuals as rows)
print(g.in [1:10, 1:10])  # Visualize 

#row.names(g.in) <- gsub("^X", "", row.names(g.in)) #accidental X added
#print(g.in [1:10, 1:10])  # Visualize 

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


# Subset g.pred to include only individuals in `pred`
pred <- setdiff(row.names(g.in), train)  # This gives # individuals not in the training set
g.pred <- g.in[pred, ]

##################################################
# Calculate the overall percentage of missing values in the entire matrix
total_missing <- sum(is.na(g.train))
total_values <- prod(dim(g.train))
percentage_missing <- (total_missing / total_values) * 100

cat("Overall percentage of missing data:", percentage_missing, "%\n")
#4.201%
##################################################

common_snps <- intersect(colnames(g.train), colnames(g.pred))

# Subset both g.train and g.pred to keep only the common SNPs
g.train <- g.train[, common_snps]
g.pred <- g.pred[, common_snps]

# Impute missing values in g.train and g.pred
# Impute missing values with the column mean
g.train[is.na(g.train)] <- apply(g.train, 2, function(x) mean(x, na.rm = TRUE)) #warning is ok
g.pred[is.na(g.pred)] <- apply(g.pred, 2, function(x) mean(x, na.rm = TRUE))#warning is ok

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
write.csv(gebv_df, '~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy_8X.csv')

########################
#Plotting Ren + Soorni dataset for 8X GEBVs 
########################
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the new data
data <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy_8X.csv')

# Reshape the data into a long format
data_long <- pivot_longer(data, 
                          cols = starts_with("wc2.1_30s"), 
                          names_to = c("variable", "month"), 
                          names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",
                          values_to = "value")

# Convert month to numeric so that ggplot orders the months correctly
data_long$month <- as.numeric(data_long$month)

# Ensure `Phylogenetic_tree` is a factor with correct levels
data_long$Phylogenetic_tree <- factor(data_long$Phylogenetic_tree, 
                                      levels = c("basal", "hemp-type", "drug-type", "drug-type feral", 
                                                 "Population_1", "Population_2"))

# Calculate the median, lower (25th), and upper (75th) quartiles for each group
summary_data <- data_long %>%
  group_by(variable, month, Phylogenetic_tree) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  )

# Remove the elevation variable if it exists
summary_data <- summary_data %>% filter(variable != "elevation")

# Specify the desired order of variables
variable_order <- c("tavg", "tmin", "tmax", "prec", "vapr", "wind", "srad")
facet_labels <- c(
  "tavg" = "Average temperature (tavg °C)",
  "tmin" = "Minimum temperature (tmin °C)",
  "tmax" = "Maximum temperature (tmax °C)",
  "prec" = "Precipitation (mm)",
  "vapr" = "Water vapor pressure (vapr kPa)",
  "wind" = "Wind speed (wind m/s)",
  "srad" = "Solar radiation (srad kJ/m²/day)"
)

# Reorder variables
summary_data$variable <- factor(summary_data$variable, levels = variable_order)

# Define color mapping explicitly
phylo_colors <- c("basal" = "orange", 
                  "hemp-type" = "green", 
                  "drug-type" = "red", 
                  "drug-type feral" = "blue",
                  "Population_1" = "#9370DB",
                  "Population_2" = "#008080")

# Define labels explicitly
phylo_labels <- c("Basal", "Hemp-type", "Drug-type", "Drug-type feral", "Population-1", "Population-2")

# Plot with lines for the median and ribbons for the quartiles
p <- ggplot(summary_data, aes(x = month, y = median_value, color = Phylogenetic_tree, fill = Phylogenetic_tree)) +
  geom_line(size = 1) +  # Median line
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), alpha = 0.2) +  # Shaded area for IQR
  facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = facet_labels)) +  # Custom facet labels
  labs(title = "Climate Variables Over Time by Phylogenetic Tree", 
       x = "Month", 
       y = "Value") + 
  scale_color_manual(name = "Phylogenetic Tree", values = phylo_colors, labels = phylo_labels) +
  scale_fill_manual(name = "Phylogenetic Tree", values = phylo_colors, labels = phylo_labels) +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 24),  # Adjust the size of the y-axis numbers
    axis.text.x = element_text(size = 24, angle = 0, vjust = 0.5, hjust = 0.5),  # Adjust x-axis numbers
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    strip.text = element_text(size = 24, face = "plain", hjust = 0),
    panel.spacing = unit(5, "lines")  # Increase space between plots
  ) + 
  scale_x_continuous(breaks = seq(1, 12, 1))  # Show months as integers

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S20_fixed.pdf", plot = p, width = 30, height = 18)

# Display the plot
p
################
#Elevation plot
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the new data
data <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/rrblup_GEBV_values_REN_SOORNI_n50_greedy_8X.csv')

# Create the boxplot
p <- ggplot(data, aes(x = Phylogenetic_tree, y = wc2.1_30s_elev.1, fill = Phylogenetic_tree)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +  # Add boxplot
  scale_fill_manual(
    name = "Phylogenetic Tree",
    values = c(
      "basal" = "orange",
      "hemp-type" = "green",
      "drug-type" = "red",
      "drug-type feral" = "blue",
      "Population_1" = "#9370DB",
      "Population_2" = "#008080"
    ),
    labels = c("Basal", "Hemp-type", "Drug-type", "Drug-type feral", "Population-1", "Population-2")
  ) +
  labs(
    title = "Elevation (m)",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "plain", hjust = 0),  # Left-align and increase title size
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_text(size = 28),  # Increase y-axis text size
    axis.title.x = element_text(size = 28, face = "bold"),  # Style x-axis title
    axis.title.y = element_text(size = 28, face = "bold"),  # Style y-axis title
    legend.title = element_text(size = 18, face = "bold"),  # Style legend title
    legend.text = element_text(size = 28)  # Style legend text
  )

# Display the plot
print(p)

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S21_elev.pdf", plot = p, width = 14, height = 8)


#########################
#Supplemental Figure 22
##########################
#plotting seasonal data by climate niche 
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the filtered data
filtered_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_8X_climate_niches_to_5above.csv")

# Reshape the data into a long format
seasonal_long <- filtered_data %>%
  pivot_longer(
    cols = starts_with("wc2.1_30s"),       # Select columns that start with "wc2.1_30s" for monthly climate variables
    names_to = c("variable", "month"),     # Split into variable (e.g., "prec") and month (e.g., "01")
    names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",  # Regex pattern to capture variable and month
    values_to = "value"
  ) %>%
  mutate(
    month = as.numeric(month)              # Convert month to numeric for ordering
  )

# Calculate median and quartiles by Climate_Class_Description
summary_data <- seasonal_long %>%
  group_by(variable, month, Climate_Class_Description) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  )

# Define custom colors for each Climate_Class_Description
custom_colors <- c(
  "BSk - Arid, steppe, cold" = "purple",
  "BSh - Arid, steppe, hot" = "purple3",
  "BWh - Arid, desert, hot" = "purple4",
  "BWk - Arid, desert, cold" = "orange",
  "Cwa - Temperate, dry winter, hot summer" = "blue",
  "Dfb - Cold, no dry season, warm summer" = "green",
  "Dsa - Cold, dry summer, hot summer" = "#1ABC9C"
)

# Facet labels for the variables
facet_labels <- c(
  "prec" = "Precipitation (mm)",
  "srad" = "Solar Radiation (kJ/m²/day)",
  "tavg" = "Average Temperature (°C)",
  "tmax" = "Max Temperature (°C)",
  "tmin" = "Min Temperature (°C)",
  "vapr" = "Water Vapor Pressure (kPa)",
  "wind" = "Wind Speed (m/s)"
)

# Reorder the variables for better presentation
# Reorder the variables for better presentation
variable_order <- c("tavg", "tmin", "tmax", "prec", "vapr", "wind", "srad")  # Desired order
facet_labels <- c(
  "tavg" = "Average Temperature (°C)",
  "tmin" = "Minimum Temperature (°C)",
  "tmax" = "Maximum Temperature (°C)",
  "prec" = "Precipitation (mm)",
  "vapr" = "Water Vapor Pressure (kPa)",
  "wind" = "Wind Speed (m/s)",
  "srad" = "Solar Radiation (kJ/m²/day)"
)

# Reorder variables in the summary_data
summary_data$variable <- factor(summary_data$variable, levels = variable_order)

# Filter out rows where the variable is NA
summary_data_filtered <- summary_data %>% 
  filter(!is.na(variable))

# Plot with updated variable order
p <- ggplot(summary_data_filtered, aes(x = month, y = median_value, color = Climate_Class_Description, fill = Climate_Class_Description)) +
  geom_line(size = 1) +  # Median line
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), alpha = 0.2) +  # Shaded area for IQR
  facet_wrap(~variable, scales = "free_y", labeller = labeller(variable = facet_labels)) +  # Custom facet labels
  labs(title = "",
       x = "Month",
       y = "") + 
  scale_color_manual(values = custom_colors) +  # Custom colors for lines
  scale_fill_manual(values = custom_colors) +   # Matching fill color for ribbon
  scale_x_continuous(breaks = seq(1, 12, 1), labels = seq(1, 12, 1)) +  # Set x-axis breaks and labels
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 24),  # Adjust the size of the y-axis numbers
    axis.text.x = element_text(size = 24, angle = 0, vjust = 0.5, hjust = 0.5),  # Adjust x-axis numbers
    axis.title.x = element_text(size = 24),  # Adjust the size of the x-axis title
    axis.title.y = element_text(size = 24),  # Adjust the size of the y-axis title
    plot.title = element_text(size = 24, face = "plain"),  # Adjust the plot title size
    legend.text = element_text(size = 24),  # Adjust the size of legend text
    legend.title = element_text(size = 24),  # Adjust the size of legend title
    strip.text = element_text(size = 24, face = "plain", hjust = 0),  # Adjust facet label size
    panel.spacing = unit(5, "lines")  # Increase space between plots
  )

# Display the plot
print(p)

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S21.pdf", plot = p, width = 30, height = 18)

####################
#elevation plot
####################
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the filtered data
filtered_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_8X_climate_niches_to_5above.csv")

# Define custom colors
custom_colors <- c(
  "BSk - Arid, steppe, cold" = "purple",
  "BSh - Arid, steppe, hot" = "purple3",
  "BWh - Arid, desert, hot" = "purple4",
  "BWk - Arid, desert, cold" = "orange",
  "Cwa - Temperate, dry winter, hot summer" = "blue",
  "Dfb - Cold, no dry season, warm summer" = "green",
  "Dsa - Cold, dry summer, hot summer" = "#1ABC9C"
)

# Create the boxplot

p <- ggplot(filtered_data, aes(x = Climate_Class_Description, y = wc2.1_30s_elev.1, fill = Climate_Class_Description)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +  # Add boxplot
  scale_fill_manual(
    name = "Climate Class",
    values = custom_colors
  ) +
  labs(
    title = "Elevation (m)",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 38, face = "plain", hjust = 0),  # Left-align and increase title size
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_text(size = 38),  # Increase y-axis text size
    axis.title.x = element_text(size = 38, face = "plain"),  # Style x-axis title
    axis.title.y = element_text(size = 38, face = "plain"),  # Style y-axis title
    legend.position = "none"  # Remove legend
  )



# Display the plot
print(p)

# Save the plot
ggsave("~/R/cannabis_GEAV/FIGS/S22_elev.pdf", plot = p, width = 14, height = 10)

##########################
#Supplemental Figure 23
##########################
############
#Figure S23 A
############
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the filtered data
filtered_data <- read.csv("~/R/cannabis_GEAV/Ren_Soorni/outputs/filtered_8X_climate_niches_to_5above.csv")

# Reshape the data into a long format
seasonal_long <- filtered_data %>%
  pivot_longer(
    cols = starts_with("wc2.1_30s"),       # Select columns that start with "wc2.1_30s" for monthly climate variables
    names_to = c("variable", "month"),     # Split into variable (e.g., "prec") and month (e.g., "01")
    names_pattern = "wc2.1_30s_([a-z]+)_(\\d+)",  # Regex pattern to capture variable and month
    values_to = "value"
  ) %>%
  mutate(
    month = as.numeric(month)              # Convert month to numeric for ordering
  )

# Filter for precipitation data only and Cwa climate class
cwa_prec_data <- seasonal_long %>%
  filter(variable == "prec", Climate_Class_Description == "Cwa - Temperate, dry winter, hot summer")

# Calculate median and quartiles for precipitation in Cwa climate class
summary_cwa_prec <- cwa_prec_data %>%
  group_by(month, Climate_Class_Description) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  )

# Plot seasonal precipitation data for Cwa climate class

summary_cwa_prec <- cwa_prec_data %>%
  group_by(month, Climate_Class_Description) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  )

# Plot with mean
p <- ggplot(summary_cwa_prec, aes(x = month, y = mean_value, fill = Climate_Class_Description)) +
  geom_line(color = "blue", size = 1) +  # Mean line in blue
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), fill = "blue", alpha = 0.2) +
  labs(
    title = "",
    x = "",
    y = "GEAV"
  ) + 
  scale_x_continuous(breaks = 1:12, labels = 1:12) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 25),  # Increase x-axis label size
    axis.text.y = element_text(size = 25),  # Increase y-axis label size
    axis.title.x = element_text(size = 25), # Increase x-axis title size
    axis.title.y = element_text(size = 25), # Increase y-axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25), # Increase legend title size
    plot.title = element_text(size = 25, face = "bold")  # Increase plot title size
  )

# Print the plot
p


ggsave("~/R/cannabis_GEAV/FIGS/S22A.pdf", plot = p, width = 15, height = 9)


############
#Figure 23 B
############
library(ggplot2)
library(dplyr)

# Load the data
individual_marker_effects_long <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_prec_per_SNP_june_aug.csv')

# Define selected individuals and chromosome order
selected_individuals <- c("SRR14708198", "SRR14708199", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258")
chromosome_labels <- c(
  "NC_044371.1" = "1 - NC_044371.1",
  "NC_044375.1" = "2 - NC_044375.1",
  "NC_044372.1" = "3 - NC_044372.1",
  "NC_044373.1" = "4 - NC_044373.1",
  "NC_044374.1" = "5 - NC_044374.1",
  "NC_044377.1" = "6 - NC_044377.1",
  "NC_044378.1" = "7 - NC_044378.1",
  "NC_044379.1" = "8 - NC_044379.1",
  "NC_044376.1" = "9 - NC_044376.1",
  "NC_044370.1" = "10 - NC_044370.1"
)

# Filter data for the selected individuals and format columns
filtered_data <- individual_marker_effects_long %>%
  filter(Individual %in% selected_individuals) %>%
  mutate(
    Chromosome = factor(sub("_[0-9]+$", "", SNP), levels = names(chromosome_labels), labels = chromosome_labels),
    Position = as.numeric(sub(".*_", "", SNP)) / 1e6,  # Convert position to Mb
    Month = factor(Month, levels = c("June", "August"))  # Ensure correct ordering of months
  )

# Plot marker effects with chromosome labels and x-axis in Mb
# Plot marker effects with chromosome labels and x-axis in Mb
p <- ggplot(filtered_data, aes(x = Position, y = Effect, color = Effect)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add dashed line at zero
  geom_hline(yintercept = 0.2, linetype = "solid", color = "red", size = 1) +  # Add red line at 0.2
  labs(title = "", x = "Genomic Position (Mb)", y = "Effect Size") +
  facet_wrap(~ Month + Chromosome, scales = "free_x", ncol = 10) +  # Separate plots for each month and chromosome
  scale_color_gradientn(colors = c("white", "blue"), 
                        name = "Effect Size", 
                        limits = c(-0.5, 0.5)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 24),  # Increase facet label size
    axis.text.x = element_text(size = 24, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
    axis.text.y = element_text(size = 24),  # Increase y-axis text size
    axis.title.x = element_text(size = 24),  # Increase x-axis title size
    axis.title.y = element_text(size = 24),  # Increase y-axis title size
    plot.title = element_text(size = 24),  # Increase plot title size
    legend.text = element_text(size = 24),  # Increase legend text size
    legend.title = element_text(size = 24)  # Increase legend title size
  ) +
  guides(color = guide_colorbar(title.position = "top"))

ggsave("~/R/cannabis_GEAV/FIGS/S22B.pdf", plot = p, width = 40, height = 20)


############
#Figure 23 C
############
# Boxplot showing significant difference in effect sizes between JUNE and august
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the data
individual_marker_effects_long <- read.csv('~/R/cannabis_GEAV/Ren_Soorni/outputs/individual_marker_effects_prec_per_SNP_june_aug.csv')
selected_individuals <- c("SRR14708198", "SRR14708199", "SRR14708228", "SRR14708248", "SRR14708251", "SRR14708258")

#selected_individuals <- c("SRR14708198", "SRR14708199", "SRR14708201", 
#"SRR14708228", "SRR14708245", "SRR14708246", 
#"SRR14708247", "SRR14708248", "SRR14708249", 
#"SRR14708250", "SRR14708251", "SRR14708252", 
#"SRR14708253", "SRR14708254", "SRR14708256", 
#"SRR14708257", "SRR14708258")


# Define chromosome order for consistent plotting
chromosome_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", "NC_044373.1", "NC_044374.1", 
                      "NC_044377.1", "NC_044378.1", "NC_044379.1", "NC_044376.1", "NC_044370.1")

# Step 1: Filter data for the selected individuals
data_long <- individual_marker_effects_long %>%
  filter(Individual %in% selected_individuals) %>%
  mutate(
    Chromosome = factor(sub("_[0-9]+$", "", SNP), levels = chromosome_order),
    Position = as.numeric(sub(".*_", "", SNP))
  )

# Step 2: Calculate the average effect size per SNP for June and August
data_avg <- data_long %>%
  group_by(SNP, Chromosome, Position, Month) %>%
  summarize(avg_effect = mean(Effect, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Month, values_from = avg_effect) %>%
  drop_na(June, August)  # Remove rows with NA in either June or August

# Step 3: Add a new column for visualization purposes
data_avg_long <- data_avg %>%
  pivot_longer(cols = c(June, August), names_to = "Month", values_to = "Effect") %>%
  mutate(Month = factor(Month, levels = c("June", "August")))  # Set order of levels

# Step 4: Perform a paired t-test across all SNPs on the average effect sizes between January and August
overall_t_test_result <- t.test(data_avg$June, data_avg$August, paired = TRUE)

# Print the overall t-test result
cat("Overall t-test on mean effect size differences between June and August:\n")
print(overall_t_test_result)

# Step 5: Create the boxplot
p<-ggplot(data_avg_long, aes(x = Month, y = Effect, fill = Month)) +
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black") +
  labs(
    title = "",
    subtitle = paste("Paired t-test p-value: 2.2e-16"),
    y = "Mean Effect Size",
    x = "Month"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 50, face = "bold"),       # Increase plot title size
    plot.subtitle = element_text(size = 50, face = "italic"), # Increase subtitle size
    axis.title.x = element_text(size = 50),                   # Increase x-axis title size
    axis.title.y = element_text(size = 50),                   # Increase y-axis title size
    axis.text.x = element_text(size = 50),                    # Increase x-axis label size
    axis.text.y = element_text(size = 50)                     # Increase y-axis label size
  ) +
  annotate("text", x = 1.5, y = max(data_avg_long$Effect), 
           label = ifelse(overall_t_test_result$p.value < 0.001, "***", "NS"),
           color = "red", size = 30)  # Increase annotation text size

p

ggsave("~/R/cannabis_GEAV/FIGS/S22C.pdf", plot = p, width = 30, height = 20)


############
#Figure 23 D
############
library(ggplot2)
library(dplyr)

file_path <- "~/R/cannabis_GEAV/Ren_Soorni/GEA/SNPs_with_Annotations_High_Effect_june_aug.csv"
data <- read.csv(file_path)
colnames(data)
unique_snp_data <- data %>% distinct(SNP_ID, Feature_Type, .keep_all = TRUE)

# Remove duplicates based on the combination of Feature_ID and Feature_Type
unique_genes_data <- data %>% distinct(Feature_ID, Feature_Type)

# Count the number of unique genes in each Feature_Type category
feature_counts_genes <- unique_genes_data %>%
  group_by(Feature_Type) %>%
  summarise(Count = n()) %>%
  ungroup()

# Display the resulting counts for verification
print(feature_counts_genes)

# Create the bar chart for unique genes
# Create the bar chart for unique genes
p<- ggplot(feature_counts_genes, aes(x = reorder(Feature_Type, -Count), y = Count, fill = Feature_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, color = "black", size = 18) +  # Increase label text size
  labs(
    title = "",  # Adjust or add a title here
    x = "Feature Type",  # Adjust x-axis title
    y = "Count"  # Adjust y-axis title
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Hide legend since labels are directly on bars
    axis.text.x = element_text(size = 40, angle = 45, hjust = 1),  # Increase x-axis text size
    axis.text.y = element_text(size = 40),  # Increase y-axis text size
    axis.title.x = element_text(size = 40),  # Increase x-axis title size
    axis.title.y = element_text(size = 40),  # Increase y-axis title size
    plot.title = element_text(size = 40, face = "bold")  # Increase title text size
  )
p
ggsave("~/R/cannabis_GEAV/FIGS/S22D.pdf", plot = p, width = 30, height = 20)


##########################
#Supplemental Figure 24
##########################
#GEA lat + temp

#All in LFMM_for_Latitude.R

##########################
#Supplemental Figure 25
##########################
#GEA ppt


#All in LFMM_for_Latitude.R

##########################
#Supplemental Figure 26
##########################
#Fst
##########################  WORKING
#step prior to R 
#vcftools --gzvcf merged_Ren_Soorni_snps.vcf.gz --weir-fst-pop Cwa_population_file.txt --weir-fst-pop BSk_population_file.txt --out Cwa_vs_BSk_fst
library(ggplot2)
library(dplyr)
library(tidyverse)

###########################
#basal versus drug-type feral
###########################
#fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/GEA/vcftools/Cwa_vs_BSk_fst.weir.fst", header = TRUE)
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Cwa_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus Cwa (drug-type feral)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")
p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_drugtype_feral.pdf", plot = p, device = "pdf", width = 15, height = 4)



###########################
#Basal versus pop2
###########################
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Dsa_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus Dsa (Population-2)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_pop2.pdf", plot = p, device = "pdf", width = 15, height = 4)



###########################
#Basal versus hemp
###########################
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Dfb_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus Dfb (Hemp)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_hemp.pdf", plot = p, device = "pdf", width = 15, height = 4)





###########################
#Basal versus Pop1 - Bsk
###########################
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Bsk_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus BSk (Population 1)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

p
#ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_pop1_bwk.pdf", plot = p, device = "pdf", width = 15, height = 4)
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_pop1_bsk.pdf", plot = p, device = "pdf", width = 15, height = 4)





###########################
#Basal versus Pop1 - Bsh
###########################
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Bsh_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus BSh (Population 1)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_pop1_bsh.pdf", plot = p, device = "pdf", width = 15, height = 4)




###########################
#Basal versus Pop1 - Bwh
###########################
fst_data <- read.table("~/R/cannabis_GEAV/Ren_Soorni/Fst/Bwk_vs_Bwh_fst.weir.fst", header = TRUE)    #basal versus dt-feral

# Define the custom chromosome order
chrom_order <- c("NC_044371.1", "NC_044375.1", "NC_044372.1", 
                 "NC_044373.1", "NC_044374.1", "NC_044377.1", 
                 "NC_044378.1", "NC_044379.1", "NC_044376.1", 
                 "NC_044370.1")
fst_data$CHROM <- factor(fst_data$CHROM, levels = chrom_order)

fst_data <- fst_data %>% 
  filter(WEIR_AND_COCKERHAM_FST >= 0)


# Define colors for each chromosome
chrom_colors <- c("NC_044371.1" = "#FF6347", "NC_044375.1" = "#4682B4", 
                  "NC_044372.1" = "#8A2BE2", "NC_044373.1" = "#3CB371", 
                  "NC_044374.1" = "#FFD700", "NC_044377.1" = "#DC143C", 
                  "NC_044378.1" = "#20B2AA", "NC_044379.1" = "#FF8C00", 
                  "NC_044376.1" = "#9370DB", "NC_044370.1" = "#32CD32")

# Calculate cumulative positions for each chromosome
chrom_offsets <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(offset = max(POS)) %>% 
  mutate(cumulative_offset = cumsum(lag(offset, default = 0)))

fst_data <- fst_data %>% 
  left_join(chrom_offsets, by = "CHROM") %>% 
  mutate(cumulative_POS = POS + cumulative_offset)

# Calculate midpoints for chromosome labels
chrom_labels <- fst_data %>% 
  group_by(CHROM) %>% 
  summarize(midpoint = mean(cumulative_POS)) %>% 
  mutate(label = case_when(
    CHROM == "NC_044371.1" ~ "1",
    CHROM == "NC_044375.1" ~ "2",
    CHROM == "NC_044372.1" ~ "3",
    CHROM == "NC_044373.1" ~ "4",
    CHROM == "NC_044374.1" ~ "5",
    CHROM == "NC_044377.1" ~ "6",
    CHROM == "NC_044378.1" ~ "7",
    CHROM == "NC_044379.1" ~ "8",
    CHROM == "NC_044376.1" ~ "9",
    CHROM == "NC_044370.1" ~ "10"
  ))

#####
# Calculate the 95th and 99th percentiles
fst_95th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
fst_99th <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)

# Add these lines to the plot
p<-ggplot(fst_data, aes(x = cumulative_POS, y = WEIR_AND_COCKERHAM_FST, color = CHROM)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_manual(values = chrom_colors, name = "Chromosome") +
  labs(title = "Fst across Genome Bwk (Basal) versus Bwh (Population 1)",
       x = "Cumulative Genomic Position",
       y = "Fst") +
  geom_hline(yintercept = fst_95th, linetype = "dashed", color = "blue", size = 0.7) +
  geom_hline(yintercept = fst_99th, linetype = "dashed", color = "red", size = 0.7) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_95th, label = "95th percentile", color = "blue", vjust = -0.5) +
  annotate("text", x = max(fst_data$cumulative_POS) * 0.9, y = fst_99th, label = "99th percentile", color = "red", vjust = -0.5) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add chromosome labels at midpoints
  geom_text(data = chrom_labels, aes(x = midpoint, y = 1.1, label = label), 
            inherit.aes = FALSE, size = 5, fontface = "bold")

p
ggsave("~/R/cannabis_GEAV/Ren_Soorni/FIGS/Fst_basal_vs_pop1_bwh.pdf", plot = p, device = "pdf", width = 15, height = 4)

###########################################################################################
#FIN
