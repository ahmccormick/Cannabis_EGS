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


