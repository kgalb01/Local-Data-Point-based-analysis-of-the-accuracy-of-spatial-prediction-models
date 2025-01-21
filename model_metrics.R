# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2025-01-20
# Description: R Script to read and model metrics

# delete env
rm(list=ls())

# incorperate libraries
library(raster)
library(terra)
library(sf)
library(CAST)
library(ggplot2)
library(caret)
library(reshape2)
library(viridis)
library(tmap)
library(grid)
library(tidyr)

# set wd - this may have to be adjusted
setwd("path/to/data")

# load data
train_fiji <- st_read("path/to/fiji-lulc-2021-test-data-modified.geojson", quiet=TRUE) # the alternated data
train_rlp <- st_read("path/to/traindata_rlp.geojson", quiet=TRUE) #load train data

# repeat with modified data
train_fiji_modified <- st_read("path/to/fiji-lulc-2021-test-data-modified_notree.geojson", quiet=TRUE)
train_rlp_modified <- st_read("path/to/traindata_rlp_modified.geojson", quiet=TRUE)


fiji_10x10 <- rast("path/to/fiji_10x10.grd") 
rlp_10x10 <- rast("path/to/rlp_10x10_uncropped.grd") 
germ_60x60 <- rast("path/to/germ_60x60_uncropped.grd")

fiji_10x10_RF_model <- readRDS("path/to/fiji_10x10_RF_model.RDS")

rlp_10x10_RF_model <- readRDS("path/to/rlp_10x10_RF_model.RDS")

# repeat for modified data
fiji_10x10_RF_model_modified <- readRDS("path/to/fiji_10x10_RF_model_modified.RDS")

rlp_10x10_RF_model_modified <- readRDS("path/to/rlp_10x10_RF_model_modified.RDS")

################################################################################
################################################################################
################################################################################
# basic aoa, di and lpd calculation
fiji_self_10x10_aoa <-   CAST::aoa(fiji_10x10, 
                                   fiji_10x10_RF_model, 
                                   LPD = TRUE,
                                   maxLPD = 1,
                                   useWeight = TRUE,
                                   method = "L2",
                                   indices = FALSE)

rlp_fiji_10x10_aoa <-   CAST::aoa(fiji_10x10, 
                                  rlp_10x10_RF_model, 
                                  LPD = TRUE,
                                  maxLPD = 1,
                                  useWeight = TRUE,
                                  method = "L2",
                                  indices = FALSE)

# repeat for rlp
rlp_self_10x10_aoa <-   CAST::aoa(rlp_10x10, 
                                  rlp_10x10_RF_model, 
                                  LPD = TRUE,
                                  maxLPD = 1,
                                  useWeight = TRUE,
                                  method = "L2",
                                  indices = FALSE)

fiji_rlp_10x10_aoa <-   CAST::aoa(rlp_10x10, 
                                  fiji_10x10_RF_model, 
                                  LPD = TRUE,
                                  maxLPD = 1,
                                  useWeight = TRUE,
                                  method = "L2",
                                  indices = FALSE)

# repeat for germany
rlp_germ_60x60_aoa <-   CAST::aoa(germ_60x60, 
                                  rlp_10x10_RF_model, 
                                  LPD = TRUE,
                                  maxLPD = 1,
                                  useWeight = TRUE,
                                  method = "L2",
                                  indices = TRUE)

fiji_germ_60x60_aoa <-   CAST::aoa(germ_60x60, 
                                   fiji_10x10_RF_model, 
                                   LPD = TRUE,
                                   maxLPD = 1,
                                   useWeight = TRUE,
                                   method = "L2",
                                   indices = TRUE)

# saving the data
# Name: Model first, area second
writeRaster(fiji_self_10x10_aoa$AOA, "path/to/fiji_self_10x10_aoa.grd", overwrite = TRUE)
writeRaster(rlp_fiji_10x10_aoa$AOA, "path/to/rlp_fiji_10x10_aoa.grd", overwrite = TRUE)

# repeat for rlp
rlp_boundary <- st_read("path/to/rlp_bord.geojson", quiet = TRUE)
rlp_self_10x10_aoa_cropped <- crop(rlp_self_10x10_aoa$AOA, vect(rlp_boundary))
rlp_self_10x10_aoa_mask <- mask(rlp_self_10x10_aoa_cropped, vect(rlp_boundary))
rlp_self_10x10_aoa$AOA <- rlp_self_10x10_aoa_mask
writeRaster(rlp_self_10x10_aoa$AOA, "path/to/rlp_self_10x10_aoa.grd", overwrite = TRUE)

# repeat for RPL & germany
# ...
# repeat for DI...
# ...
# repeat for LPD ...
# ...

# repeat for modified data
# basic aoa, di and lpd calculation
fiji_self_10x10_modified_aoa <-   CAST::aoa(fiji_10x10, 
                                            fiji_10x10_RF_model_modified, 
                                            LPD = TRUE,
                                            maxLPD = 1,
                                            useWeight = TRUE,
                                            method = "L2",
                                            indices = FALSE)

# ...

# Load AOA and prediction data into R if necessary
# ...
fiji_self_10x10_aoa <- NULL
fiji_self_10x10_aoa$AOA <- rast("results/raster/fiji/aoa/fiji_self_10x10_aoa.grd")
fiji_self_10x10_aoa$DI <- rast("results/raster/fiji/aoa/fiji_self_10x10_di.grd")
fiji_self_10x10_aoa$LPD <- rast("results/raster/fiji/aoa/fiji_self_10x10_lpd.grd")
# ...
prediction_fiji_10x10 <- rast("results/raster/fiji/prediction/prediction_fiji_10x10.tif")
# ...


# Visualisation

#' calculate_rejection_rate:
#' This function calculates the rejection rate of pixels outside the Area of Applicability (AOA).
#' @param aoa A raster object containing the AOA values.
#'
#' @return A list containing the number of rejected pixels and the rejection rate as a percentage.
#' @export
#'
#' @examples
#' result <- calculate_rejection_rate(aoa)
#' rejection_rate <- result$rejection_rate
calculate_rejection_rate <- function(aoa) {
  total_pixels <- ncell(aoa$AOA)  # Total number of pixels in the AOA raster
  rejected_pixels <- sum(values(aoa$AOA) == 0, na.rm = TRUE)  # Count of pixels outside AOA
  rejection_rate <- (rejected_pixels / total_pixels) * 100  # Calculate rejection rate as a percentage
  return(list(rejected_pixels = rejected_pixels, rejection_rate = rejection_rate))  # Return results as a list
}


#' create_combined_map:
#' This function creates a combined map showing prediction results and areas rejected by the Area of Applicability (AOA).
#' @param aoa A raster object containing AOA values.
#' @param prediction A raster object with prediction results.
#' @param title The title of the map.
#' @param legend_colors A vector of colors for the prediction map legend.
#' @param rejection_rate The percentage of rejected pixels.
#' @param kappa The kappa statistic of the prediction model.
#'
#' @return A tmap object displaying the prediction map with rejected areas highlighted.
#' @export
#'
#' @examples
#' map <- create_combined_map(aoa, prediction, "My Map", legend_colors, rejection_rate, kappa)
#' tmap_save(map, "combined_map.png")
create_combined_map <- function(aoa, prediction, title, legend_colors, rejection_rate, kappa) {
  
  # Remove areas accepted by the AOA from the prediction map
  prediction[aoa$AOA == 0] <- NA 
  
  # Create a layer for rejected areas and fill them with magenta
  rejected_areas <- aoa$AOA
  rejected_areas[rejected_areas == 1] <- NA  # Set accepted areas to NA
  rejected_areas[rejected_areas == 0] <- 1   # Set rejected areas to 1
  
  # Create the map
  map <- tm_shape(prediction, raster.downsample = FALSE) +
    tm_raster("class", palette = legend_colors, title = "") + 
    tm_shape(rejected_areas) +  # Add the rejected areas
    tm_raster(style = "cat", palette = "black", legend.show = FALSE) +
    tm_scale_bar(bg.color = "white", bg.alpha = 0.5) +
    tm_layout(
      legend.position = c("right", "top"),
      legend.bg.color = "white",
      legend.bg.alpha = 0.6,
      outer.margins = c(0, 0, 0, 0), 
      inner.margins = c(0.05, 0.05, 0.15, 0.05),
      title = paste0("Rejected Area: ", round(rejection_rate, 2), "\nKappa: ", round(kappa, 2)),
      title.position = c("left", "top"),
      title.size = 0.4,
      title.bg.alpha = 0.5,
      main.title = title,
      main.title.size = 0.5,
      main.title.position = c("left", "top")
    ) +
    tm_legend(outside = FALSE) +
    tm_compass(position = c("left", "bottom"))
  
  return(map)
}



# Define color schemes for different prediction maps
prediction_colorscheme_fiji <- c("Agriculture" = "bisque1", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                 "Rock" = "azure4", "Shrubland" = "darkgoldenrod2", "Tree" = "darkgreen", 
                                 "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp <- c("Agriculture" = "bisque1","Urban" = "red","Vegetation" = "darkgreen","Water" = "cyan")

prediction_colorscheme_fiji_modified <- c("Agriculture" = "bisque1", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                          "Rock" = "azure4", "Shrubland" = "darkgoldenrod2", "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp_modified <- c("Agriculture" = "bisque1","Urban" = "red","Water" = "cyan")

# Create the maps for the four datasets
# Calculate the rejection rates for each dataset
rejection_rlp <- calculate_rejection_rate(fiji_self_10x10_aoa)
rejection_rlp_germ <- calculate_rejection_rate(rlp_germ_60x60_aoa)
rejection_rlp_fiji <- calculate_rejection_rate(rlp_fiji_10x10_aoa)


rlp_map <- create_combined_map(rlp_self_10x10_aoa, prediction_rlp_10x10, 
                               "AOA of Rhineland-Palatinate, Germany. Model of Rhineland-Palatinate (original data)", 
                               prediction_colorscheme_rlp, rejection_rate = rejection_rlp$rejection_rate, 
                               kappa = rlp_10x10_RF_model$results$Kappa)

rlp_germ_map <- create_combined_map(rlp_germ_60x60_aoa, prediction_rlp_germ_60x60, 
                                    "AOA of Germany. Model of Rhineland-Palatinate (original data)", 
                                    prediction_colorscheme_rlp, rejection_rate = rejection_rlp_germ$rejection_rate, 
                                    kappa = rlp_10x10_RF_model$results$Kappa)

rlp_fiji_map <- create_combined_map(rlp_fiji_10x10_aoa, prediction_rlp_fiji_10x10, 
                                    "AOA of Viti Levu, Fiji. Model of Rhineland-Palatinate (original data)", 
                                    prediction_colorscheme_rlp, 
                                    rejection_rate = rejection_rlp_fiji$rejection_rate, 
                                    kappa = rlp_10x10_RF_model$results$Kappa)


# Combine the four maps into a single overview map
combined_map_rlp_model <- tmap_arrange(rlp_map, rlp_germ_map, rlp_fiji_map, ncol = 3)

# Save the combined map to a file
tmap_save(combined_map_rlp_model, "comparison_aoa_rlp_model.png", width = 20, height = 15, units = "cm")

# repeat for the other models
# ...


#' create_di_map_from_raster:
#' This function creates a map visualizing the Dissimilarity Index (DI) from a raster object.
#' @param di_raster A raster object containing DI values.
#' @param main_title The title of the map.
#' @param legend_title The title for the legend.
#' @param min_value The minimum DI value to cap the raster at.
#' @param max_value The maximum DI value to cap the raster at.
#'
#' @return A tmap object displaying the DI map with a continuous color scale.
#' @export
#'
#' @examples
#' di_map <- create_di_map_from_raster(di_raster, "DI Map", "Dissimilarity Index", 0, 1)
#' tmap_save(di_map, "di_map.png")
create_di_map_from_raster <- function(di_raster, main_title, legend_title, min_value, max_value) {
  
  # Cap DI values at the global min and max values
  di_capped <- clamp(di_raster, lower = min_value, upper = max_value)
  
  # Define the map shape and add a raster layer with a continuous color scale
  map <- tm_shape(di_capped) +
    tm_raster(palette = "viridis", title = legend_title, style = "cont", 
              breaks = c(min_value, max_value), legend.is.portrait = TRUE) +
    tm_layout(
      legend.position = c("right", "top"),
      legend.bg.color = "white",
      legend.bg.alpha = 0.2,
      outer.margins = c(0, 0, 0, 0),
      inner.margins = c(0.05, 0.05, 0.15, 0.05),
      main.title = main_title,                  
      main.title.size = 0.345,
      main.title.position = c("center", "top"), 
      title.size = 0.3                          
    ) +
    tm_compass(position = c("left", "bottom")) +
    tm_scale_bar(bg.color = "white", bg.alpha = 0.1)
  
  return(map)
}


# Combine all DI rasters to compute the global min and max values
all_di_rasters <- list(
  rlp_self_10x10_aoa$DI,
  rlp_germ_60x60_aoa$DI,
  rlp_fiji_10x10_aoa$DI,
  fiji_self_10x10_aoa$DI,
  fiji_rlp_10x10_aoa$DI,
  rlp_self_10x10_modified_aoa$DI,
  rlp_germ_60x60_modified_aoa$DI,
  fiji_rlp_10x10_modified_aoa$DI,
  fiji_self_10x10_modified_aoa$DI,
  rlp_fiji_10x10_modified_aoa$DI
)

# Calculate the global min and max DI values across all datasets
min_di_value <- min(sapply(all_di_rasters, function(r) min(values(r), na.rm = TRUE)), na.rm = TRUE)
max_di_value <- max(sapply(all_di_rasters, function(r) quantile(values(r), 0.95, na.rm = TRUE)), na.rm = TRUE)

# Generate DI maps using the global scale
# Generate DI maps using the global scale
rlp_di_map <- create_di_map_from_raster(rlp_self_10x10_aoa$DI, 
                                        "95th Percentile DI of Rhineland-Palatinate, Germany. Model of Rhineland-Palatinate (original data)", 
                                        "DI", 
                                        min_value = min_di_value, max_value = max_di_value)

rlp_germ_di_map <- create_di_map_from_raster(rlp_germ_60x60_aoa$DI, 
                                             "95th Percentile DI of Germany. Model of Rhineland-Palatinate (Original Data)", 
                                             "DI", 
                                             min_value = min_di_value, max_value = max_di_value)

rlp_fiji_di_map <- create_di_map_from_raster(rlp_fiji_10x10_aoa$DI, 
                                             "95th Percentile DI of Viti Levu, Fiji - Rhineland-Palatinate Model (original data)", 
                                             "DI", 
                                             min_value = min_di_value, max_value = max_di_value)

# Arrange the four DI maps into a combined layout and save to a file
combined_di_map_rlp_model <- tmap_arrange(rlp_di_map, rlp_germ_di_map, rlp_fiji_di_map, ncol = 3)


tmap_save(combined_di_map_rlp_model, "comparison_di_rlp_model.png", width = 20, height = 15, units = "cm")

# repeat for RLP and GERM
# ...



#' create_lpd_map_from_raster:
#' This function creates a map visualizing the Local Data Point Density (LPD) from a raster object.
#' @param lpd_raster A raster object containing LPD values.
#' @param main_title The title of the map.
#' @param legend_title The title for the legend.
#'
#' @return A tmap object displaying the LPD map with a continuous color scale.
#' @export
#'
#' @examples
#' lpd_map <- create_lpd_map_from_raster(lpd_raster, "LPD Map", "Local Data Point Density")
#' tmap_save(lpd_map, "lpd_map.png")
create_lpd_map_from_raster <- function(lpd_raster, main_title, legend_title) {
  # Define the map shape using the LPD raster data
  map <- tm_shape(lpd_raster) +
    
    # Add the raster layer with a continuous color scale (viridis palette)
    tm_raster(palette = "viridis", title = legend_title, style = "cont", 
              breaks = c(0, 50, 100, 150, 200), legend.is.portrait = TRUE) +
    
    # Customize the layout of the map
    tm_layout(
      legend.position = c("right", "top"),  
      legend.bg.color = "white",  
      legend.bg.alpha = 0.2,  
      outer.margins = c(0, 0, 0, 0),  
      inner.margins = c(0.05, 0.05, 0.15, 0.05),  
      main.title = main_title,
      main.title.size = 0.345,
      main.title.position = c("center", "top")
    ) +
    
    # Add a compass to the map
    tm_compass(position = c("left", "bottom")) +
    
    # Add a scale bar to the map
    tm_scale_bar(bg.color = "white", bg.alpha = 0.1)
  
  return(map)
}


rlp_lpd_map <- create_lpd_map_from_raster(rlp_self_10x10_aoa$LPD, 
                                          "LPD of Rhineland-Palatinate, Germany. Model of Rhineland-Palatinate (original data)", 
                                          "LPD")

rlp_germ_lpd_map <- create_lpd_map_from_raster(rlp_germ_60x60_aoa$LPD, 
                                               "LPD of Germany. Model of Rhineland-Palatinate (original data)", 
                                               "LPD")

rlp_fiji_lpd_map <- create_lpd_map_from_raster(rlp_fiji_10x10_aoa$LPD, 
                                               "LPD of Viti Levu, Fiji. Model of Rhineland-Palatinate (original data)", 
                                               "LPD")


# Arrange the four LPD maps into a combined layout
combined_lpd_map_rlp_model <- tmap_arrange(rlp_lpd_map, rlp_germ_lpd_map, rlp_fiji_lpd_map, ncol = 3)


# Save the combined map to a file
tmap_save(combined_lpd_map_rlp_model, "comparison_lpd_rlp_model.png", width = 20, height = 15, units = "cm")

# repeat for RLP and GERM
# ...
