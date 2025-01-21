# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2025-01-20
# Description: R Script to process Sentinel-2 data from Fijis and Germany
# data downloaded from Google Earth Engine: https://code.earthengine.google.com/

# delete env
rm(list=ls())

# incorperate libraries
library(raster)
library(terra)
library(sf)

# set working directory
setwd("path/to/data")
list.files()

#' calculate_indices:
#' This function calculates several spectral indices from satellite imagery bands.
#' @param image A list or data frame containing the bands of the satellite image (e.g., B2, B3, B4, B8, B11).
#'
#' @return A named vector of calculated indices: NDVI, NDWI, MNDWI, NDBI, GCVI, and EVI.
#' @export
#'
#' @examples
#' indices <- calculate_indices(image)
#' plot(indices["NDVI"])
calculate_indices <- function(image) {
  ndvi <- (image$B8 - image$B4) / (image$B8 + image$B4)
  ndwi <- (image$B3 - image$B8) / (image$B3 + image$B8)
  mndwi <- (image$B3 - image$B11) / (image$B3 + image$B11)
  ndbi <- (image$B11 - image$B8) / (image$B11 + image$B8)
  gcvi <- image$B8 / image$B3 - 1
  evi <- 2.5 * ((image$B8 - image$B4) / (image$B8 + 6 * image$B4 - 7.5 * image$B2 + 1))
  
  indices <- c(ndvi, ndwi, mndwi, ndbi, gcvi, evi)
  names(indices) <- c("NDVI", "NDWI", "MNDWI", "NDBI", "GCVI", "EVI")
  
  return(indices)
}

# load data for Fiji
fiji_10x10_1 <- rast("path/to/Sen2Fiji_10x10_cropped_1.tif")
fiji_10x10_1 <- aggregate(fiji_10x10_1, 3) # aggregate the higher resolution for faster processing
fiji_10x10_2 <- rast("path/to/Sen2Fiji_10x10_cropped_2.tif")
fiji_10x10_2 <- aggregate(fiji_10x10_2, 3) # aggregate the higher resolution for faster processing
fiji_10x10_3 <- rast("path/to/Sen2Fiji_10x10_cropped_3.tif")
fiji_10x10_3 <- aggregate(fiji_10x10_3, 3) # aggregate the higher resolution for faster processing
fiji_10x10_4 <- rast("path/to/Sen2Fiji_10x10_cropped_4.tif")
fiji_10x10_4 <- aggregate(fiji_10x10_4, 3) # aggregate the higher resolution for faster processing
fiji_10x10 <- terra::mosaic(fiji_10x10_1, fiji_10x10_2, fiji_10x10_3, fiji_10x10_4, overwrite=TRUE)
fiji_10x10 <- aggregate(fiji_10x10, 3) # aggregate the higher resolution for faster processing

fiji_ndvi <- ((fiji_10x10$B8 - fiji_10x10$B4)/(fiji_10x10$B8 + fiji_10x10$B4))

# the indices of the fiji data set seem to be broken
fiji_10x10 <- fiji_10x10[[1:12]]  # Assuming the first 12 bands are the original bands
fiji_indices <- calculate_indices(fiji_10x10)
fiji_10x10 <- c(fiji_10x10, fiji_indices)
# save the progress so far
writeRaster(fiji_10x10, "path/to/fiji_10x10.grd", overwrite = TRUE)

# load data for RLP and GERM
# ...
