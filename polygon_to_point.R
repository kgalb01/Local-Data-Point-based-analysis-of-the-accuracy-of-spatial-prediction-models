# Local Data Point Density-basierte Analyse der Genauigkeit von räumlichen Vorhersagemodellen
# Author: Kieran Galbraith
# Date: 2024-06-08
# Description: R Script to convert DLM to trainingsdata for Rhineland-Palatinate
# Source for DLM: https://lvermgeo.rlp.de/produkte/geotopografie/digitale-landschaftsmodelle-dlm/digitales-basislandschaftsmodell-basis-dlm

# delete env
rm(list=ls())

# incorporate lib
library(sf)
library(dplyr)
library(sp)
library(lwgeom)

# set global variables
seed = 123 # seed for reproducibility

#' polygon_to_point:
#' This function extracts a specific number of random points out of polygons
#' from shapefiles
#' @param data Polygons, n number of random points
#'
#' @return A specified number of random selected points of a shapefile
#' @export
#'
#' @examples
#' point <- polygon_to_point(polygon, n = 100, seed = 123)
polygon_to_point <- function(polygon, n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  polygon <- st_set_crs(polygon, 32632) # set the current CRS
  polygon <- st_transform(polygon, 4326) # transform CRS to EPSG:4326 for later use
  points <- st_sample(polygon, n)
  return(points)
}

# set wd - this may have to be adjusted
setwd("path/to/data")

# Class: Water
water_rlp <- st_read("path/to/DLM_RP_gew01f.shp", quiet = TRUE) # water shapes rhineland palatinate
# select 100 random points
water_combined_points_rlp <- polygon_to_point(water_rlp, n = 100, seed)
# labeling the classes properly to get valid geojson
water_rlp <- st_sf(geometry = water_combined_points_rlp)
water_rlp <- water_rlp %>%
  mutate(Label = "Water",
         ClassID = 1,
         fid = 1:nrow(water_rlp))

# save the data
water_combined_points_output <- "water_rlp.geojson"
st_write(water_rlp, water_combined_points_output, driver = "GeoJSON")

# Class: Agriculture
agriculture_rlp <- st_read("path/to/DLM_RP_veg01f.shp", quiet = TRUE) # rhineland palatinate
# select 259 random points
agriculture_rlp_points <- polygon_to_point(agriculture_rlp, n = 308, seed)
# labeling the classes properly to get valid geojson
agriculture_rlp <- st_sf(geometry = agriculture_rlp_points)
agriculture_rlp <- agriculture_rlp %>%
  mutate(Label = "Agriculture",
         ClassID = 2,
         fid = nrow(water_rlp)+1:nrow(agriculture_rlp))

# save the data
agriculture_combined_points_output <- "agriculture_rlp.geojson"
st_write(agriculture_rlp, agriculture_combined_points_output, driver = "GeoJSON")

# Class: Urban (combining Urban and Traffic)
urban_rlp <- st_read("path/to/.shp", quiet = TRUE) # aparment, mining, graveyards, ... rhineland
urban_2_rlp <- st_read("path/to/.shp", quiet = TRUE) # sports facilities, historical buildings, ... rhineland
urban_3_rlp <- st_read("path/to/.shp", quiet = TRUE) # habours rhineland
urban_4_rlp <- st_read("path/to/.shp", quiet = TRUE) # towers rhineland
urban_rlp <- bind_rows(urban_rlp, urban_2_rlp, urban_3_rlp, urban_4_rlp) # urban shapes rhineland

trafic_rlp <- st_read("path/to/.shp", quiet = TRUE) # streets, sealed area rhineland palatinate
trafic_2_rlp <- st_read("path/to/.shp", quiet = TRUE) # train tracks rhineland palatinate
trafic_rlp <- bind_rows(trafic_rlp, trafic_2_rlp) # trafic shapes rhineland palatinate

urban_rlp_combined <- bind_rows(urban_rlp, trafic_rlp)

# select 130 random points
urban_combined_points_rlp <- polygon_to_point(urban_rlp_combined, n = 118, seed)
# labeling the classes properly to get valid geojson
urban_rlp <- st_sf(geometry = urban_combined_points_rlp)
urban_rlp <- urban_rlp %>%
  mutate(Label = "Urban",
         ClassID = 3,
         fid = nrow(agriculture_rlp)+1:nrow(urban_rlp))

# save the data
urban_combined_points_output <- "urban_rlp.geojson"
st_write(urban_rlp, urban_combined_points_output, driver = "GeoJSON")

# Class: Vegetation
vegetation_1_rlp <- st_read("path/to/DLM_RP_veg02f.shp", quiet = TRUE) # forest rhineland palatinate
vegetation_2_rlp <- st_read("path/to/DLM_RP_veg03f.shp", quiet = TRUE) # wood, bog, swamp, ... rhineland palatinate
vegetation_rlp <- bind_rows(vegetation_1_rlp, vegetation_2_rlp) # vegetation shapes rhineland palatinate
# select 311 random points
vegetation_combined_points_rlp <- polygon_to_point(vegetation_rlp, n = 308, seed)
# labeling the classes properly to get valid geojson
vegetation_rlp <- st_sf(geometry = vegetation_combined_points_rlp)
vegetation_rlp <- vegetation_rlp %>%
  mutate(Label = "Vegetation",
         ClassID = 4,
         fid = nrow(urban_rlp)+1:nrow(vegetation_rlp))

# save the data
vegetation_combined_points_output <- "vegetation_rlp.geojson"
st_write(vegetation_rlp, vegetation_combined_points_output, driver = "GeoJSON")

# combine all extracted points into one - this is the training data
water_combined_points_rlp <- st_read("water_rlp.geojson", quiet = TRUE)
agriculture_combined_points_rlp <- st_read("agriculture_rlp.geojson", quiet = TRUE)
urban_combined_points_rlp <- st_read("urban_rlp.geojson", quiet = TRUE)
vegetation_combined_points_rlp <- st_read("vegetation_rlp.geojson", quiet = TRUE)

combined_rlp <- bind_rows(water_combined_points_rlp, agriculture_combined_points_rlp, urban_combined_points_rlp, vegetation_combined_points_rlp)

# save the data
combined_output <- "traindata_rlp.geojson"
st_write(combined_rlp, combined_output, driver = "GeoJSON")
