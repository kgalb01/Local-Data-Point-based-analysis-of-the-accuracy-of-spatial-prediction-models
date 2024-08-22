# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2024-06-08
# Description: R Script to slightly modify trainingdata from fiji
# for easier processing during modeltraining
# Link to the data: https://pacificdata.org/data/dataset/fiji-land-use-land-cover-test-dataset

# delete env
rm(list=ls())

# incorperate lib
library(dplyr)
library(jsonlite)
library(sf)

# set wd - this may have to be adjusted
setwd("path/to/data")

# input data
input_file <- "fiji-lulc-2021-test-data.geojson"
fiji_data <- fromJSON(input_file)

# map the ref_class values to the corresponding name of the class
class_map <- c("1" = "Water", "2" = "Mangrove", "3" = "Rock", 
               "4" = "Urban", "5" = "Agriculture", 
               "6" = "Grassland", "7" = "Shrubland", "8" = "Tree")

# modify fiji data
# Note: Might get in trouble because of the "[]" around each "Label" and "ClassID"
modified <- lapply(seq_along(fiji_data$features$type), function(i) {
  feature <- list(
    type = fiji_data$features$type[i],
    properties = list(
      Label = class_map[as.character(fiji_data$features$properties$ref_class[i])], # replace "ref_class" with corresponding class name
      ClassID = fiji_data$features$properties$ref_class[i], # include "ClassID" incase it is needed
      fid = i # attach "fid" to each point incase it is needed
    ),
    geometry = list(
      type = fiji_data$features$geometry$type[i], # include information about type of geometry
      coordinates = fiji_data$features$geometry$coordinates[[i]] # include information about coords
    )
  )
  return(feature) # return modified feature list
})

# create feature collection out of existing points to be valid geojson
modified_fiji_data <- list(type = "FeatureCollection", features = modified)

# save the modified data
# Note: can't use "st_write" because at this point "modified_fiji_data" is a list of features
# converting that to a dataframe would be possible but this way it's easier
output_file <- "modified_fiji_training_data.geojson"
write_json(modified_fiji_data, path = output_file)

