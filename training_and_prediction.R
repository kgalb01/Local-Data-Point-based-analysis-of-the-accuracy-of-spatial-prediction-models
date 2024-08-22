# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2024-06-09
# Description: R Script to predict data on basis of random forest trained model

# delete env
rm(list=ls())

# incorperate libraries
library(raster)
library(terra)
library(caret)
library(sf)
library(CAST)
library(tmap)
library(cluster)

# set wd - this may have to be adjusted
setwd("path/to/data")

# load data
train_fiji <- st_read("path/to/fiji-lulc-2021-test-data-modified.geojson", quiet=TRUE) # the alternated data
train_rlp <- st_read("path/to/traindata_rlp.geojson", quiet=TRUE) #load train data

# repeat with modified data
train_fiji_modified <- st_read("path/to/fiji-lulc-2021-test-data-modified_notree.geojson", quiet=TRUE)
train_rlp_modified <- st_read("path/to/traindata_rlp_modified.geojson", quiet=TRUE)

train_fiji_extr <- readRDS("path/to/train_fiji_extr.rds")
train_fiji_extr <- train_fiji_extr[, !colnames(train_fiji_extr) %in% ".geo"]
train_fiji_extr <- train_fiji_extr[, !colnames(train_fiji_extr) %in% "system:index"]

# repeat for RLP 
# ...

# repeat with modified data
# ...



fiji_10x10 <- rast("path/to/fiji_10x10.grd")

rlp_10x10 <- rast("path/to/rlp_10x10_uncropped.grd")

germ_60x60 <- rast("path/to/germ_60x60_uncropped.grd")




#' train_model:
#' This function trains a Random Forest model using spatial cross-validation based on spatial clusters.
#' @param extr_data A data frame containing predictor variables and labels for training.
#' @param training_data A spatial object with coordinates used for clustering and cross-validation.
#' @param ntree The number of trees to grow in the Random Forest model (default is 500).
#' @param k The number of cross-validation folds (default is 10).
#' @param seed A seed for reproducibility (default is 321).
#' @param num_clusters The number of spatial clusters to create for cross-validation (default is 5).
#'
#' @return A trained Random Forest model object.
#' @export
#'
#' @examples
#' model <- train_model(extr_data, training_data, ntree = 500, k = 10, seed = 321, num_clusters = 5)
train_model <- function(extr_data, training_data, ntree = 500, k = 10, seed = 321, num_clusters = 5) {
  
  # Set a seed for reproducibility
  set.seed(seed)
  
  # Extract spatial coordinates from the training data
  coords <- st_coordinates(training_data)
  
  # Perform k-means clustering on the spatial coordinates to create spatial clusters
  clusters <- kmeans(coords, centers = num_clusters, nstart = 25)
  
  # Assign the spatial clusters as a new factor variable in the training data
  training_data$spatial_cluster <- as.factor(clusters$cluster)
  
  # Create spatial cross-validation folds based on the spatial clusters
  folds <- CreateSpacetimeFolds(training_data, spacevar = "spatial_cluster", k = k)
  
  # Set up training control parameters for cross-validation using the spatial folds
  train_control <- trainControl(method = "cv", index = folds$index, savePredictions = "final")
  
  # Define predictor variables by excluding the label and geometry columns
  predictors <- setdiff(names(extr_data), c("Label", "geometry"))
  
  # Train a Random Forest model using the extracted data and specified parameters
  model <- train(extr_data[, predictors], 
                 extr_data$Label,
                 method = "rf",
                 trControl = train_control,
                 tuneLength = 10,       
                 ntree = ntree,         
                 importance = TRUE,     
                 metric = "Kappa",      
                 seed = seed)           
  
  # Return the trained model
  return(model)
}


# Beispielaufruf der Funktion
fiji_10x10_RF_model <- train_model(train_fiji_extr, train_fiji, num_clusters = 10)
rlp_10x10_RF_model <- train_model(train_rlp_extr, train_rlp, num_clusters = 10)


# repeat with modified data
# ...

# save the created models
saveRDS(fiji_10x10_RF_model,file="path/to/fiji_10x10_RF_model_spatCV.RDS") 
# ...

# repeat with modified data
# ...

# predict
# NOTE: In file name model location is named first, then prediction location
# NOTE: The resolution named is also the resolution of the used model
fiji_10x10[is.na(fiji_10x10)] <- 0
prediction_fiji_10x10 <- predict(fiji_10x10, fiji_10x10_RF_model) # Model: Fiji; Prediction: Fiji
prediction_rlp_fiji_10x10 <- predict(fiji_10x10, rlp_10x10_RF_model) # Model: Rhineland; Prediction: Fiji

# repeat for RLP and GERM
# ...


# save the prediction data
writeRaster(prediction_fiji_10x10, "path/to/prediction_fiji_10x10_spatCV.tif", overwrite = TRUE)
writeRaster(prediction_rlp_fiji_10x10, "path/to/prediction_rlp_fiji_10x10_spatCV.tif", overwrite = TRUE)

# repeat the other way around
rlp_boundary <- st_read("path/to/rlp_bord.geojson", quiet = TRUE)
prediction_rlp_10x10_cropped <- crop(prediction_rlp_10x10, vect(rlp_boundary))
prediction_rlp_10x10_mask <- mask(prediction_rlp_10x10_cropped, vect(rlp_boundary))
prediction_rlp_10x10 <- prediction_rlp_10x10_mask
writeRaster(prediction_rlp_10x10, "path/to/prediction_rlp_10x10_spatCV.tif", overwrite = TRUE)

# repeat for RLP and GERM
# ...

#################
# set prediction color scheme for better plotting
rediction_colorscheme_fiji <- c("Agriculture" = "yellow", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                 "Rock" = "azure4", "Shrubland" = "brown", "Tree" = "darkgreen", 
                                 "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp <- c("Agriculture" = "yellow","Urban" = "red","Vegetation" = "darkgreen","Water" = "cyan")

# Ensure the colors match the classes
labels <- levels(factor(fiji_10x10_RF_model$levels))

if (length(labels) > 0) {
  # Match the colors to the labels
  legend_colors <- prediction_colorscheme_fiji[labels]
  
  # Convert the raster values to factors and add the labels
  prediction_fiji_10x10_raster <- raster(prediction_fiji_10x10)
  prediction_fiji_10x10_raster <- ratify(prediction_fiji_10x10_raster)
  rat <- data.frame(ID = 1:length(labels), class = labels)
  levels(prediction_fiji_10x10_raster) <- list(rat)
  
  # Create a map with tmap
  map <- tm_shape(prediction_fiji_10x10_raster, raster.downsample = FALSE) +
    tm_raster("class", palette = legend_colors, title = "") + 
    tm_scale_bar(bg.color = "white", bg.alpha = 0.5) +
    tm_layout(
      legend.position = c("right", "top"),
      legend.bg.color = "white",
      legend.bg.alpha = 0.8,
      outer.margins = c(0.1, 0.1, 0.1, 0.1), 
      title = "Land Use Classification of Viti Levu, Fiji, using a Spatially Validated Random Forest Model",
      title.position = c("left", "top"),
      title.size = 0.4
    ) +
    tm_compass(position = c("left", "bottom")) 
  
  # Save the map to a file
  tmap_save(map, "prediction_fiji_10x10.png")
} else {
  print("No labels found. Check the model output.")
}

# Save the map to a file
tmap_save(map, "prediction_fiji_10x10_spatCV.png")

# repeat for RLP and GERM
# ...

##################### repeat with modified data ######################
# ...

# set prediction color scheme for better plotting
prediction_colorscheme_fiji_modified <- c("Agriculture" = "yellow", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                          "Rock" = "black", "Shrubland" = "brown", "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp_modified <- c("Agriculture" = "yellow","Urban" = "red","Water" = "cyan")

# ...