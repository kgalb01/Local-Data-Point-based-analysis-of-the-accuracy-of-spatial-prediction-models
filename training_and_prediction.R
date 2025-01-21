# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2025-20-01
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
library(blockCV)
library(mapac)

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
#' This function trains a Random Forest model using spatial cross-validation with kNNDM.
#' @param extr_data A data frame containing predictor variables and labels for training.
#' @param training_data A spatial object with coordinates used for creating spatial folds.
#' @param raster_data A raster object representing the study area for spatial fold creation.
#' @param ntree The number of trees to grow in the Random Forest model (default is 500).
#' @param k The number of cross-validation folds to create using kNNDM (default is 10).
#' @param seed A seed for reproducibility (default is 321).
#'
#' @return A trained Random Forest model object.
#' @export
#'
#' @examples
#' model <- train_model(extr_data, training_data, raster_data, ntree = 500, k = 10, seed = 321)
train_model <- function(extr_data, training_data, raster_data, ntree = 500, k = 10, seed = 321) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Create spatial folds using kNNDM
  knndm_folds <- knndm(
    tpoints = training_data,          # Training points (sf object)
    modeldomain = raster_data,       # Prediction area (SpatRaster)
    k = k                            # Number of folds
  )
  
  # Extract folds for cross-validation
  folds <- knndm_folds$indx_train
  
  # Define predictor variables by excluding the label column
  predictors <- colnames(extr_data)[colnames(extr_data) != "Label"]
  
  # Train Random Forest model
  model <- train(
    x = extr_data[, predictors],     # Predictor variables
    y = extr_data$Label,             # Target variable
    method = "rf",                   # Random Forest algorithm
    importance = TRUE,               # Calculate variable importance
    ntree = ntree,                   # Number of trees in the forest
    trControl = trainControl(
      method = "cv",                 # Cross-validation method
      index = folds,                 # Folds from kNNDM
      savePredictions = "final"      # Save final predictions
    )
  )
  
  # Return trained model
  return(model)
}



# Beispielaufruf der Funktion
fiji_10x10_RF_model <- train_model(train_fiji_extr, train_fiji, num_clusters = 10)
rlp_10x10_RF_model <- train_model(train_rlp_extr, train_rlp, num_clusters = 10)


# repeat with modified data
# ...

# save the created models
saveRDS(fiji_10x10_RF_model,file="path/to/fiji_10x10_RF_model.RDS") 
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
writeRaster(prediction_fiji_10x10, "path/to/prediction_fiji_10x10.tif", overwrite = TRUE)
writeRaster(prediction_rlp_fiji_10x10, "path/to/prediction_rlp_fiji_10x10.tif", overwrite = TRUE)

# repeat the other way around
rlp_boundary <- st_read("path/to/rlp_bord.geojson", quiet = TRUE)
prediction_rlp_10x10_cropped <- crop(prediction_rlp_10x10, vect(rlp_boundary))
prediction_rlp_10x10_mask <- mask(prediction_rlp_10x10_cropped, vect(rlp_boundary))
prediction_rlp_10x10 <- prediction_rlp_10x10_mask
writeRaster(prediction_rlp_10x10, "path/to/prediction_rlp_10x10.tif", overwrite = TRUE)

# repeat for RLP and GERM
# ...



#################
# evaluate model accuracy for stratified models according to Pflugmacher, 2024

#' evaluate_model_accuracy:
#' This function processes spatial training and prediction data for stratified accuracy assessment.
#' @param prediction_raster A raster object containing predicted classes (e.g., multispectral prediction raster).
#' @param training_points A spatial vector containing training points with strata and reference labels.
#' @param strata_ids A vector of unique strata IDs (e.g., 1:8).
#' @param stratum_pixel_counts A numeric vector containing the number of pixels per stratum.
#' @param seed A seed for reproducibility (default is 321).
#'
#' @return A list containing stratified accuracy statistics, overall accuracy, User's and Producer's accuracy, and estimated area proportions per class.
#' @export
#'
#' @examples
#' accuracy_results <- evaluate_model_accuracy(prediction_raster = prediction_fiji_10x10, 
#'                                             training_points = train_fiji, 
#'                                             strata_ids = 1:8, 
#'                                             stratum_pixel_counts = c(1779768, 3549325, 541204, 687659, 14279258, 15115599, 4972515, 116131948))
evaluate_model_accuracy <- function(prediction_raster, training_points, strata_ids, stratum_pixel_counts, seed = 456) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Convert training points to spatial vector and filter by raster extent
  training_vect <- vect(training_points)
  raster_extent <- ext(prediction_raster)
  filtered_points <- crop(training_vect, raster_extent)
  
  # Extract predictions for filtered training points
  predictions <- terra::extract(prediction_raster, filtered_points, ID = FALSE)
  
  # Extract strata, reference labels, and predicted classes
  strata <- filtered_points$strata
  reference <- filtered_points$Label
  map <- predictions$class
  
  # Ensure class consistency between reference and map
  all_classes <- union(unique(reference), unique(map))
  reference <- factor(reference, levels = all_classes)
  map <- factor(map, levels = all_classes)
  
  # Validate lengths of strata, reference, and map
  stopifnot(length(strata) == length(reference))
  stopifnot(length(strata) == length(map))
  
  # Perform stratified accuracy assessment
  accuracy_stats <- aa_stratified(
    stratum = strata,                # Stratum IDs
    reference = reference,           # Reference classes
    map = map,                       # Predicted classes
    h = strata_ids,                  # Unique strata IDs
    N_h = stratum_pixel_counts       # Pixel counts per stratum
  )
  
  # Return stratified accuracy results
  return(accuracy_stats)
}

accuracy_results <- evaluate_model_accuracy(prediction_raster = prediction_fiji_10x10, 
                                            training_points = train_fiji, 
                                            strata_ids = 1:8, 
                                            stratum_pixel_counts = c(1779768, 3549325, 541204, 687659, 14279258, 15115599, 4972515, 116131948))

# Save the accuracy results
saveRDS(accuracy_fiji, file = "path/to/accuracy_fiji_results.RDS")



#################
# set prediction color scheme for better plotting
prediction_colorscheme_fiji <- c("Agriculture" = "bisque1", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                 "Rock" = "azure4", "Shrubland" = "darkgoldenrod2", "Tree" = "darkgreen", 
                                 "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp <- c("Agriculture" = "bisque1","Urban" = "red","Vegetation" = "darkgreen","Water" = "cyan")

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
tmap_save(map, "prediction_fiji_10x10.png")

# repeat for RLP and GERM
# ...

##################### repeat with modified data ######################
# ...

# set prediction color scheme for better plotting
prediction_colorscheme_fiji_modified <- c("Agriculture" = "bisque1", "Grassland" = "green", "Mangrove" = "royalblue2", 
                                          "Rock" = "azure4", "Shrubland" = "darkgoldenrod2", "Urban" = "red", "Water" = "cyan")

prediction_colorscheme_rlp_modified <- c("Agriculture" = "bisque1","Urban" = "red","Water" = "cyan")

# ...
