# Local Data Point Density-based analysis of spatial prediction models
# Author: Kieran Galbraith
# Date: 2025-01-20
# Description: R Script to convert .CSV file containing the extraction data from Google Earth Engine to .RDS File


# Load necessary libraries
library(readr)

# Define the path to the .CSV file
csv_file_path_fiji <- "path/to/fiji_training_extr.csv"

csv_file_path_rlp <- "path/to/rlp_training_extr.csv"


# repeat with modified data
csv_file_path_fiji_modified <- "path/to/fiji_training_extr_modified.csv"

csv_file_path_rlp_modified <- "path/to/rlp_training_extr_modified.csv"


# Define the path to the .RDS file
rds_file_path_fiji <- "path/to/train_fiji_extr.rds"

rds_file_path_rlp <- "path/to/train_rlp_extr.rds"



rds_file_path_fiji_modified <- "path/to/train_fiji_extr_modified.rds"

rds_file_path_rlp_modified <- "path/to/train_rlp_extr_modified.rds"




# Read the .CSV file
training_data_fiji <- read_csv(csv_file_path_fiji)
training_data_rlp <- read_csv(csv_file_path_rlp)


training_data_fiji_modified <- read_csv(csv_file_path_fiji_modified)
training_data_rlp_modified <- read_csv(csv_file_path_rlp_modified)

# Save the data as a .RDS file
saveRDS(training_data_fiji, rds_file_path_fiji)
saveRDS(training_data_rlp, rds_file_path_rlp)

saveRDS(training_data_fiji_modified, rds_file_path_fiji_modified)
saveRDS(training_data_rlp_modified, rds_file_path_rlp_modified)
