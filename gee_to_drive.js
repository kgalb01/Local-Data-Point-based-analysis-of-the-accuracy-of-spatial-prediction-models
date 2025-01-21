/*
* Local Data Point Density-based analysis of spatial prediction models
* Author: Kieran Galbraith
* Date: 2025-01-20
* Description: JS Script to download Sentinel-2 data from Google Earth Engine and export to Google Drive
* NOTE: Only works in Google Earth Engine Code Editor
*/

// Define the area of interest (AOI) for Rheinland-Pfalz or Fiji
var aoi = geometry

// Function to mask clouds using the Sentinel-2 Scene Classification Layer (SCL) band
function cloudMask(image) {
    var scl = image.select('SCL');
    var mask = scl.eq(3).or(scl.gte(7).and(scl.lte(11)));
    return image.updateMask(mask.not());
}

// Retrieve the Sentinel-2 image collection for the specified date range and filter by cloud percentage
var image = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterDate('2023-05-01', '2023-10-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .filterBounds(aoi)
    .map(cloudMask)
    .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12'])
    .median();

// Retrieve another Sentinel-2 image collection for a different date range and filter by cloud percentage
// to fill the gaps in the first image
var image2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate('2016-01-01', '2022-12-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1.0))
    .map(cloudMask)
    .select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12'])
    .median()
    .clip(aoi);

var combined_img = ee.ImageCollection([image2, image]).mosaic();

// Set visualization parameter for Sentinel-2 image
var visParamsTrue = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3500, gamma: 1.1};
Map.addLayer(combined_img.clip(aoi), visParamsTrue, "Rheinland-Pfalz"); // might need to change the name for clarity
Map.centerObject(aoi, 8);

// Export the image to Google Drive
Export.image.toDrive({
  image: combined_img,
  description: 'Sen2RLP', // might need to change the name for clarity
  scale: 10,
  region: aoi,
  maxPixels: 1e13
});
