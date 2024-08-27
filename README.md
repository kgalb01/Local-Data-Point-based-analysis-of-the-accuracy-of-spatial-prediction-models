# Local Data Point-based analysis of the accuracy of spatial prediction models 
This work aims to evaluate the accuracy of spatial prediction models using the "Local Data Point Density" (F. Schumacher et al., 2024) and therefore proving the functionality of the "LPD" and it's advantages as a extension of the "Disimilarity Index" (H. Meyer, E. Pebesma, 2021). The purpose of this GitHub Page is to grant further insides into the methodology of the study. It is supposed to give a quick overview of the project, summarising the results. The uploaded code is slightly modified and for the sake of reproducibility every used function is seeded. In order to improve readability of the code, some parts will be removed (e.g. where the same command is repeated several times to get the same result for each case study area). For further questions feel free to contact me: k_galb01@uni-muenster.de

NOTE: At this point the bibliography is incomplete and it's supposed to be completed once the study is done!

## Motivation
- Machine learning is now a mainstream technology, used beyond just scientific applications.
- In geoscience, it helps predict spatial processes and calculate land use classifications.
- The "Area of Applicability" (AOA) determines where models can be applied reliably, by setting a threshold for the "Dissimilarity Index" (DI) (Source).
- The DI calculates the Euclidian distance between a new data point and a training data point, assuming that for more similarities between a model is to the new data, the model is able to create more reliable predictions (Source).
- The "Local Data Point Density" (LPD) improves AOA by considering both distance and data density, enhancing model accuracy in clustered areas (Source).
- This study evaluates LPD's effectiveness in improving spatial prediction models.

## Case Study
- Chosen for their contrasting geographical characteristics to test the LPD approach across diverse landscapes the Fiji Islands and Rhineland-Palatinate act as case study areas.
- One of the reasons for chosing Fiji is the availability of high-quality training data with ground truthing with diverse land cover classes (Source).
- The plan is to be using these existing training data to model land use and land cover (LULC) classifications in R, leveraging the LPD for better accuracy. 
- Rhineland-Palatinate was picked for its geographical contrast to Fiji but at the same time its similar size, which allows for comparable data distribution.
  * The training data for Rhineland-Palatinate was created using a R script that picks n random points out of the official DLM for Rhineland-Palatinate (Source).
- The goal with these training data is to create a model in a different environmental context, testing the LPD’s applicability.
- As a way of validating the results both sets of training data are modified to remove an entire class ('Tree' for Fiji and 'Vegetation for Rhineland-Palatinate).
- At the end the study is supposed to demonstrate LPD’s effectiveness in improving spatial prediction accuracy across different regions by assessing its performance in both tropical and temperate environments.

## Results

## Bibliography
