library(raster)
bat_ensemble<- raster("binomial_icar_pred_bat.tif")
rat_ensemble<- raster("binomial_icar_pred_rat.tif")
bat_coords <- read.csv("bat_X.csv")
coordinates(bat_coords)<-~x+y
rat_coords <- read.csv("rat_X.csv")
coordinates(rat_coords)<-~x+y
# Extract values for these coordinates
bat_coords$bat_extracted_values <- extract(bat_ensemble, bat_coords)
rat_coords$rat_extracted_values <- extract(rat_ensemble, rat_coords)
# Generate random points for comparison
set.seed(123)
library(dismo)
# Number of random points
num_points <- 1000 

# Generate random points while avoiding NA values
random_pts <- randomPoints(bat_ensemble, num_points)
random_df <- as.data.frame(random_pts)
random_risk_values <- raster::extract(bat_ensemble, random_pts)
random_risk_values <- raster::extract(rat_ensemble, random_pts)

# Wilcoxon rank-sum test (Comparing outbreak locations vs random locations)
wilcox.test(bat_coords$bat_extracted_values, random_risk_values, alternative = "greater")
wilcox.test(rat_coords$rat_extracted_values, random_risk_values, alternative = "greater")
