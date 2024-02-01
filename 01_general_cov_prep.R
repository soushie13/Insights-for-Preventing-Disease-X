### Disease X modelling- general covariates####
### General environmental and demographic covariates:temparature (max,min), precipitation, population density, 
tmax <- getData(name = "worldclim",var = "tmax", res = 5, path = "covars_10k")
tmax_m <- mean(tmax)
res(tmax)
e <- extent(covars)
tmax_m <- crop(tmax_m, e)
writeRaster(tmax_m, "covars_10k/tmax_m.tif",format = 'GTiff', overwrite = T)
tmin <- getData(name = "worldclim",var = "tmin", res = 5, path = "covars_10k")
tmin_m <- mean(tmin)
tmin_m <- crop(tmin_m, e)
writeRaster(tmin_m, "covars_10k/tmin_m.tif",format = 'GTiff', overwrite = T)
ppt <- getData(name = "worldclim",var = "prec", res = 5, path = "covars_10k")
ppt_m <- mean(ppt)
ppt_m <- crop(ppt_m, e)
writeRaster(ppt_m, "covars_10k/ppt_m.tif",format = 'GTiff', overwrite = T)
clim = stack(tmax_m,tmin_m,ppt_m)
names(clim) = c("tmax","tmin","ppt")

### Deforestation and reforestation data from GFW ###
tiles <- calc_gfc_tiles(test_poly)
plot(tiles)
plot(test_poly, lt=2, add=TRUE)

download_tiles(
  tiles,
  output_folder,
  images = c( "gain"),
  dataset = "GFC-2022-v1.10"
)

download_tiles(
  tiles,
  output_folder,
  images = c( "lossyear"),
  dataset = "GFC-2022-v1.10"
)



### check covariate correlation #
correlation.matrix<-cor(envtrain, method=c("spearman"))

hc <- hclust(as.dist(1-abs(correlation.matrix)),
             method="centroid")
par(mfrow=c(1,1))
plot (hc, sub="", xlab="Variables", hang=-1, ylab="1-
     abs(correlation)")
abline(0.3,0)
rect.hclust(hc,h=0.3)
##Variance Inflation Factor (VIF) 
#load the car library
library(car)
model <- lm(y ~x1+  x2+x3+x4
            + x5 + x6 +x7+x8+x9, data = envtrain)

#view the output of the regression model
summary(model)
#calculate the VIF for each predictor variable in the model
vif(model)
#create vector of VIF values
vif_values <- vif(model)