output_folder <- "~/Documents/GitHub/Disease_X/covar_10k/bat_covars"

tmax <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/tmax_m.tif")
tmin <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/tmin_m.tif")
ppt <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/ppt_m.tif")
bushmeat <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/dist_BA_10k.tif")
bushmeat <- resample(bushmeat, tmin, method= "bilinear") 
writeRaster(bushmeat, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/bushmeat_global.tif",format = 'GTiff', overwrite = T)
bushmeat <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/bushmeat_global.tif")

bushmeat[is.na(bushmeat)] <- 0.00001
plot(bushmeat_global)
# Create a global extent raster
global_extent <- extent(-180, 180, -90, 90)
global_raster <- raster(global_extent, res = res(bushmeat), crs = crs(bushmeat))
values(global_raster) <- 0  # Set all values to 0

# Resample your original raster to fit the global extent
bushmeat_global <- resample(bushmeat, global_raster, method = "bilinear")
bushmeat <- resample(bushmeat_global, tmin, method= "bilinear") 

deforestation <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/deforestation.tif")
deforestation <- resample(deforestation, tmin, method= "bilinear") 
writeRaster(deforestation, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/treeloss.tif",format = 'GTiff', overwrite = T)
deforestation <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/treeloss.tif")


lulc <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/lc_modi.tif")
lulc <- resample(lulc, tmin, method= "bilinear") 
writeRaster(lulc, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/lulc.tif",format = 'GTiff', overwrite = T)
lulc <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/lulc.tif")


pop <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/pop.tif")
pop <- resample(pop, tmin, method= "bilinear") 
writeRaster(pop, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/pop.tif",format = 'GTiff', overwrite = T)
pop <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/pop.tif")


travel <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/travel_time_to_cities_5.tif")
travel <- resample(travel, tmin, method= "bilinear") 
writeRaster(travel, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/travel.tif",format = 'GTiff', overwrite = T)
travel <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/travel.tif")


bats <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/bats_merged.tif")
bats <- resample(bats, tmin, method= "ngb") 
writeRaster(bats, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/bats.tif",format = 'GTiff', overwrite = T)
bats <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/bats.tif")


mammals <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/all_mammals.tif")
mammals <- resample(mammals, tmin, method= "bilinear") 
writeRaster(mammals, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/mammals.tif",format = 'GTiff', overwrite = T)
mammals <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/mammals.tif")


lc <- raster("~/Documents/GitHub/Disease_X/covar_10k/lc.tif")
lc <- resample(lc, tmin, method= "ngb") 
writeRaster(lc, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/lc.tif",format = 'GTiff', overwrite = T)

alt <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/wc2.1_5m_elev.tif")
alt <- resample(alt, tmin, method= "bilinear") 
writeRaster(alt, "~/Documents/GitHub/Disease_X/covar_10k/bat_covars/alt.tif",format = 'GTiff', overwrite = T)
alt <- raster("~/Documents/GitHub/Disease_X/covar_10k/bat_covars/alt.tif")


covars <- stack(alt, bats, bushmeat, deforestation, lc, lulc, mammals, pop, ppt,  tmin, travel)
plot(covars)
### open presence data##
bat_xy <- read.csv("~/Documents/GitHub/Disease_X/bat_X.csv", header = TRUE)
obs.data <- bat_xy[, c("x", "y")]
plot(obs.data, pch = 20)

### define study extent ##
e <- extent(covars)
##create the background##
data(wrld_simpl)
bw<- wrld_simpl[wrld_simpl@data$UN!="10",]
## training and testing data ##
set.seed(0)
k<-4
group <- kfold(obs.data, 4)

pres_train <- obs.data[group != 1, ]
pres_test <- obs.data[group == 1, ]

##create the pseudo-absences##
set.seed(10)
background <- spsample(bw,n= 1000,"random", iter= 5000) 
plot(background)
background <- background@coords
group <- kfold(background, 4)
backg_train <- background[group != 1, ]
backg_test <- background[group == 1, ]

#extract data
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(covars, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain[is.na(envtrain)] <- 0.00000001

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
model <- lm(pa ~ wc2.1_5m_elev+ bats_merged+dist_BA_10k+deforestation+ lc+ lc_modi
            + all_mammals+ ppp_2020_1km_Aggregated+ ppt_m+ tmin_m+ travel_time_to_cities_5, data = envtrain)

#view the output of the regression model
summary(model)


#calculate the VIF for each predictor variable in the model
vif(model)
#create vector of VIF values
vif_values <- vif(model)

#create horizontal bar chart to display each VIF value
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")

#add vertical line at 5
abline(v = 10, lwd = 3, lty = 2)

covars <- dropLayer(covars, tmax)
stackSave(covars, "mystack")
