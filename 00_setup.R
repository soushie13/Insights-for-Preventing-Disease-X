
##covariates
library(geodata) #download worldclim data
library(gfcanalysis) #download reforestation and deforestation data
##spatial mapping and plotting packages
library(ggplot2)
library(rasterVis)
library(maptools)
library(maps)
library(terra)
library(dplyr)
library(raster)
library(rgeos)
library(rgdal)
library(sf)
library(tmap)
library(plyr)
library(doParallel)
library(foreach)
library(coda)
##spatial modeling packages
library(seegSDM)      #generation of background points biased on population density
library(dismo)
library(rJava)        #MaxENT
library(SDMtune)
library(zeallot)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(kableExtra)   # Compile ROC reports
library(blockCV)      # spatial cross-validation
library(randomForest)
library(embarcadero)
detach("package:embarcadero", unload=TRUE)
library(hSDM)
