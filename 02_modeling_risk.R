##===============================================
##
## 1. Random Forest
##
##===============================================

###RF training
#extract data
coordinates(obs.data)<-~x+y
p <- obs.data
presence = p@coords
##create the pseudo-absences##
pop_d <- raster("pop_density.tif")
remotes::install_github("SEEG-Oxford/seegSDM")
set.seed(10)
bg <- bgSample(pop_d,
               n = 1000,
               prob = TRUE,
               replace = TRUE,
               spatial = FALSE)

colnames(bg) <- c('x', 'y')
absence<- data.frame(bg)


### random forest  with spatial cross-validation
predictors <- terra::rast(covars)
# Create SWD object
swd_data <- prepareSWD(species = "bat_pathogens", p = presence, a = absence,
                       env = predictors)
# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(swd_data, test = 0.2, seed = 25)
train <- datasets[[1]]
test <- datasets[[2]]

rf_model <- train(method = "RF", data = train, ntree=1000, nodesize=10, importance =T)
auc(rf_model)
tss(rf_model)
library(plotROC)
plotROC(rf_model)
map <- predict(rf_model, data = predictors, type = "cloglog")
plot(map)


c(train, test) %<-% trainValTest(swd_data, test = 0.2, only_presence = FALSE, seed = 25)
rf_model <- train(method = "RF", data = train, ntree=1000, nodesize=10, importance =T)
cat("Training auc: ", auc(rf_model))
cat("Testing auc: ", auc(rf_model, test = test))

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(25)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = seeds[i]) # Make the train/test split
  m <- train("RF", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
# Print the output
output
# compute the range of the testing AUC
range(output[, 3])
#Cross Validation
folds <- randomFolds(swd_data, k = 4,  seed = 25)
cv_model <- train("RF", data = swd_data, folds = folds)
cv_model
cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))
cat("Training TSS: ", tss(cv_model))
cat("Testing TSS: ", tss(cv_model, test = TRUE))
#Spatial Cross Validation
# Create spatial points data frame
sp_df <- SpatialPointsDataFrame(swd_data@coords, data = as.data.frame(swd_data@pa), proj4string = crs(covars))
e_folds <- envBlock(rasterLayer = covars, speciesData = sp_df, species = "swd_data@pa", k = 4, standardization = "standard", rasterBlock = FALSE, numLimit = 100)
scv_model <- train(method = "RF", data = swd_data, fc = "l", reg = 0.8, folds = e_folds)
cat("Training AUC: ", auc(scv_model))
cat("Testing AUC: ", auc(scv_model, test = TRUE))
cat("Training TSS: ", tss(scv_model))
cat("Testing TSS: ", tss(scv_model, test = TRUE))
plotResponse(cv_model, var = "deforestation", type = "cloglog", marginal = TRUE, fun = mean, rug = TRUE)
#tune the model hyperparameters
c(train, val, test) %<-% trainValTest(swd_data, val = 0.2, test = 0.2,  seed = 61516)
cat("# Training  : ", nrow(train@data))
cat("# Validation: ", nrow(val@data))
cat("# Testing   : ", nrow(test@data))
rf_model <- train("RF", data = train)

getTunableArgs(rf_model)
#[1] "mtry"     "ntree"    "nodesize"
# Define the values for bg
h <- list(ntree=seq(500))
# Call the gridSearch function
exp_1 <- gridSearch(rf_model, hypers = h, metric = "auc", test = val)
exp_1 <- gridSearch(rf_model, hypers = h, metric = "tss", test = val)
#Random search
h <- list(nodesize = seq(10), ntree = seq(500))
exp_6 <- randomSearch(rf_model, hypers = h, metric = "auc", test = val, pop = 10, seed = 65466)
exp_6 <- randomSearch(rf_model, hypers = h, metric = "tss", test = val, pop = 10, seed = 65466)
exp_6@results
exp_7 <- optimizeModel(rf_model, hypers = h, metric = "tss", test = val, pop = 15, gen = 2, keep_best = 0.4, keep_random = 0.2, mutation_chance = 0.4, seed = 798)

## fit final random forest model
#merge training and validation data
merged_data <- mergeSWD(train, val)
set.seed(2024)
RFfinal_model <- train("RF", data = merged_data, ntree = 500, nodesize = 4)
auc(RFfinal_model, test = test)
tss(RFfinal_model, test = test)

#predict distribution
pr <- predict(RFfinal_model, data = predictors, type = "cloglog")
##plot
par(mfrow=c(1,1))
plot(pr, main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
##save prediction raster
pr_bat_pop <- raster::raster(pr)
writeRaster(pr_bat_pop, "rf_pred_bat_pop.tif",format = 'GTiff', overwrite = T)
plotROC(RFfinal_model)
modelReport(RFfinal_model,type = "cloglog", test = test,
            jk = TRUE, permut = 10, folder = "RF")

varImp(RFfinal_model)
##===============================================
##
## 2. Boosted regression trees - Gradient Boosted Machine (GBM)
##
##===============================================

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(swd_data, test = 0.2,  seed = 345)
train <- datasets[[1]]
test <- datasets[[2]]

brt_model <- train(method = "BRT", data = train, family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
auc(brt_model)
tss(brt_model)


c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = 345)
brt_model <- train(method = "BRT", data = train, family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
cat("Training auc: ", auc(brt_model))
cat("Testing auc: ", auc(brt_model, test = test))

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(345)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = seeds[i]) # Make the train/test split
  m <- train("BRT", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
# Print the output
output
# compute the range of the testing AUC
range(output[, 3])
#Cross Validation
folds <- randomFolds(swd_data, k = 4,  seed = 345)
cv_model <- train("BRT", data = swd_data, folds = folds, family = "bernoulli", tree.complexity = 5,
                  learning.rate = 0.01, bag.fraction = 0.5)
cv_model
cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))
cat("Training TSS: ", tss(cv_model))
cat("Testing TSS: ", tss(cv_model, test = TRUE))
#Spatial Cross Validation
# Create spatial points data frame
sp_df <- SpatialPointsDataFrame(swd_data@coords, data = as.data.frame(swd_data@pa), proj4string = crs(covars))
e_folds <- envBlock(rasterLayer = covars, speciesData = sp_df, species = "swd_data@pa", k = 4, standardization = "standard", rasterBlock = FALSE, numLimit = 100)
scv_model <- train(method = "BRT", data = swd_data, folds = e_folds,family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
cat("Training AUC: ", auc(scv_model))
cat("Testing AUC: ", auc(scv_model, test = TRUE))
cat("Training TSS: ", tss(scv_model))
cat("Testing TSS: ", tss(scv_model, test = TRUE))
#tune the model hyperparameters
c(train, val, test) %<-% trainValTest(swd_data, val = 0.2, test = 0.2,  seed = 345)
cat("# Training  : ", nrow(train@data))
cat("# Validation: ", nrow(val@data))
cat("# Testing   : ", nrow(test@data))
brt_model <- train("BRT", data = train,family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.0025, bag.fraction = 0.5)

getTunableArgs(brt_model)
# "distribution"      "n.trees"           "interaction.depth" "shrinkage"   "bag.fraction"
# Define the values for bg
h <- list(n.trees=seq(500))
# Call the gridSearch function
exp_1 <- gridSearch(brt_model, hypers = h, metric = "auc", test = val)
exp_1 <- gridSearch(brt_model, hypers = h, metric = "tss", test = val)
exp_1@results
#Random search
h <- list(n.trees=seq(500), shrinkage=seq(0.5))
exp_6 <- randomSearch(brt_model, hypers = h, metric = "auc", test = val, pop = 10, seed = 345)
exp_6 <- randomSearch(brt_model, hypers = h, metric = "tss", test = val, pop = 10, seed = 345)
exp_6@results

## fit final BRT model
#merge training and validation data
merged_data <- mergeSWD(train, val)
set.seed(2022)
brtfinal_model <- train("BRT", data = merged_data, n.trees=500, learning.rate = 0.0025, bag.fraction = 0.5)
auc(brtfinal_model, test = test)
tss(brtfinal_model, test = test)
#predict distribution
pbrt <- predict(brtfinal_model, data = predictors, type = "cloglog")

##plot
plot(pbrt, main='BRT prediction')
##save prediction raster
pbrt_bat <- raster::raster(pbrt)
writeRaster(pbrt_bat, "brt_pred_bat.tif",format = 'GTiff', overwrite = T)
plotROC(brtfinal_model)
modelReport(brtfinal_model,folder = "BRT",type = "cloglog", test = test,
            response_curves = TRUE,jk = TRUE, permut = 10)

varImp(brtfinal_model)

##===============================================
##
## 4. BART
##
##===============================================
remotes::install_github("cjcarlson/embarcadero")
set.seed(12345)
pb_xy <- rbind(presence, absence)
pb_pred <- raster::extract(covars, pb_xy)
pb <- c(rep(1, nrow(presence)), rep(0, nrow(absence)))
data <- data.frame( cbind(pa=pb, pb_pred) )
data <- na.omit(data)
xvars <- names(data)[!names(data)=="pa"]

BART.model <- bart.step(x.data= data[,xvars], 
                        y.data= data[,"pa"], 
                        full=TRUE,
                        quiet = TRUE)
# Saving in RData format
save(BART.model, file = "bart.RData")
summary(BART.model)


bart.map <- predict2.bart(object= BART.model,
                          x.layers = covars,
                          quantiles = c(0.025, 0.975),
                          splitby=20)
plot(bart.map[[1]])
plot(bart.map[[3]])
q <- quantile(values(bart.map[[3]]- bart.map[[2]]), 0.75, na.rm=TRUE)
uncertinity <- (bart.map[[3]]- bart.map[[2]]) > q
sp <- spartial(bart.map, preds_10, x.vars="bushmeat", equal=TRUE)
plot(uncertinity)
writeRaster(bart.map[[1]], "bart_bat.tif" ,format = 'GTiff', overwrite = T)
writeRaster(bart.map[[2]], "bart2.5_bat.tif" ,format = 'GTiff', overwrite = T)
writeRaster(bart.map[[3]], "bart97.5_bat.tif" ,format = 'GTiff', overwrite = T)
writeRaster(uncertinity, "uncertainitybart_bat.tif" ,format = 'GTiff', overwrite = T)


###geographical null model using geographical distance
null_model2 <- convHull(pres_train, lonlat =TRUE)
null_e <- evaluate(null_model2, p= pres_test, a=backg_test)

