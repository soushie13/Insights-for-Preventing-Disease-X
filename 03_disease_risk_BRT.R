##===============================================
##
## Boosted regression trees - Gradient Boosted Machine (GBM)
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
pbrt <- predict(brtfinal_model, data = covars, type = "cloglog")

##plot
plot(pbrt, main='BRT prediction')
##save prediction raster
writeRaster(pbrt, "brt_pred.tif",format = 'GTiff', overwrite = T)
plotROC(brtfinal_model)
modelReport(brtfinal_model,folder = "BRT",type = "cloglog", test = test,
            response_curves = TRUE,jk = TRUE, permut = 10)

varImp(brtfinal_model)
