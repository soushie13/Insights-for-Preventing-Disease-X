#covariates bat
covars <- stack(alt, bats, bushmeat, deforestation, mammals, lulc, pop,  ppt,  tmin,travel)
#covariates rodent 
covars <- stack(alt, rats, crops, deforestation,  lulc,  pop, ppt,  tmin, travel)

plot(covars)
### open presence data##
bat_xy <- read.csv("~/Documents/GitHub/Disease_X/bat_X.csv", header = TRUE)
obs.data <- bat_xy[, c("x", "y")]
plot(obs.data, pch = 20)

rat_xy <- read.csv("~/Documents/GitHub/Disease_X/rat_X.csv", header = TRUE)
obs.data <- rat_xy[, c("x", "y")]
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
envtrain <- raster::extract(covars, train)
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
model <- lm(pa ~ alt+ rodent_density+treeloss+  lulc
            + pop+ ppt_m+ tmin_m+ travel+cropland, data = envtrain)

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
