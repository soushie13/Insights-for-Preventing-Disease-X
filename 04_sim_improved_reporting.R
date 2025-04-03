prf <- raster::raster("rf_pred_bat.tif")
pbrt <- raster::raster("brt_pred_bat.tif")
prob.p.bi <- raster::raster("Binomial_pred_bats.tif")

prf <- raster::raster("rf_pred_rat.tif")
pbrt <- raster::raster("brt_pred_rat.tif")
prob.p.bi <- raster::raster("Binomial_pred_rats.tif")
preds= stack(prf, pbrt,prob.p.bi)
names(preds) = c("RandomForest","BRT", "biniomal")

#preds=covars
#extract data
obs.data <- read.csv("bat_X.csv", header = TRUE)
obs.data <- read.csv("rat_X.csv", header = TRUE)
obs.data <- na.omit(obs.data)
coordinates(obs.data)<-~x+y
p <- obs.data
presence = p@coords
presvals <- extract(preds, presence, cellnumber=TRUE)
prevals_coords =cbind(presence,presvals)
##create the pseudo-absences##
pop_d <- raster("pop_density.tif")
Sys.setenv(GITHUB_PAT = "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")

Sys.setenv(GEOS_CONFIG = "/opt/homebrew/bin/geos-config")
remotes::install_version("rgeos", type = "source", version = "0.6-4", configure.args = "--with-geos-config=/opt/homebrew/bin/geos-config")


remotes::install_version("rgdal", type = "source",version = "1.6-7",configure.args = c("--with-gdal-config=/opt/homebrew/bin/gdal-config", 
                                                                                       "--with-proj-lib=/opt/homebrew/lib", 
                                                                                       "--with-proj-include=/opt/homebrew/include"))

 remotes::install_github("SEEG-Oxford/seegSDM")
set.seed(2024)
bg <- bgSample(pop_d,
               n = 1000,
               prob = TRUE,
               replace = TRUE,
               spatial = FALSE)

colnames(bg) <- c('x', 'y')
absence<- data.frame(bg)
absvals <- raster::extract(preds, absence, cellnumber=TRUE)
absvals_coords =cbind(absence,absvals)
pb <- c(rep(1, nrow(prevals_coords)), rep(0, nrow(absvals_coords)))
pa <- data.frame(cbind(pb, rbind(prevals_coords, absvals_coords)))

## Extract predicted probability of presence
## Suitable sites
prob.p.b<- raster(".tif")
pa$Suit <- 0
pa$Suit[pa$pb>0] <- 1
pa.obs.data <- pa[, c("x", "y")]

Index.fun <- function(obs,prob.p,thresh) {
  ## Transform probabilities into {0,1} given threshold
  pred <- ifelse(prob.p>=thresh,1,0)
  ## Contingency table (pred/obs)
  n00 <- sum(pred==0 & obs==0,na.rm=TRUE)
  n11 <- sum(pred==1 & obs==1,na.rm=TRUE)
  n01 <- sum(pred==0 & obs==1,na.rm=TRUE)
  n10 <- sum(pred==1 & obs==0,na.rm=TRUE)
  ## Threshold  dependent indexes
  OA <- (n11+n00)/(n11+n10+n00+n01) ## Overall accuracy
  Sensitivity <- n11/(n11+n01)
  Specificity <- n00/(n00+n10)
  TSS <- Sensitivity+Specificity-1
  return(list(OA=OA,TSS=TSS,Sens=Sensitivity,Spe=Specificity))
}
## Extract predicted probability of presence
pa$prob.p <- raster::extract(prob.p.b,pa.obs.data)

## TSS as a function of the threshold
OA <- vector()
TSS <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(pa$Suit,pa$prob.p,thresh.seq[i])
  OA[i] <- Index$OA
  TSS[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS
maxTSS <- max(TSS,na.rm=TRUE)
maxTSS
w.th <- which(TSS==maxTSS)
thresh.maxTSS <- mean(thresh.seq[w.th],na.rm=TRUE)
OA.thresh.maxTSS <- Index.fun(pa$Suit,pa$prob.p,thresh.maxTSS)$OA 
tss.df <- data.frame(masTSS=maxTSS,OA=OA.thresh.maxTSS,prob=thresh.maxTSS)

## Plot evolution of TSS with threshold
pdf(file="TSS.pdf")
plot(thresh.seq,TSS,type="l",xlab=c("probability threshold"),ylab="TSS")
abline(v=thresh.maxTSS)
dev.off()
