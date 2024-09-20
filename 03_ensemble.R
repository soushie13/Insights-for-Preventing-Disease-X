##================================================================================
##
## 5. Ensemble model - binomial regression with hierarchical Bayesian framework
##
##================================================================================
# stack meta-covariates
prf <- raster("rf_pred_bat.tif")
pbrt <- raster("brt_pred_bat.tif")
prob.p.bi <- raster("Binomial_pred_bats.tif")

prf <- raster("rf_pred_rat.tif")
pbrt <- raster("brt_pred_rat.tif")
prob.p.bi <- raster("Binomial_pred_rats.tif")
preds= stack(prf, pbrt,prob.p.bi)
names(preds) = c("RandomForest","BRT", "biniomal")

#preds=covars
#extract data
obs.data <- read.csv("bat_X.csv", header = TRUE)
obs.data <- read.csv("rat_X.csv", header = TRUE)
coordinates(obs.data)<-~x+y
p <- obs.data
presence = p@coords
presvals <- raster::extract(preds, presence, cellnumber=TRUE)
prevals_coords =cbind(presence,presvals)
##create the pseudo-absences##
pop_d <- raster("pop_density.tif")
Sys.setenv(GITHUB_PAT = "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")

GITHUB_PAT="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"
remotes::install_version("rgeos", version = "0.6-4")
remotes::install_version("rgdal", version = "1.6-7")

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

## Extract environmental values and cell number for observations
pa$Presences <- pa$pb
pa$Trials <- c(1)

## omit rows with missing data
pa.cc=na.omit(pa)


## Normalized continuous covariates
pa.norm <- pa.cc
Mean <- vector()
Sd <- vector()
for (i in c(5:7)) {
  m <- mean(pa.cc[,i],na.rm=TRUE)
  s <- sd(pa.cc[,i],na.rm=TRUE)
  Mean <- c(Mean,m)
  Sd <- c(Sd,s)
  pa.norm[,i] <- (pa.cc[,i]-m)/s
}
## Data-frame with mean and sd for each variable
df.mean.sd <- as.data.frame(rbind(Mean,Sd))
names(df.mean.sd) <- names(pa.norm)[c(5:7)]

## Raster stack for predictions (with normalized covariates)
env <- preds

for (i in c(5:7)) {
  var.name <- names(pa.norm)[i] ## Variable name
  w <- which(names(env)==var.name) ## Position in the stack 
  m <- df.mean.sd[1,var.name] ## Mean
  s <- df.mean.sd[2,var.name] ## Sd
  orig <- values(subset(env,w)) ## Original values
  trans <- (orig-m)/s ## Transformed values
  env[[w]][] <- trans
}

## Select only grid cells with no NA
env.df.pred <- as.matrix(env)
w <- complete.cases(env.df.pred) ## Note: w will be used to obtain the cell identifier for predictions in iCAR model
env.df.pred.complete <- as.data.frame(env.df.pred[w,])

## Make a cluster for parallel MCMCs
nchains <- 2
ncores <- nchains ## One core for each MCMC chains
cores<-detectcores()
clust <- makeCluster(2)
registerDoParallel(clust)

## Starting values and random seed
seed <- 1234
set.seed(seed)
beta.start <- runif(nchains,-1,1)
gamma.start <- runif(nchains,-1,1)
Vrho.start <- runif(nchains,0,10)
seed.mcmc <- round(runif(nchains,0,1e6))
pa.norm$Trials <- c(1)
pa.norm$Presences <- pa.norm$pb




##=============================================================================
##
## 5.2. Binomial with iCAR (spatial autocorrelation)  
##=============================================================================


## Landscape and neighbors
ncells <- ncell(preds)
neighbors.mat <- adjacent(preds, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]
cells.pred <- which(w) ## Vector w indicates the cells with environmental information (without NA)

## binomial icar model
## hSDM model using Binomial icar for perfect detection
mod.binomial.icar <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.binomial.iCAR(presences=pa.norm$Presences,
                            trials=pa.norm$Trials,
                            suitability=~ RandomForest + BRT + biniomal,
                            data=pa.norm,
                            ## Spatial structure
                            spatial.entity=pa.norm$cells,
                            n.neighbors=n.neighbors,
                            neighbors=adj,
                            suitability.pred=env.df.pred.complete,
                            spatial.entity.pred=cells.pred,
                            burnin=1000,
                            mcmc=1000, thin=5,
                            beta.start=beta.start[i],
                            Vrho.start=Vrho.start[i],
                            ## Priors
                            priorVrho="Uniform",
                            mubeta=0, Vbeta=1.0E6,
                            Vrho.max=10,
                            seed=seed.mcmc[i], verbose=1,
                            save.p=1) ## save the post. distributions for each pixel to calculate uncertainty
  return(mod)
}

## Extract list of MCMCs from output
binomial.icar.mcmc <- mcmc.list(lapply(mod.binomial.icar,"[[","mcmc"))

## Outputs summary
bionomial.icar.stat <- summary(binomial.icar.mcmc)$statistics
sink(file="binomial.icar_mcmc_summary_bat.txt")
summary(binomial.icar.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.icar.mcmc)
sink()

bionomial.icar.stat <- summary(binomial.icar.mcmc)$statistics
sink(file="binomial.icar_mcmc_summary_rat.txt")
summary(binomial.icar.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.icar.mcmc)
sink()
## Deviance
deviance.bionomial.icar <- bionomial.icar.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("bionomial.icar_mcmc_trace_bat.pdf")
plot(binomial.icar.mcmc)
dev.off()
pdf("bionomial.icar_mcmc_trace_rat.pdf")
plot(binomial.icar.mcmc)
dev.off()
## Spatial random effects
rho <- subset(preds,1) ## create a raster
values(rho) <- mod.binomial.icar[[1]]$rho.pred
pdf(file="binomial.iCAR_random_effects_bat.pdf")
plot(rho)
dev.off()

rho <- subset(preds,1) ## create a raster
values(rho) <- mod.binomial.icar[[1]]$rho.pred
pdf(file="binomial.iCAR_random_effects_rat.pdf")
plot(rho)
dev.off()
## Prediction on the landscape
prob.p.b <- subset(preds,1) ## create a raster for predictions
values(prob.p.b)[w] <- mod.binomial.icar[[1]]$theta.pred
values(prob.p.b)[w] <- apply(mod.binomial.icar[[1]]$theta.pred,2,mean) ## assign predicted values
values(prob.p.b)[w] <- mod.binomial.icar[[1]]$theta.pred
values(prob.p.b)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
plot(prob.p.b)
plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
plot(pa.norm[pa.norm$pb>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p.b,filename="binomial_icar_pred_bat.tif",overwrite=TRUE)

writeRaster(prob.p.b,filename="binomial_icar_pred_rat.tif",overwrite=TRUE)

#Matrix with CI for each pixel
#prob.p.quant <- apply(mod.binomial.icar[[1]]$theta.pred,2,quantile, c(0.025,0.975))
#prob.p.mean <- apply(mod.binomial.icar[[1]]$theta.pred,2,mean)
prob.p.sd <- apply(mod.binomial.icar[[1]]$theta.pred,2,sd)

#map uncertainity
prob.p.stdev <- subset(preds,1)
#assign values
values(prob.p.stdev)[w]<- prob.p.sd
uncertainity <- prob.p.stdev
plot(uncertainity)

writeRaster(uncertainity,filename="uncertainity_bat.tif",overwrite=TRUE)

writeRaster(uncertainity,filename="uncertainity_rat.tif",overwrite=TRUE)
