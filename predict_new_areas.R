## Script to evaluate PJ spread model in new areas
# N2 - a site close by model training area in central OR dominated by western Juniper
# N3 - a site in Western NV dominated by singleleaf pinyon

# load packages
library(tidyverse)
library(abind) # arrays
library(spatstat) # distance matrix w/ pairwise dist. function
library(raster)
library(cowplot)

# load functions
source("forecastFunctions.R")

# setwd 
setwd("R:/Shriver_Lab/PJspread/evaluate_out_of_sample")

# convert rasters to dataframes
raster2df <- function(raster, na.rm = T){
  vals <- data.frame(raster::getValues(raster))
  cellnum <- raster::cellsFromExtent(raster, raster)
  xy <- raster::xyFromCell(raster, cellnum)
  df <- cbind.data.frame(cellnum, xy, vals)
  if(na.rm == T) df <- df[complete.cases(df),]
  return(df)
}

# Function to re-rasterize dataframe
df2raster <- function(df) { # df should be x, y, vals only
  val <- names(df)[3]
  spg <- df %>% rename(values = val)
  sp::coordinates(spg) <- ~ x + y
  sp::gridded(spg) <- TRUE # coerce to SpatialPixelsDataFrame
  r <- raster::raster(spg, values = TRUE) # coerce to raster
  return(r)
}
# ------------- Observed Data -------------------
# load of observed data for in sample locations
load("R:/Shriver_Lab/PJspread/data_prepping/data_prepped.RData")

## load observed data for new locations
# cover data
N2 <- read.csv("RAPtreecoverData_oos1.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

N3 <- read.csv("RAPtreecoverData_oos2.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

# ---covariate data----

# in-sample un-standardized covariates
N.clim <- read.csv("PrismClimateDataV3.csv")
N.topo <- read.csv("topographic_data.csv")

## N2
N2.clim <- read.csv("PrismClimateDataV3_oos1.csv") %>% 
  
  #standardize covariates based on mean and sd from in-sample data
  dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%
  
  dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%
  
  dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N2.topo <- read.csv("topographic_data_oos1.csv") %>% 
  
  #standardize covariates based on mean and sd from in-sample data
  dplyr::mutate(heatload = (heatload-mean(N.topo$heatload))/sd(N.topo$heatload)) %>%
  
  dplyr::mutate(elev = (elev-mean(N.topo$elev))/sd(N.topo$elev)) 

# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N2.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N2.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N2.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N2.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N2 <- abind(ppt, tmean, heatload, elev, along = 3)

# ## N3
N3.clim <- read.csv("PrismClimateDataV3_oos2.csv") %>%

   #standardize covariates based on mean and sd from in-sample data
   dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%

   dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%

   dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N3.topo <- read.csv("topographic_data_oos2.csv") %>%

  #standardize covariates based on mean and sd from in-sample data
  dplyr::mutate(heatload = (heatload-mean(N.topo$heatload))/sd(N.topo$heatload)) %>%

  dplyr::mutate(elev = (elev-mean(N.topo$elev))/sd(N.topo$elev))

# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N3.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N3.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N3.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N3.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N3 <- abind(ppt, tmean, heatload, elev, along = 3)

# # ** test plotting variation in environmental covariates for each study area **
# par(mfrow = c(2,2))
# matplot(enviro.var[,,1], type = 'l') # ppt
# matplot(enviro.var[,,2], type = 'l') # tmean
# matplot(enviro.var[,,4], type = 'l') # heatload
# matplot(enviro.var[,,5], type = 'l') # elev
# 
# matplot(enviro.var.N2[,,1], type = 'l') # ppt
# matplot(enviro.var.N2[,,2], type = 'l') # tmean
# matplot(enviro.var.N2[,,3], type = 'l') # heatload
# matplot(enviro.var.N2[,,4], type = 'l') # elev
# 
# matplot(enviro.var.N3[,,1], type = 'l') # ppt
# matplot(enviro.var.N3[,,2], type = 'l') # tmean
# matplot(enviro.var.N3[,,3], type = 'l') # heatload
# matplot(enviro.var.N3[,,4], type = 'l') # elev

#------------- Build distance matrix ---------------

# create empty rasters w same extent as each study area
r2 <- raster(nrows = 42, ncols = 56, xmn = -120.5686+0.334, xmx = -120.5535+0.334, ymn = 43.87138-0.16, ymx = 43.8827-0.16) # extent of N2
r3 <- raster(nrows = 42, ncols = 56, xmn = -120.5686+1.834, xmx = -120.5535+1.834, ymn = 43.87138-5.5324, ymx = 43.8827-5.5324) # extent of N3

# function that returns squared pairwise distance matrix for a given raster
create_dist_mat <- function(raster) {
  # to points
  points <- rasterToPoints(r)
  # specify window 
  rowin <- as.owin(as.vector(extent(r)))
  # calculate pairwise distances, and output into a matrix 
  D <- pairdist(as.ppp(X = points, W = rowin))*1000
  # Dsc = D*1000 # rescale the distance matrix so values are close to 1-100
  Dsq = D*D # square D for input into dispersal equation
  return(Dsq)
}

# create distance matrices
Dsq2 <- create_dist_mat(r2)
Dsq3 <- create_dist_mat(r3)

# ------------ Retrieve parameter estimates for each model ------------------

# load model objects generated from:
# 'eval_model_insample.R'
# Latent states and parameter estimates are post burnin and thinning 

mod1 <- readRDS(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/base_model.rds") # base

mod2 <- readRDS(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/topo_model.rds") # topo

mod3 <- readRDS(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/clim_model.rds") # clim

mod4 <- readRDS(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/topoclim_model.rds") # topoclim

# ----------------- Forecast across new locations -------------------

# fixed/single origin forecast
# allow error to propagate through time
forecast_new_loc <- function(obs, pars, covars, Dsq, mod) {
  
  pars <- as.data.frame(pars)
  
  pixels = length(obs[1,])
  years = length(obs[,1])
  
  # empty array to fill
  predOut <- array(NA, c(years, pixels, length(pars$tauOut))) 
  predst <- abind(list(replicate(length(pars$tauOut), t(as.matrix(colMeans(obs[1:5,])))))) 
  predOut[1,,] <- predst[1,,]
  
  # evaluation
  rmseTotOut<-biasOut<-denseOut<-matrix(NA, years, length(pars$tauOut)) # evaluation metrics
 
   # empty matrix for dispersal matrix
  M1 <- matrix(NA, pixels, pixels)
  
  # calculate predicted values
  for(i in 1:length(pars$tauOut)) { 
    
    M <- dispersal(Dsq, pars$tauOut[i]) 
    
    for(p in 1:pixels){
      M1[,p]=M[,p]/sum(M[,p])
    }
    
    Nt <- colMeans(obs[1:5,])# set initial cover value to average cover for first 5 years
    # Nt <- obs[1,] # set initial cover value as actual observed cover
    
    for(t in 2:years) {
      
      if(mod=='base') {
        G <- growthsimple(a=pars$alphaOut[i],b=pars$betaOut[i],nt=as.vector(Nt))
      }
      
      if(mod=='clim') {
        G <- growthClim(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],b0=pars$betaOut[i],X=covars[t,,],nt=as.vector(Nt))
      }
      
      if(mod=='topo') {
        G <- growthTopo(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2[i],X=covars[t,,],nt=as.vector(Nt))
        
      }
      
      if(mod=='topoclim') {
        G <- growthTopoClim(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],a3=pars$alpha3[i],a4=pars$alpha4[i],b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2[i],X=covars[t,,],nt=as.vector(Nt))
      }
      
      # if(mod=='base_q') {
      #   G <- growth(b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2,as.vector(Nt))
      # }

      # if(mod=='clim_q') {
      #   G <- 
      # }
      # 
      # if(mod=='topo_q') {
      #   G <- 
      # }
      # 
      # if(mod=='topoclim_q') {
      #   G <- 
      # }
      
      Nmean <-M1%*%(diag(G)%*%Nt)
      
      Nt <- rnorm(pixels, Nmean, pars$sig.pOut[i])
      
      # Nt <- replace(Nt, Nt<0, 0) # convert negative values to zeroes
      
      predOut[t,,i] <- Nt
      
      rmseTotOut[t,i] <- rmsefunc(pred=Nt, obs=obs[t,])
      # biasOut[i] <- biasfunc(pred=predOut[,,i], obs=obs[2:years,])
      # denseOut[i] <- densefunc(pred=predOut[,,i], obs=obs[2:years,], sig_o=1)
      
    }
    
  if(i %in% seq(1, length(pars$tauOut), 100)) {print(i)} 
  }
  
  
  to_return <- list()
  to_return$predOut <- predOut
  to_return$rmseTotOut <- rmseTotOut
  # to_return$biasOut <- biasOut
  # to_return$denseOut <- denseOut
  
  return(to_return)
}

# # ** plot predOut to see what is going on with models breaking down **
# n2.xy <- read.csv("RAPtreecoverData_oos1.csv") %>% dplyr::filter(year==1986) %>% dplyr::select(c(x,y))
# par(mfrow = c(3,4))
# for (i in 18:29) {
#   weirdyear <- predOut[i,,1]
#   test <- cbind(n2.xy, weirdyear)
#   plot(main = as.character(i), df2raster(test), zlim = c(-20, 30))
# }

## ----- Build forecasts for each area with each model ----
# (save predictions as a workspace file

# # N
# for.base.N <- forecast_new_loc(obs = N, pars = mod1$pars, Dsq = Dsq, mod = 'base')
# for.topo.N <- forecast_new_loc(obs = N, pars = mod2$pars, covars = enviro.var[,,-3], Dsq = Dsq, mod = 'topo')
# for.clim.N <- forecast_new_loc(obs = N, pars = mod3$pars, covars = enviro.var[,,-3], Dsq = Dsq, mod = 'clim')
# for.topoclim.N <- forecast_new_loc(obs = N, pars = mod4$pars, covars = enviro.var[,,-3], Dsq = Dsq, mod = 'topoclim')
# 
# save.image(file = "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/35y_insample_predictions_5y_average_initial.RData")

# # N2
# for.base.N2 <- forecast_new_loc(obs = N2, pars = mod1$pars, Dsq = Dsq2, mod = 'base')
# for.topo.N2 <- forecast_new_loc(obs = N2, pars = mod2$pars, covars = enviro.var.N2, Dsq = Dsq2, mod = 'topo')
# for.clim.N2 <- forecast_new_loc(obs = N2, pars = mod3$pars, covars = enviro.var.N2, Dsq = Dsq2, mod = 'clim')
# for.topoclim.N2 <- forecast_new_loc(obs = N2, pars = mod4$pars, covars = enviro.var.N2, Dsq = Dsq2, mod = 'topoclim')
# 
# save.image(file = "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_1_near_predictions/35y_OOS_1_near_predictions_5y_average_initial.RData")

# N3
for.base.N3 <- forecast_new_loc(obs = N3, pars = mod1$pars, Dsq = Dsq3, mod = 'base')
for.topo.N3 <- forecast_new_loc(obs = N3, pars = mod2$pars, covars = enviro.var.N3, Dsq = Dsq3, mod = 'topo')
for.clim.N3 <- forecast_new_loc(obs = N3, pars = mod3$pars, covars = enviro.var.N3, Dsq = Dsq3, mod = 'clim')
for.topoclim.N3 <- forecast_new_loc(obs = N3, pars = mod4$pars, covars = enviro.var.N3, Dsq = Dsq3, mod = 'topoclim')

save.image(file = "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/35y_OOS_2_far_predictions_5y_average_initial_v2.RData")
