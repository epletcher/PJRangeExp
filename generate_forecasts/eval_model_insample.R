# evaluate model ft for in sample forecasts
# In sample evaluation evaluates forecast over 5 withheld years
# Produces figure 3 from main text

# load packages
library(tidyverse)
library(abind)

# --------- functions -----------
source("forecastFunctions.R")

## set working dir. to sampler outputs folder
setwd("FILEPATH/sampleroutput")

# ---- Load model outputs from file ----

## Pull parameter and latent state estimates based on model version, base, topo, clim or topoclim
retrieve_mod <- function(model) {
  if(model=="topo") {model="topo_"}
  if(model=="base") {model="base_v"}
  if(model=="clim") {model="r_clim"}
  if(model %in% c("base_v","topo_","topoClim","r_clim","quad")) {
    
    files <- list.files("FILEPATH/sampleroutput", pattern=model, include.dirs = TRUE)
      
    # chain 1
      load(files[1])
      
      burninone <- burnin+1
      
      alphaOut <- alphaOut[burninone:Niter,]
      betaOut <- betaOut[burninone:Niter,]
      tauOut <- tauOut[burninone:Niter,]
      sig.pOut <- sig.pOut[burninone:Niter,]
      
      if(is.matrix(alphaOut)) {
        a_names <- dim(alphaOut)[2] 
        for(i in 1:dim(alphaOut)[2]) {
          a_names[i] <- paste("alpha",as.character(i-1),sep='')
        }
        colnames(alphaOut) <- a_names
      }
      
      if(is.matrix(betaOut)) {
        b_names <- dim(betaOut)[2]
        for(i in 1:dim(betaOut)[2]) {
          b_names[i] <- paste("beta",as.character(i-1),sep='')
        }
        colnames(betaOut) <- b_names
      }
      
      mod <- cbind(alphaOut,betaOut,tauOut,sig.pOut)
      
      #rmseTotOut_comb <- rmseTotOut
      latburn <- length(tenIter)/2+1
      
      Nlat_comb <- NlatOut[,,latburn:length(tenIter)]
      
      # chain 2 +
      for(c in 2:length(files)) {
        
        load(files[c])
        
        alphaOut <- alphaOut[burninone:Niter,]
        betaOut <- betaOut[burninone:Niter,]
        tauOut <- tauOut[burninone:Niter,]
        sig.pOut <- sig.pOut[burninone:Niter,]
        
        # combine all chains
        mod <- rbind(mod, cbind(alphaOut,betaOut,tauOut,sig.pOut)) 
        
        #rmseTotOut_comb <- cbind(rmseTotOut_comb, rmseTotOut) 
        Nlat_comb <- abind(Nlat_comb, NlatOut[,,latburn:length(tenIter)], along = 3) 
      }
    
  
  # thinning outputs here
  to_return <- list()
  thin.iter <- seq(10,dim(mod)[1], by = 10) # vector of every 10th iteration
  to_return$pars <- mod[thin.iter,]
  #to_return$rmse <- rmseTotOut_comb[,thin.iter]
  to_return$Nlat <- Nlat_comb # an array of estimated latent state for every pixel-year-iteration
  }
  
  return(to_return)
}

## retrieve model estimats for all four models:
mod1 <- retrieve_mod("base") # model type here
mod2 <- retrieve_mod("topo")
mod3 <- retrieve_mod("clim")
mod4 <- retrieve_mod("topoClim")

# ----------------- Forecast across withheld years -------------------
# fixed/single origin forecast
# allow error to propagate through time, starting at year 31, through 36

forecast_is_loc <- function(obs, mod, covars, Dsq, mod.name) {
  
  pars <- as.data.frame(mod$pars)
  
  pixels = length(obs[1,])
  years = length(obs[,1])
  
  # empty array to fill
  predOut <- array(NA, c(years, pixels, length(pars$tauOut))) 
  
  # evaluation
  rmseTotOut<-biasOut<-denseOut<-rep(0,length(pars$tauOut)) # evaluation metrics
  
  # empty matrix for dispersal matrix
  M1 <- matrix(NA, pixels, pixels)
  
  # calculate predicted values
  for(i in 1:length(pars$tauOut)) { 
    
    M <- dispersal(Dsq, pars$tauOut[i])
    
    for(p in 1:pixels){
      M1[,p]=M[,p]/sum(M[,p])
    }
    
    Nt <- mod$Nlat[31,,i] # set initial cover value as actual latent value for the last year of model fit
    
    for(t in 1:years) {
      
      if(mod.name=='base') {
        G <- growthsimple(a=pars$alphaOut[i],b=pars$betaOut[i],nt=Nt)
      }
      
      if(mod.name=='clim') {
        G <- growthClim(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],b0=pars$betaOut[i],X=covars[t,,],nt=Nt)
      }
      
      if(mod.name=='topo') {
        G <- growthTopo(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2[i],X=covars[t,,],nt=Nt)
      }
      
      if(mod.name=='topoclim') {
        G <- growthTopoClim(a0=pars$alpha0[i],a1=pars$alpha1[i],a2=pars$alpha2[i],a3=pars$alpha3[i],a4=pars$alpha4[i],b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2[i],X=covars[t,,],nt=Nt)
      }
      
      Nmean <-M1%*%(diag(G)%*%Nt)
      
      Nt <- rnorm(pixels, Nmean, pars$sig.pOut[i])
      
      Nt <- replace(Nt, Nt<0, 0) # convert negative values to zeroes
      
      predOut[t,,i] <- Nt
      
    }
    
    rmseTotOut[i] <- rmsefunc(pred=predOut[,,i], obs=obs)
    # biasOut[t,,i] <- biasfunc(pred=predOut[t,,i], obs=obs[t,])
    # denseOut[t,,i] <- densefunc(pred=predOut[t,,i], obs=obs[t,], sig_o=1)
    
    print(i)  
  }
  to_return <- list()
  to_return$predOut <- predOut
  to_return$rmseTotOut <- rmseTotOut
  # to_return$biasOut <- biasOut
  # to_return$denseOut <- denseOut
  
  return(to_return)
}

# ----------- RUN FORECASTS ----------------------------
obs <- N[32:36,] # starting at year 31, through 36

for.base.N <- forecast_is_loc(obs = obs, mod = mod1, Dsq = Dsq, mod.name = 'base')
for.topo.N <- forecast_is_loc(obs = obs, mod = mod2, covars = enviro.var[,,-3], Dsq = Dsq, mod.name = 'topo')
for.clim.N <- forecast_is_loc(obs = obs, mod = mod3, covars = enviro.var[,,-3], Dsq = Dsq, mod.name = 'clim')
for.topoclim.N <- forecast_is_loc(obs = obs, mod = mod4, covars = enviro.var[,,-3], Dsq = Dsq, mod.name = 'topoclim')

# save forecast to file 
save.image(file = "FILEPATH/eval_model_insample_5y.RData")

# ---------- SUMMARIZE AND PLOT RMSE ------------------
# average rmse functions
average_rmse <- function(mod) {
  rmse <- mod$rmseTotOut

  to_return <- data.frame(
  avg = median(rmse),
  low.cred = quantile(rmse, 0.1),
  up.cred = quantile(rmse, 0.9))
  return(to_return)
}

average_rmse(for.base.N)
average_rmse(for.topo.N)
average_rmse(for.clim.N)
average_rmse(for.topoclim.N)

normalized_rmse <- function(mod) {
  rmse <- mod$rmseTotOut
  nrmse <- rmse/(max(N)-min(N))
  
  to_return <- data.frame(
    avg = median(nrmse),
    low.cred = quantile(nrmse, 0.1),
    up.cred = quantile(nrmse, 0.9))
  return(to_return)
}

normalized_rmse(for.base.N)
normalized_rmse(for.topo.N)
normalized_rmse(for.clim.N)
normalized_rmse(for.topoclim.N)


## plot RMSE
BASE = for.base.N$rmseTotOut
TOPO = for.topo.N$rmseTotOut
CLIM = for.clim.N$rmseTotOut
TOPOCLIM = for.topoclim.N$rmseTotOut

CLIM = c(CLIM, rep(NA, length(BASE) - length(CLIM)))
TOPOCLIM = c(TOPOCLIM, rep(NA, length(BASE) - length(TOPOCLIM)))

rmsedat <- data.frame(BASE, TOPO, CLIM,TOPOCLIM)

# plot
# colors
group.cols <- c("#F8766D","#00BFC4","#C77Cff","#7CAE00")

rmsedat %>% pivot_longer(everything(), names_to = "model", values_to = "rmse") %>% 
  mutate(rmse = rmse/max(N)) %>% # normalize rmse values by range of observed data (rmse/max(N)) or leave un normalized
  ggplot(aes(x = model, y = rmse, col = model)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  scale_color_manual(values=group.cols) +
  labs(x = "MODEL", y = "NRMSE") + 
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=16))