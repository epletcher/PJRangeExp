# evaluate model ft for in sample forecasts
# In sample evaluation evaluates forecast over 5 withheld years

# ** script takes a while to run, load 'eval_model_in_sample_5y.RData' for model predictions and rmse's **
# In R:/Shriver_Lab/PJspread/evaluate_in_sample

# load packages
library(tidyverse)
library(abind) # arrays
library(spatstat) # distance matrix w/ pairwise dist. function
library(raster)
library(cowplot)

# --------- functions -----------
source("forecastFunctions.R")

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

## set working dir
setwd("R:/Shriver_Lab/PJspread/sampleroutput")

# ---- Load model outputs from file ----

# Pull parameter and latent state estimates based on model version, base, topo, clim or topoclim
retrieve_mod <- function(model) {
  if(model=="topo") {model="topo_"}
  if(model=="base") {model="base_v"}
  if(model=="clim") {model="r_clim"}
  if(model %in% c("base_v","topo_","topoClim","r_clim","quad")) {
    
    files <- list.files("R:/Shriver_Lab/PJspread/sampleroutput", pattern=model, include.dirs = TRUE)
      
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

mod1 <- retrieve_mod("base") # model type here
mod2 <- retrieve_mod("topo")
mod3 <- retrieve_mod("clim")
mod4 <- retrieve_mod("topoClim")

# # save model objects as RDS for out of sample forecasts
# saveRDS(mod1, file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/base_model.rds")
# saveRDS(mod2, file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/topo_model.rds")
# saveRDS(mod3, file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/clim_model.rds")
# saveRDS(mod4, file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/topoclim_model.rds")

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
  
  no.z.Nlat <- replace(mod$Nlat, mod$Nlat<0, 0) # replace negative latent values with zero
  
  # calculate predicted values
  for(i in 1:length(pars$tauOut)) { 
    
    M <- dispersal(Dsq, pars$tauOut[i])
    
    for(p in 1:pixels){
      M1[,p]=M[,p]/sum(M[,p])
    }
    
    Nt <- no.z.Nlat[31,,i] # set initial cover value as actual latent value for the last year of model fit
    
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
      
      # if(mod=='base_q') {
      #   G <- growth(b0=pars$beta0[i],b1=pars$beta1[i],b2=pars$beta2,as.vector(Nt))
      # }
      # 
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

save.image(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/eval_model_insample_5y_rm_z.RData")
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
  mutate(rmse = rmse/max(N)) %>% # normalize rmse values by range of observed data
  ggplot(aes(x = model, y = rmse, col = model)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  scale_color_manual(values=group.cols) +
  labs(x = "MODEL", y = "NORMALIZED RMSE") + 
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=20))

# save to server
ggsave("R:/Shriver_Lab/PJspread/figures/rmse-plot_insample_5year_forecast_normalized.png", plot = last_plot(), dpi = 400)

# save to Google drive project folder
ggsave("G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Figures/Final_model_run/rmse-plot_insample_5year_forecast_normalized.png", plot = last_plot(), dpi = 400)

## ---- Calculate median and 95% Credible interval ----

med.pred.base <- apply(for.base.N$predOut, MARGIN = c(1,2), FUN = median)
write.csv(med.pred.base, "R:/Shriver_Lab/PJspread/evaluate_in_sample/insample_predictions/medpredpixel_base.csv")
up.pred.base <- apply(for.base.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.95)
low.pred.base <- apply(for.base.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.05)

med.pred.topo <- apply(for.topo.N$predOut, MARGIN = c(1,2), FUN = median)
up.pred.topo <- apply(for.topo.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.95)
low.pred.topo <- apply(for.topo.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.05)
write.csv(med.pred.topo, "R:/Shriver_Lab/PJspread/evaluate_in_sample/insample_predictions/medpredpixel_topo.csv")

med.pred.clim <- apply(for.clim.N$predOut, MARGIN = c(1,2), FUN = median)
low.pred.clim <- apply(for.clim.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.05)
up.pred.clim <- apply(for.clim.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.95)
write.csv(med.pred.clim, "R:/Shriver_Lab/PJspread/evaluate_in_sample/insample_predictions/medpredpixel_clim.csv")

med.pred.topoclim <- apply(for.topoclim.N$predOut, MARGIN = c(1,2), FUN = median)
write.csv(med.pred.topoclim, "R:/Shriver_Lab/PJspread/evaluate_in_sample/insample_predictions/medpredpixel_topoclim.csv")
up.pred.topoclim <- apply(for.topoclim.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.95) 
low.pred.topoclim <- apply(for.topoclim.N$predOut, MARGIN = c(1,2), FUN = quantile, 0.05)

# ------------------ Predicted vs. Observed plots -----------------------

# ** these plots had an error and are not currently used to generate any results **

# base
#**** weird thing happening here, obs[2:5,]-ob[1,] i not equal to change in cover at 1,2,3,4 years out! seems like the matrix vector subtraction wasn't working the way we thought, but this seems to work:
# sweep function subtracts vector from matrix
sweep(obs[2:5,], 2, obs[1,])

## use sweep to calc pred vs. obs
plot(sweep(obs[2:5,], 2, obs[1,]),sweep(med.pred.base[2:5,], 2, med.pred.base[1,]), main = "BASE MODEL", xlab = "OBSERVED CHANGE IN COVER (%)", ylab= "PREDICTED CHANGES IN COVER (%)", col = "aquamarine3", pch = 20, cex.lab = 1.4)
abline(0,1, col = "purple", lwd = 2.5)

# # old way from Bob:
# plot(obs[5,]-obs[1,],med.pred[5,]-med.pred[1,], main = "BASE MODEL", xlab = "OBSERVED CHANGE IN COVER (%)", ylab= "PREDICTED CHANGES IN COVER (%)", col = "aquamarine3", pch = 20, cex.lab = 1.4)
# abline(0,1, col = "purple", lwd = 2.5)
# # TOPO
# plot(obs[2:5,]-obs[1,],med.pred.topo[2:5,]-med.pred.topo[1,], main = "Topographic (spatially varying)", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3", pch = 20)
# abline(0,1, col = "purple", lwd = 2)
# # clim
# plot(obs[2:5,]-obs[1,],med.pred.clim[2:5,]-med.pred.clim[1,], main = "Climatic (temporally varying)", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3", pch = 20)
# abline(0,1, col = "purple", lwd = 2)
# # base
# plot(obs[2:5,]-obs[1,],med.pred.topoclim[2:5,]-med.pred.topoclim[1,], main = "Topoclimatic (spatial and temporally varying)", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3", pch = 20)
# abline(0,1, col = "purple", lwd = 2)

# ---- plot model predictions over time ----

# ** plots of predictions through time we use the predictions forecasted across entire 36 year time period, not just the last five years **

#par(mfrow = c(2,2))

plot_pix <- function(pix,obs) {
  # color palette
  cryo.pal <- c("#C96B95","#5bb2dd","#7ABD61","#E2AC3B","#747BA0","#EC5D4E")
  
  plot(med.pred.base[,pix], type = "l", ylim = c(0,25), ylab = "% cover", xlab = "years out", col = "coral", lwd = 2)
  lines(up.pred.base[,pix],lty = 2,col = "coral")
  lines(low.pred.base[,pix],lty = 2,col = "coral")
  # 
  # lines(med.pred.topo[,pix],col = cryo.pal[4],lwd = 2)
  # # lines(low.topo[,pix],lty = 2,col = cryo.pal[4])
  # # lines(up.topo[,pix],lty = 2,col = cryo.pal[4])
  # # 
  # lines(med.pred.clim[,pix], col = cryo.pal[5],lwd = 2)
  # # lines(low.clim[,pix],lty = 2, col = cryo.pal[5])
  # # lines(up.clim[,pix],lty = 2,col = cryo.pal[5])
  
  lines(med.pred.topoclim[,pix],col = "aquamarine3",lwd=2)
  lines(low.pred.topoclim[,pix], lty = 2, col = "aquamarine3")
  lines(up.pred.topoclim[,pix], lty = 2, col = "aquamarine3")

  points(obs[,pix], pch = 20)
  lines(obs[,pix], lty = 3)
}

plot_pix(68,N)

## ggplot version
plot_pix_gg <- function(pix) {
  
  px <- as.character(pix)
  pxv <- paste("V",as.character(pix), sep = "")

#base
med.base <- med.pred.base %>% as.data.frame() 
med.base$year <- row.names(med.base) 
med <- med.base %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "med.pred.base", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)

low.base <- low.pred.base %>% as.data.frame() # low
low.base$year <- row.names(low.base)
low <- low.base %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "low.pred.base", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)

up.base <- up.pred.base %>% as.data.frame() # high
up.base$year <- row.names(up.base)
up <- up.base %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "up.pred.base", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)

# #topoclim
# med.topoclim <- med.pred.topoclim %>% as.data.frame() 
# med.topoclim$year <- row.names(med.topoclim) 
# med.tc <- med.topoclim %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "med.pred.topoclim", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)
# 
# low.topoclim <- low.pred.topoclim %>% as.data.frame() # low
# low.topoclim$year <- row.names(low.topoclim)
# low.tc <- low.topoclim %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "low.pred.topoclim", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)
# 
# up.topoclim <- up.pred.topoclim %>% as.data.frame()  # high
# up.topoclim$year <- row.names(up.topoclim)
# up.tc <- up.topoclim %>% mutate(year = as.numeric(year)+31) %>% pivot_longer(!year, values_to = "up.pred.topoclim", names_to = "pixel") %>% filter(pixel == pxv) %>% dplyr::select(-pixel)

# observed 
Ndat <- N %>% as.data.frame() %>% rownames_to_column("year") %>% pivot_longer(!year, values_to = "obs", names_to = "pixel") %>% filter(pixel == px) %>% mutate(year = as.numeric(year)) %>% dplyr::select(-pixel)

# combine predictions and observations into one dataframe
plot.dat <- left_join(Ndat, med) %>% left_join(.,low) %>% left_join(.,up) 
#%>% left_join(.,up.tc) %>% left_join(.,low.tc) %>% left_join(.,med.tc)

# plot
(test <- plot.dat %>% ggplot(aes(x = year, y = obs)) + 
    geom_point(aes(x = year, y = obs)) + 
    geom_line(aes(x = year, y = obs), lty = 2) + 
    scale_y_continuous(limits = c(0, 25))) + 
  geom_ribbon(aes(ymin = low.pred.base, ymax = up.pred.base), alpha=0.4, fill = "#F8766D") + 
  #geom_ribbon(aes(ymin = low.pred.topoclim, ymax = up.pred.topoclim), alpha=0.3, fill = "#C77Cff") + 
  geom_line(aes(x = year, y = med.pred.base), lwd = 1.5, col = "#F8766D") + 
  #geom_line(aes(x = year, y = med.pred.topoclim), lwd = 1.5, col = "#C77Cff") +      theme_bw() + 
  labs(x = "YEARS OUT", y = "TREE COVER (%)") + theme(legend.position="none", text = element_text(size=25))
}

plot_pix_gg(43)

# test out different pixels to plot based on 36-yr trends
N %>% as.data.frame() %>% rownames_to_column("year") %>% pivot_longer(!year, values_to = "obs", names_to = "pixel") %>% filter(pixel == "94") %>% mutate(year = as.numeric(year)) %>% ggplot(aes(x = year, y = obs)) + geom_point(aes(x = year, y = obs)) + geom_line(aes(x = year, y = obs), lty = 2) + scale_y_continuous(limits = c(0, 25)) 

# plot focal pixel
tree.raster <- raster::raster("R:/Shriver_Lab/PJspread/data_prepping/Data/climate_and_tree_cover_data/RAP_treeCover_1986.tif") %>% raster::aggregate(fact = 30)

raster.dat <- raster2df(tree.raster)
raster.dat[94,4] <- -150 # re-assign focal pixel value to map

focR <- raster.dat %>% as.data.frame() %>% dplyr::select(-cellnum) %>% df2raster()

breakr <- c(-150,-1,5,10,15,20,25)
col5 = c("hotpink", brewer.pal(7,"YlGn"))
plot(focR, breaks = breakr, col = col5, legend = F)