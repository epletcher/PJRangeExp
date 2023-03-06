#.libPaths("C:/Rpackages/R/win-library/4.1")

# load packages
library(maptools)
library(raster)
library(rgdal)
library(spatstat)
library(tidyverse)
library(rstan)
library(abind)
library(shinystan)
library(msm) # truncated normal 
library(OneR) # binning function

# load model workspace
load("G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Scripts/Sampler/sampler_output_worskpace_files/sampler_base_v3.RData")

burnin = .5*Niter

#------------------- Inspect model convergence ----------------------

# trace plot
par(mfrow = c(2, 3))
#tau
plot(tauOut[burnin:Niter], type = "l", xlab = "Iteration",
     col = "#a6bddb", ylab = expression(tau))
abline(h = mean(tauOut[burnin:Niter]), col = "#1c9099")
# beta 1
plot(betaOut[burnin:Niter,1], type = "l", xlab = "Iteration",
     col = "#a6bddb", ylab = expression(beta[0]))
abline(h = mean(betaOut[burnin:Niter,1]), col = "#1c9099")
#beta 2
plot(betaOut[burnin:Niter,2], type = "l", xlab = "Iteration",
     col = "#a6bddb", ylab = expression(beta[1]))
abline(h = mean(betaOut[burnin:Niter,2]), col = "#1c9099")
# sigma obs
plot(sig.oOut[burnin:Niter], type = "l", xlab = "Iteration",
     col = "#a6bddb", ylab = expression(sigma[o]))
abline(h = mean(sig.oOut[burnin:Niter]), col="#1c9099")
# sigma process
plot(sig.pOut[burnin:Niter], type = "l", xlab = "Iteration",
     col = "#a6bddb", ylab = expression(sigma[p]))
abline(h = mean(sig.pOut[burnin:Niter]), col="#1c9099")

# density plots
par(mfrow = c(2, 3))
#tau
hist(tauOut[burnin:Niter], breaks = 40, xlab = "", main = expression(tau))
abline(v = mean(tauOut[burnin:Niter]), col = "#1c9099", lwd = 2)
# beta 1
hist(betaOut[burnin:Niter,1], breaks = 40, xlab = "", main = expression(beta[1]))
abline(v = mean(betaOut[burnin:Niter,1]), col = "#1c9099", lwd = 2)
# beta 2
hist(betaOut[burnin:Niter,2], breaks = 40, xlab = "", main = expression(beta[2]))
abline(v = mean(betaOut[burnin:Niter,2]), col = "#1c9099", lwd = 2)
# sigma obs
hist(sig.oOut[burnin:Niter], breaks = 40, xlab = "", main = expression(sigma[o]))
abline(v = mean(sig.oOut[burnin:Niter]), col = "#1c9099", lwd = 2)
# sigma process
hist(sig.pOut[burnin:Niter], breaks = 40, xlab = "", main = expression(sigma[p]))
abline(v = mean(sig.pOut[burnin:Niter]), col = "#1c9099", lwd = 2)

# ------------ Inspect latent states ------------------

## all latent cover values for each iteration plotted
rep.pix <- c(115:145, 910:940, 1865:1895) # representative pixels (high,med,low density)
par(mfrow = c(2,2))
# high cover pixels
matplot(NlatOut[,115:145,Niter/10], type = "l")
#medium cover pixels
matplot(NlatOut[,910:940,Niter/10], type = "l")
#low cover pixels
matplot(NlatOut[,1865:1895,Niter/10], type = "l")

## mean of param estimates
mean.lat.cover <- apply(NlatOut[,,1000:Niter/10], MARGIN = c(1,2), FUN = mean)

# credible intervals for latent states
lower.lat <- apply(NlatOut[,,burnin/10:Niter/10], MARGIN = c(1,2), FUN = quantile, 0.975)
upper.lat <- apply(NlatOut[,,burnin/10:Niter/10], MARGIN = c(1,2), FUN = quantile, 0.025)

# change in latent cover
mean.lat.change <- mean.lat.cover[2:31,] - mean.lat.cover[1:30,]

# # relative change in latent cover
mean.rel.lat.change <- log(mean.lat.cover[2:31,]) - log(mean.lat.cover[1:30,])

## all mean latent cover values plotted over time
# grouped by pixel type
par(mfrow = c(2,2))
# high cover pixels
matplot(mean.lat.cover[,115:145], type = "l")
#medium cover pixels
matplot(mean.lat.cover[,910:940], type = "l")
#low cover pixels
matplot(mean.lat.cover[,1865:1895], type = "l")


# ------------ Parameter estimates --------------------
# # pull out parameter estimates with specified burnin
# pars <- output %>% dplyr::select(c("tau", "beta0", "a.1", "a.2", "a.3", "beta1", "sigma")) # "a.1", "a.2", "a.3",
# pars <- pars[burnin:niter,]

# combine all pars after burnin into dataframe
# remove burnin and name
pars <- as.data.frame(cbind(tauOut[burnin:Niter],betaOut[burnin:Niter,1],beta1 = betaOut[burnin:Niter,2],sig.p = sig.pOut[burnin:Niter]))
colnames(pars) <- c("tau", "beta0", "beta1","sig.p")

# Thin parameter estimates to be the same as the subset of iterations in the Nlatent output (NlatOut)
latburnin <- tenIter[1001:2000] # vector of every 10th iteration, minus burnin

latPars <- as.data.frame(cbind(tauOut[latburnin],betaOut[latburnin,1],beta1 = betaOut[latburnin,2],sig.p = sig.pOut[latburnin]))
colnames(latPars) <- c("tau", "beta0", "beta1","sig.p")

# check mean parameter estimates
pars %>% as.data.frame() %>% summarize_all(quantile, 0.025)
pars %>% as.data.frame() %>% summarize_all(mean)
pars %>% as.data.frame() %>% summarize_all(quantile, 0.975)

# ------------ Calculate predictions using parameter estimates

# reload observed data (so that all years are present)
N <- treedat %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

#------------ Functions ---------------
## dispersal
# D = matrix of distances between grid cells (should already be squared)
# tau = dispersal constant
dispersal = function(D, tau) { # D should be distance matrix already squared
  M = exp(-(D/tau))# output is dispersal matrix
  return(M)
}


## growth
# # parameterization #1
# a + b*T is regression line for growth rate with no density dependence
# k is carrying capacity
# nt is vector of data, autoregressive term

# parameterization #2
growthsimple = function(a, b, nt) {
  G = exp(a+b*nt) # output is growth matrix
  return(G)
}

# ----------------------- evaluate 1-year predictability ------------------------
# predicted cover based off estimated latent state
# this will demonstrate how well the process model is working



pixels = 2352
years = 36
ylat = 31 # number of years for latent state
latIter = 1000 # number of latent state iterations, not including burning 2000-1000=1000 

## calculate predicted cover
# empty matrix for predicted cover at one time step out
pred.cover <- array(NA, c(ylat, pixels, latIter))
pred.start <- abind(list(replicate(latIter, t(as.matrix(N[1,]))))) # fill in the actual data for 1st year
pred.cover[1,,] <- pred.start[1,,]

NlatBurn <- NlatOut[,,1001:2000]
# empty matrix for dispersal matrix
M1 <- matrix(NA, pixels, pixels)

for(i in 1:latIter) {
  
  M <- dispersal(Dsq, latPars$tau[i])
  
  for(p in 1:pixels){
    M1[,p]=M[,p]/sum(M[,p])
  }
  for(t in 2:ylat) {
    
    G <- growthsimple(latPars$beta0[i],latPars$beta1[i], as.vector(NlatBurn[t-1,,i]))
    
    H <- G * M1
    
    pred.cover[t,,i] <- H %*% NlatBurn[t-1,,i]
    
  }
  
}

## Calculate predictions using posterior means
# mean predicted cover per pixel
mean.pred.cover <- apply(pred.cover, MARGIN = c(1,2), FUN = mean)

# mean predicted change in cover per pixel
mean.pred.change <- mean.pred.cover[2:ylat,] - N[1:ylat-1,]

# # mean predicted relative change
mean.pred.rel.change <- log(mean.pred.cover[2:ylat,]) - log(N[1:ylat-1,])

# ----------------- Change in cover -----------------------
# 
# ## calculate observed change in cover
# 
# obs.change <- matrix(NA, ylat, numpix)
# obs.change[1,] <- 0
# 
# for(t in 2:ylat)  {
#   obs.change[t,] <- N[t,rep.pix]-N[t-1,rep.pix]
# }
# 
# 
# ## calculate log transformed observed change in *relative* cover
# # note: just let zeroes turn into negative inifinity and get filtered out for now
# # Nn <- N
# # Nn[Nn==0] <- min(Nn[Nn>0]) # find the minimum observed value and reassign zero vals
# 
# obs.rel.change <- matrix(NA, ylat, numpix)
# obs.rel.change[1,] <- 0
# 
# for(t in 2:ylat)  {
#   obs.rel.change[t,] <- log(N[t,rep.pix])-log(N[t-1,rep.pix])
# }
# 

# ----------------------Model coverage-------------------------
# # Test for coverage - how much of the true variability in the dataset is captured by the model's estimate of sigma?
# outofbounds <- 0
# for(t in 2:ylat) {
#   for(i in 1:numpix) {
#     error <- mean.pred.cover[t,i] - mean.lat.cover[t,i]
#     if(error < qtnorm(0.025, 0, mean(pars$sig.p))) {outofbounds <- outofbounds+1}
#     if(error > qtnorm(0.975, 0, mean(pars$sig.p))) {outofbounds <- outofbounds+1}
#   }
# }
# print((numpix*ylat-outofbounds)/(numpix*ylat))
# 

# -------------- plot pred vs. obs (ggplot version) ---------------
# change in cover
pred.dat <- t(mean.pred.change) %>% as.data.frame()
colnames(pred.dat) <- (as.character(1987:2016))
pred.dat.long <- pred.dat %>% pivot_longer(cols = everything(), names_to = "years2", values_to = "predChange") 

# # obs
# obs.dat <- t(obs.change) %>% as.data.frame() 
# colnames(obs.dat) <- (as.character(1986:2016))
# obs.dat.long <- obs.dat%>% dplyr::select(-"1986") %>% pivot_longer(cols = everything(), names_to = "years", values_to = "obsChange") %>% mutate(obs.bin = bin(obsChange, nbins=20, method = "content")) %>% mutate(bin.range = gsub("[(]|\\]","", as.character(obs.bin))) %>% separate(bin.range, c("min", "max"), sep = ",") %>% mutate(max = as.numeric(max)) %>% mutate(min = as.numeric(min)) %>% mutate(bin.range = max-min) %>% select(-c(min, max))

## lat
# total cover
lat.dat.cover <- t(mean.lat.cover) %>% as.data.frame()
colnames(lat.dat.cover) <- as.character(1986:2016)
cover.dat.long <- lat.dat.cover %>% select(-"1986") %>% pivot_longer(cols = everything(), names_to = "years", values_to = "latCover")
#change in cover
lat.dat <- t(mean.lat.change) %>% as.data.frame() 
colnames(lat.dat) <- (as.character(1987:2016))

lat.dat.long <- lat.dat %>% pivot_longer(cols = everything(), names_to = "years2", values_to = "latChange") %>% mutate(lat.bin = bin(latChange, nbins=20, method = "content")) %>% mutate(bin.range = gsub("[(]|\\]","", as.character(lat.bin))) %>% separate(bin.range, c("min", "max"), sep = ",") %>% mutate(max = as.numeric(max)) %>% mutate(min = as.numeric(min)) %>% mutate(bin.range = max-min) %>% select(-c(min, max)) %>% cbind(.,cover.dat.long) %>% select(-"years2")

# box plot
ggplot(lat.dat.long) + geom_point(aes(x = latChange, y = pred.dat.long$predChange, col = latCover), alpha = 0.3) + scale_colour_gradient(low = "#ffcc61", high = "#134e1d", name = "% latent cover") + geom_abline(intercept = 0, slope = 1, lwd = 1) + geom_boxplot(aes(x = latChange, y = pred.dat.long$predChange, group = lat.bin), width = lat.dat.long$bin.range, alpha = 0.5, outlier.shape = NA) + labs(title = "Assessing process model fit", x = "latent change in cover", y = "predicted change in cover") + theme_bw() 

#------------------------- plot Relative change in cover --------------------------

# lat
lat.rel.dat <- t(mean.rel.lat.change) %>% as.data.frame()
colnames(lat.rel.dat) <- (as.character(1987:2016))
lat.rel.dat.long <- lat.rel.dat %>% pivot_longer(cols = everything(), names_to = "years", values_to = "latChange")

# pred
pred.rel.dat <- t(mean.pred.rel.change) %>% as.data.frame() 
colnames(pred.rel.dat) <- (as.character(1987:2016))
pred.rel.dat.long <- pred.rel.dat %>% pivot_longer(cols = everything(), names_to = "years2", values_to = "predChange") %>% cbind(.,lat.rel.dat.long) %>% na.omit %>% select(-"years2") %>% mutate(lat.bin = bin(latChange, nbins=10, method = "content")) %>% mutate(bin.range = gsub("[(]|\\]","", as.character(lat.bin))) %>% separate(bin.range, c("min", "max"), sep = ",") %>% mutate(max = as.numeric(max)) %>% mutate(min = as.numeric(min)) %>% mutate(bin.range = max-min) %>% select(-c(min, max))

# # obs
# obs.dat <- t(obs.change) %>% as.data.frame() 
# colnames(obs.dat) <- (as.character(1986:2016))
# obs.dat.long <- obs.dat%>% dplyr::select(-"1986") %>% pivot_longer(cols = everything(), names_to = "years", values_to = "obsChange") %>% mutate(obs.bin = bin(obsChange, nbins=20, method = "content")) %>% mutate(bin.range = gsub("[(]|\\]","", as.character(obs.bin))) %>% separate(bin.range, c("min", "max"), sep = ",") %>% mutate(max = as.numeric(max)) %>% mutate(min = as.numeric(min)) %>% mutate(bin.range = max-min) %>% select(-c(min, max))

# plot with year grouping - only training years
ggplot(pred.rel.dat.long) + geom_point(aes(x = latChange, y = predChange, col = years), alpha = 0.4) + labs(title = "Assessing model fit", x = "latent change in cover", y = "predicted change in cover") + geom_abline(intercept = 0, slope = 1, lty = 2, lwd = 1.2) 

# box plot
ggplot(pred.rel.dat.long) + geom_point(aes(x = latChange, y = predChange), alpha = 0.1, col = "purple") + geom_abline(intercept = 0, slope = 1, lwd = 1) + geom_boxplot(aes(x = latChange, y = predChange, group = lat.bin), width = pred.rel.dat.long$bin.range, alpha = 0.5, outlier.shape = NA) + labs(title = "Assessing process model fit", x = "relative change in latent cover", y = "predicted change in relative cover") + theme_bw() 

# ---------------------- Calculate RMSE ------------------------
# # function to calculate RMSE
# rmsefunc <- function(pred, data) {
#   n <- length(pred)
#   RMSE <- sqrt((sum((pred-data)^2))/n)
#   return(RMSE)
# }
# 
# 
# # RMSE across years
# rmsedat <- rep(NA, length(mean.pred.change[2:36,512]))
# 
# for(i in 1:length(rmsedat)) {
#   pred <- mean.pred.change[2:36,512]
#   data <- obs.change[2:36,512]
#   rmsedat[i] <- rmsefunc(pred[i], data[i])
# }
# 
# plot(rmsedat, ylab = ("RMSE"), col = "red", type = "l")

# -------------- Plot latent cover and obs cover at each time step ---------
# rep.pix = representative pixels

# observed representative pixels
repN <- N

# one pixel obs and pred at each time step, for each category

par(mfrow = c(2,2))

# assign rep pixels for each class of pxel
hp=126
mp=975
lp=1880
par(mfrow = c(2,2))  
# high
plot(repN[,hp], col = "#d95f02", lwd = 2, ylab = "% cover", xlab = "year", ylim = c(5,30))
arrows(x0=c(1:36),y0=repN[,hp]-1, x1=c(1:36), y1=repN[,hp]+1,code=3,angle=0,length=0.05, col = "#d95f02")
lines(mean.lat.cover[,hp], col = "#1b9e77", lwd = 2)
lines(lower.lat[,hp], lty = 2, col = "#1b9e77", lwd = 2)
lines(upper.lat[,hp], lty = 2, col = "#1b9e77", lwd = 2)

# med
plot(repN[,mp], col = "#d95f02", lwd = 2, ylab = "% cover", xlab = "year", ylim = c(-3,15))
arrows(x0=c(1:36),y0=repN[,mp]-1, x1=c(1:36), y1=repN[,mp]+1,code=3,angle=0,length=0.05, col = "#d95f02")
lines(mean.lat.cover[,mp], col = "#1b9e77", lwd = 2)
lines(lower.lat[,mp], lty = 2, col = "#1b9e77", lwd = 2)
lines(upper.lat[,mp], lty = 2, col = "#1b9e77", lwd = 2)

# low
plot(repN[,lp], col = "#d95f02", lwd = 2, ylab = "% cover", xlab = "year", ylim = c(-3,7))
arrows(x0=c(1:36),y0=repN[,lp]-1, x1=c(1:36), y1=repN[,lp]+1,code=3,angle=0,length=0.05, col = "#d95f02")
lines(mean.lat.cover[,lp], col = "#1b9e77", lwd = 2)
lines(lower.lat[,lp], lty = 2, col = "#1b9e77", lwd = 2)
lines(upper.lat[,lp], lty = 2, col = "#1b9e77", lwd = 2)

  # ----------------- Forecast training years and withheld years -------------------
# fixed/single origin forecast
# allow error to propagate through time

# empty array to fill
# c(years, pixels, iterations)
pixels = 2352
years = 36
ylat = 31
numpix = length(rep.pix) # length of NLat values we saved for a subset of the pixels


predOut <- array(NA, c(years, pixels, length(latPars$tau))) # changing this to all 36 years, the first 31 will assess model fit, the last will be the forecast
predst <- abind(list(replicate(length(latPars$tau), t(as.matrix(N[1,]))))) # fill in the actual data for only the very first year???
predOut[1,,] <- predst[1,,]

# empty matrix for dispersal matrix
M1 <- matrix(NA, pixels, pixels)

# calculate predicted values
for(i in 1:length(latPars$beta0)) {
  
  M <- dispersal(Dsq, latPars$tau[i])
  
  for(p in 1:pixels){
    M1[,p]=M[,p]/sum(M[,p])
  }
  
  Nt <- NlatBurn[1,,i] # set initial cover value as actual latent value
  
  for(t in 2:years) {
    
    #G <- growthsimple(pars$beta0[i],pars$b[i], as.vector(Nt))
    G <- growthsimple(latPars$beta0[i],latPars$beta1[i], as.vector(Nt))
    
    H <- G * M1
    
    Nt <- rnorm(pixels, (H %*% Nt), latPars$sig.p[i]) # add process error, for each pixel
    #Nt <- H %*% Nt # calculates population at next time step without process variability
    predOut[t,,i] <- Nt
    
  }
  
}

# inspect output
predOut[, 200:208, 100]
predOut <-predOut[,,1:1000]

# # Calculate True mean change in cover for all data
# meantrue <- apply(N, MARGIN = 1, FUN = mean)
# 
# # plot predicted change in cover for withheld data
# meanpred <- apply(predOut, MARGIN = 1, FUN = mean)
# lowerpred <- apply(predOut, MARGIN = 1, FUN = quantile, 0.025)
# upperpred <- apply(predOut, MARGIN = 1, FUN = quantile, 0.975)
# 
# # plot
# plot(meanpred, type = "l", ylim = c(0, 25), col = "#C06C84", lwd = 2, ylab = "percent tree cover", xlab = "years out")
# lines(lowerpred, lty = 2, col = "#F8B195", lwd = 2,)
# lines(upperpred, lty = 2, col = "#F8B195", lwd = 2,)
# points(meantrue, col = "#355C7D", pch = 20)
# legend("topleft", legend = c("mean predicted cover","credible interval", "observations"), col = c("#C06C84", "#F8B195", "#355C7D"), lty = c(1, 2, 0), pch = c(NA, NA, 20), bty = "n")

# predOut mean by pixel
medpredpixel <- apply(predOut, MARGIN = c(1,2), FUN = median)
write.csv(medpredpixel, "RangeExp_model_predictions/medpredpixel_base.csv")

# calculate 95% credible intervals for each pixel
# low bound
lowboundpredpixel <- apply(predOut, MARGIN = c(1,2), FUN = quantile, 0.025)
write.csv(lowboundpredpixel, "RangeExp_model_predictions/lowboundpredpixel_base.csv")
# high bound
upboundpredpixel <- apply(predOut, MARGIN = c(1,2), FUN = quantile, 0.975)
write.csv(upboundpredpixel, "RangeExp_model_predictions/upboundpredpixel_base.csv")
