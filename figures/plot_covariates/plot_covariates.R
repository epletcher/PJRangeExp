# load packages
library(abind)
library(tidyverse)

# ---covariate data----

## ----- N -----

# clim data (read in and standardize)
climdat <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/PrismClimateDataV3.csv") 
#%>% plyr::mutate(across(c("ppt", "tmean", "vpdmax"), scale))

# topo data (read in and standardize)
topodat <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/topographic_data.csv") 
#%>% dplyr::mutate(across(c("heatload", "elev"), scale))

# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- climdat %>% dplyr::select(c(cellnum, year, ppt)) %>% 
  pivot_wider(names_from = cellnum, values_from = ppt) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

tmean <- climdat %>% dplyr::select(c(cellnum, year, tmean)) %>% 
  pivot_wider(names_from = cellnum, values_from = tmean) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

vpdmax <- climdat %>% dplyr::select(c(cellnum, year, vpdmax)) %>% 
  pivot_wider(names_from = cellnum, values_from = vpdmax) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

heatload <- topodat %>% dplyr::select(c(cellnum, year, heatload)) %>% 
  pivot_wider(names_from = cellnum, values_from = heatload) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

elev <- topodat %>% dplyr::select(c(cellnum, year, elev)) %>% 
  pivot_wider(names_from = cellnum, values_from = elev) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

#stack
enviro.var <- abind(ppt, tmean, vpdmax, heatload, elev, along = 3)

# in-sample un-standardized covariates
N.clim <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/PrismClimateDataV3.csv")
N.topo <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/topographic_data.csv")

## ----- N2 -----
N2.clim <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/PrismClimateDataV3_oos1.csv") #%>% 
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%
  # 
  # dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%
  # 
  # dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N2.topo <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/topographic_data_oos1.csv") #%>% 
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(heatload = (heatload-mean(N.topo$heatload))/sd(N.topo$heatload)) %>%
  # 
  # dplyr::mutate(elev = (elev-mean(N.topo$elev))/sd(N.topo$elev)) 

# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N2.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N2.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N2.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N2.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N2 <- abind(ppt, tmean, heatload, elev, along = 3)

# ## N3
N3.clim <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/PrismClimateDataV3_oos2.csv") #%>%
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%
  # 
  # dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%
  # 
  # dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N3.topo <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/topographic_data_oos2.csv") #%>%
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(heatload = (heatload-mean(N.topo$heatload))/sd(N.topo$heatload)) %>%
  # 
  # dplyr::mutate(elev = (elev-mean(N.topo$elev))/sd(N.topo$elev))

# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N3.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N3.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N3.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N3.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N3 <- abind(ppt, tmean, heatload, elev, along = 3)

# median ppt
ppt_N <- apply(enviro.var[,,1], MARGIN = 1, FUN = mean)
ppt_N2 <- apply(enviro.var.N2[,,1], MARGIN = 1, FUN = mean)
ppt_N3 <- apply(enviro.var.N3[,,1], MARGIN = 1, FUN = mean)

plot(ppt_N, type = 'l', col = 'red', ylim = c(50,550), lwd = 1.5, ylab = "Annual Precipitation (mm)")
lines(ppt_N2, col = '#00B0F0', lwd = 1.5)
lines(ppt_N3, col = 'purple', lwd = 1.5)

# mean tmean
tmean_N <- apply(enviro.var[,,2], MARGIN = 1, FUN = mean)
tmean_N2 <- apply(enviro.var.N2[,,2], MARGIN = 1, FUN = mean)
tmean_N3 <- apply(enviro.var.N3[,,2], MARGIN = 1, FUN = mean)

plot(tmean_N, type = 'l', col = 'red', ylim = c(4,9), lwd = 1.6, ylab = "Mean Temp (c)")
lines(tmean_N2, col = '#00B0F0', lwd = 1.6)
lines(tmean_N3, col = 'purple', lwd = 1.6)

# mean heatload
mean(enviro.var[,,4]) # n
mean(enviro.var.N2[,,3]) # n2
mean(enviro.var.N3[,,3]) # n3

# mean elev
mean(enviro.var[,,5]) # n
sd(enviro.var[,,5])
mean(enviro.var.N2[,,4]) # n2
sd(enviro.var.N2[,,4])
mean(enviro.var.N3[,,4]) # n3
sd(enviro.var.N3[,,4])

# ** test plotting variation in environmental covariates for each study area **
par(mfrow = c(2,2))
matplot(enviro.var[,,1], type = 'l') # ppt
matplot(enviro.var[,,2], type = 'l') # tmean
matplot(enviro.var[,,4], type = 'l') # heatload
matplot(enviro.var[,,5], type = 'l') # elev

matplot(enviro.var.N2[,,1], type = 'l') # ppt
matplot(enviro.var.N2[,,2], type = 'l') # tmean
matplot(enviro.var.N2[,,3], type = 'l') # heatload
matplot(enviro.var.N2[,,4], type = 'l') # elev

matplot(enviro.var.N3[,,1], type = 'l') # ppt
matplot(enviro.var.N3[,,2], type = 'l') # tmean
matplot(enviro.var.N3[,,3], type = 'l') # heatload
matplot(enviro.var.N3[,,4], type = 'l') # elev

# median precipitation for each area
med.precip
