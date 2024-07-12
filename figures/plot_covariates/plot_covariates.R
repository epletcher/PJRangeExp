# load packages
library(abind)
library(tidyverse)
library(raster)
## set working directory to PJ_photo

# functions

# Function to re-rasterize dataframe
df2raster <- function(df) { # df should be x, y, vals only
  val <- names(df)[3]
  spg <- df %>% rename(values = val)
  sp::coordinates(spg) <- ~ x + y
  sp::gridded(spg) <- TRUE # coerce to SpatialPixelsDataFrame
  r <- raster::raster(spg, values = TRUE) # coerce to raster
  return(r)
}


# --- load covariate data----

# clim data (read in and standardize)
climdat <- read.csv("cover_spread/Data/PrismClimateDataV3.csv") 
#%>% plyr::mutate(across(c("ppt", "tmean", "vpdmax"), scale))

# topo data (read in and standardize)
topodat <- read.csv("cover_spread/Data/topography_data/topographic_data.csv") 
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

# # in-sample un-standardized covariates
# N.clim <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/PrismClimateDataV3.csv")
# N.topo <- read.csv("R:/Shriver_Lab/PJspread/data_prepping/topographic_data.csv")

## ----- N2 -----
N2.clim <- read.csv("cover_spread/Data/out_of_sample/PrismClimateDataV3_oos1.csv") #%>% 
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%
  # 
  # dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%
  # 
  # dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N2.topo <- read.csv("cover_spread/Data/out_of_sample/topographic_data_oos1.csv") #%>% 
  
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
N3.clim <- read.csv("cover_spread/Data/out_of_sample/PrismClimateDataV3_oos2.csv") #%>%
  
  # #standardize covariates based on mean and sd from in-sample data
  # dplyr::mutate(ppt = (ppt-mean(N.clim$ppt))/sd(N.clim$ppt)) %>%
  # 
  # dplyr::mutate(tmean = (tmean-mean(N.clim$tmean))/sd(N.clim$tmean)) %>%
  # 
  # dplyr::mutate(vpdmax = (vpdmax-mean(N.clim$vpdmax))/sd(N.clim$vpdmax))

N3.topo <- read.csv("cover_spread/Data/out_of_sample/topographic_data_oos2.csv") #%>%
  
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


## ----- plot climate covariates ----

## ppt

# cols
linecols <- c("in sample" = "#bec991", "oos (nearby)" = "#907350", "oos (distant)" = "#5f8971")

# median ppt
ppt_N <-  data.frame("annppt" = apply(enviro.var[,,1], MARGIN = 1, FUN = mean), 
                     "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
ppt_N2 <- data.frame("annppt" = apply(enviro.var.N2[,,1], MARGIN = 1, FUN = mean), 
                     "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
ppt_N3 <- data.frame("annppt" = apply(enviro.var.N3[,,1], MARGIN = 1, FUN = mean), 
                     "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)

(ppt.plot <- ppt_N %>% ggplot(aes(x = year, y = annppt)) + 
    geom_point(aes(col = "in sample"),size = 2) + 
    geom_line(aes(col = "in sample")) + 
    geom_line(aes(col = "oos (nearby)"), data = ppt_N2) + 
    geom_point(aes(col = "oos (nearby)"), data = ppt_N2, size = 2) +
    geom_line(aes(col = "oos (distant)"), data = ppt_N3) +
    geom_point(aes(col = "oos (distant)"), data = ppt_N3, size = 2) +
    scale_color_manual(name = "landscape", values = linecols) +
    ylab("Annual Precipitation (mm)") +
    xlab("Year") +
    theme_classic())

## mean tmean
tmean_N <-  data.frame("tmean" = apply(enviro.var[,,2], MARGIN = 1, FUN = mean), 
                       "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
tmean_N2 <- data.frame("tmean" = apply(enviro.var.N2[,,2], MARGIN = 1, FUN = mean), 
                       "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
tmean_N3 <- data.frame("tmean" = apply(enviro.var.N3[,,2], MARGIN = 1, FUN = mean), 
                       "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)

(tmean.plot <- tmean_N %>% ggplot(aes(x = year, y = tmean)) + 
    geom_point(aes(col = "in sample"),size = 2) + 
    geom_line(aes(col = "in sample")) + 
    geom_line(aes(col = "oos (nearby)"), data = tmean_N2) + 
    geom_point(aes(col = "oos (nearby)"), data = tmean_N2, size = 2) +
    geom_line(aes(col = "oos (distant)"), data = tmean_N3) +
    geom_point(aes(col = "oos (distant)"), data = tmean_N3, size = 2) +
    scale_color_manual(name = "landscape", values = linecols) +
    ylab("Mean Temperature (c)") +
    xlab("Year") +
    theme_classic())

##  heatload

# map
par(mfrow = c(1,3))

topodat %>% dplyr::select(c(x,y,heatload)) %>% 
   df2raster(.) %>% raster::plot(., zlim = c(0.3,1), col.axis="white", tck = 0)
N2.topo %>% dplyr::select(c(x,y,heatload)) %>% 
  df2raster(.) %>% raster::plot(., zlim = c(0.3,1), col.axis="white", tck = 0)
N3.topo %>% dplyr::select(c(x,y,heatload)) %>% 
  df2raster(.) %>% raster::plot(., zlim = c(0.3,1), col.axis="white", tck = 0)

# mean
mean(enviro.var[1,,4]) # n
sd(enviro.var[1,,4])

mean(enviro.var.N2[1,,3]) # n2
sd(enviro.var.N2[1,,3])

mean(enviro.var.N3[1,,3]) # n3
sd(enviro.var.N3[1,,3])

## elevation

# map 
par(mfrow = c(1,3))

topodat %>% dplyr::select(c(x,y,elev)) %>% df2raster(.) %>% raster::plot(., zlim = c(1200,1600), col.axis="white", tck = 0)
N2.topo %>% dplyr::select(c(x,y,elev)) %>% df2raster(.) %>% raster::plot(., zlim = c(1200,1600), col.axis="white", tck = 0)
N3.topo %>% dplyr::select(c(x,y,elev)) %>% df2raster(.) %>% raster::plot(., zlim = c(2100,2500), col.axis="white", tck = 0)


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

