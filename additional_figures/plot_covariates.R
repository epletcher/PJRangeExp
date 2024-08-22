## Code to generate supplemental figures 

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
climdat <- read.csv("FILEPATH/PrismClimateDataV3.csv") 
#%>% plyr::mutate(across(c("ppt", "tmean", "vpdmax"), scale))

# topo data (read in and standardize)
topodat <- read.csv("FILEPATH/topographic_data.csv") 
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

## ----- N2 -----
N2.clim <- read.csv("FILEPATH/PrismClimateDataV3_oos1.csv") #%>% 
  
N2.topo <- read.csv("FILEPATH/topographic_data_oos1.csv") #%>% 
  
# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N2.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N2.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N2.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N2.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N2 <- abind(ppt, tmean, heatload, elev, along = 3)

# ## N3
N3.clim <- read.csv("FILEPATH/PrismClimateDataV3_oos2.csv") #%>%
  
N3.topo <- read.csv("cover_spread/Data/out_of_sample/topographic_data_oos2.csv") #%>%
  
# 3-dimensional array of environmental covariates by cellnum, year, and variable
ppt <- N3.clim %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- N3.clim %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

heatload <- N3.topo %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- N3.topo %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var.N3 <- abind(ppt, tmean, heatload, elev, along = 3)

## ----- plot climate covariates ----
## generates figure SA2

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