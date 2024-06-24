# This script is modifying and paring down the 'forecasts_2_rasters' script to plot/map base model predictions across each study area (in addition to initial cover values)

# **** in progress !!!

# load packages
library(raster)
library(sf)
library(rasterVis)
library(RColorBrewer)
library(tidyverse)
library(cowplot)

# --------- functions -----------
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

# function for calculating bias for a raster forecast
biasfunc <- function(pred,obs) {
  n <- dim(pred)[3]
  bias <- (sum(pred-obs))/n
  return(bias)
}
# -------- load data ------------

## setwd to PJ_Photo/cover_spread/data in G drive
setwd("E:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data")

## load observed cover for each landscape (tabular)
N <- read.csv("RAPtreecoverData.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()
  
N2 <- read.csv("out_of_sample/RAPtreecoverData_oos1.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

N3 <- read.csv("out_of_sample/RAPtreecoverData_oos2.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()


## load observed cover for 1986 for each landscape (raster)

N.rast <- raster("RAP_treeCover_1986.tif") %>% 
  raster::aggregate(., fact = 30)

N2.rast <- raster("out_of_sample/oos_raster_data/RAP_treeCover_oos1_1986.tif") %>% 
  raster::aggregate(., fact = 30)

N3.rast <- raster("out_of_sample/oos_raster_data/RAP_treeCover_oos2_1986.tif") %>% 
  raster::aggregate(., fact = 30)


## load tabular predictions for each landscape
# Forecasted cover for each pixel (TABULAR). Columns are pixel number, rows are year

# N
N.pp <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")

# N2
N2.pp <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_1_near_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")

# N3
N3.pp <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")

# -------- calculate raster time series of bias ----------
# *** this need to be checked to make sure its working right
mapping_bias <- function(obs.name, pred, obs) { # obs is character name of landscape, pred is data.frame obect name of loaded predictions
 
  ## 1. set file paths for location of observed raster data for each respective landscape
  if(obs.name=="N") {
    # load list (all RAP rasters)
    RAP.list <- list.files("E:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data", pattern = 'RAP_treeCover*', full.names = TRUE)
  }
  
  if(obs.name=="N2") {
    # load list (all RAP rasters)
    RAP.list <- list.files("E:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data/out_of_sample/oos_raster_data", pattern = 'RAP_treeCover_oos1_*', full.names = TRUE)
  }
  
  if(obs.name=="N3") {
    # load list (all RAP rasters)
    RAP.list <- list.files("E:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data/out_of_sample/oos_raster_data", pattern = 'RAP_treeCover_oos2_*', full.names = TRUE)
  }
  
  # load observed data rasters 
  tru.rasters <- lapply(lapply(RAP.list, FUN = raster::raster), FUN = raster::aggregate, fact = 30)
  # stack
  truthstack <- raster::stack(tru.rasters)
  
  ## 2. Extract xy info from observed rasters for adding to forecast rasters
  # convert RAP rasters to df
  tree.df <- raster2df(truthstack) %>% 
    pivot_longer(cols = !c("cellnum", "x", "y"), names_to = "years", values_to = "cover") %>%
    separate(years, c("name", "year"), sep = "E.") %>% 
    dplyr::select(-"name") %>% 
    dplyr::mutate('year' = as.numeric(year)+1985) %>%
    pivot_wider(names_from = year, values_from = cover)
  
  # Extract values for 1986 and xy data to join with forecasted values
  start.xy <- tree.df %>% dplyr::select(c("cellnum", "x", "y", "1986")) 
  
  ## 3. Reformat forecasted cover data and join w/ start.xy dataframe to recover x,y data
  t.pp <- t(pred) # row = pixel, column = year (transpose matrix)
  t.pp <- t.pp[-1,]
  
  # add years as column names
  oldnames <- t.pp %>% as.data.frame() %>% colnames()
  tree.for <- t.pp %>% as.data.frame() %>% rename_with(~ as.character(1986:2021)[which(oldnames == .x)], .cols = oldnames)
  
  # join
  forecast.dat <- cbind(start.xy, tree.for) %>% dplyr::select(-c("1986","cellnum"))
  forecast.dat.cells <- cbind(start.xy, tree.for) %>% dplyr::select(-c("1986"))
  
  # Create an array where each row is pixel, retain x,y, and year columns, 3rd dimension is year
  cols <- c("x", "y", "cover")
  for.array <- array(NA, c(length(forecast.dat$x),3, 36), dimnames = list(NULL,cols,NULL))
  
  years <- 1986:2021
  for(t in 1:36) {
    for.array[,,t] <- forecast.dat %>% dplyr::select(c(x, y, as.character(years[t]))) %>% as.matrix()
  }
  
  # Convert forecast array into a raster stack 
  for.rasters <- list()
  for(t in 1:36) {
    for.rasters[[t]] <- df2raster(as.data.frame(for.array[,,t]))
  }
  forecast.stack <- raster::stack(for.rasters)
  
  ## 4. calculate forecast bias and nrmse for each pixel
  n <- 35 # minus one because remove 1st year
  
  bias <- calc((forecast.stack[[2:36]] - truthstack[[2:36]]), sum)/n # bias based on observed, exclude 1st year
  
  #mae <- calc(abs(forecast.stack[[2:36]] - truthstack[[2:36]]), sum)/n
  
  #nmae <- calc((abs(forecast.stack[[2:36]] - truthstack[[2:36]])/(max(obs)-min(obs))), sum)/n
  
  #nrmse <- sqrt(calc(((forecast.stack[[2:36]] - truthstack[[2:36]])^2), sum)/n)/(max(obs)-min(obs))
  
  #nrmse2 <- sqrt(calc((((forecast.stack[[2:36]] - truthstack[[2:36]])^2)/(max(obs)-min(obs))), sum)/n)
  
  # error between predicted last year's cover and a 5-yr average
  avgerr <- forecast.stack[[36]] - calc(truthstack[[31:36]], mean)
  
  eval.raster <- raster::stack(bias, avgerr) # choose which metrics to print here
  
  return(eval.raster)
}

# calc mapped bias for each landscape
N.bias <- mapping_bias("N", N.pp, N)

N2.bias <- mapping_bias("N2", N2.pp, N2)

N3.bias <- mapping_bias("N3", N3.pp, N3) 

col1 = colorRampPalette(rev(brewer.pal(7,"RdBu")))(50)

# -------- base plot ----------------
par(mfrow = c(3,2),mai = c(1,0.1,0.1,0.1))
plot(N.bias[[1]], col = col1, zlim = c(-16, 16), yaxt="n", xaxt="n") # bias
plot(N.bias[[2]], col = col1, zlim = c(-24, 24), yaxt="n", xaxt="n") # 5-yr average error

plot(N2.bias[[1]], col = col1, zlim = c(-16, 16), yaxt="n", xaxt="n")
plot(N2.bias[[2]], col = col1, zlim = c(-24, 24), yaxt="n", xaxt="n")

plot(N3.bias[[1]], col = col1, zlim = c(-16, 16), yaxt="n", xaxt="n")
plot(N3.bias[[2]], col = col1, zlim = c(-24, 24), yaxt="n", xaxt="n")


# -------- level plot --------------

# ## plot average cover over time
tiff("R:/Shriver_lab/PJspread/figures/mapping_error_across_landscapes.tif",width = 8,height=2.75,units="in", res=300)

# panel
a.1 <- levelplot(N.bias[[2]], col.regions = col1, scales = list(draw=F), margin = F, ylab = "", 
                 xlab = "", colorkey = list(width = 0.75), at=seq(-24,24), 
                 main = "In sample")

a.1$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding',
     'main.key.padding') ] <- 1
a.1$aspect.fill <- TRUE


a.2 <- levelplot(N2.bias[[2]], col.regions = col1, scales = list(draw=F), margin = F, ylab = "", 
                 xlab = "", colorkey = list(width = 0.75), at=seq(-24,24), 
                 main = "Out of sample (nearby)")

a.2$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding',
     'main.key.padding') ] <- 1
a.2$aspect.fill <- TRUE

a.3 <- levelplot(N3.bias[[2]], col.regions = col1, scales = list(draw=F), margin = F, ylab = "", 
                 xlab = "", colorkey = list(width = 0.75), at=seq(-24,24),
                 main = "Out of sample (far away)")

a.3$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding',
     'main.key.padding') ] <- 1
a.3$aspect.fill <- TRUE

(a <- plot_grid(a.1,a.2,a.3, rel_widths = c(1,1,1), ncol = 3))


dev.off()

