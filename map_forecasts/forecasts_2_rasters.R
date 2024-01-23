# This script use outputs predicted outputs from different Samplers to map predicted shifts in cover in raster format. Outputs are saved as TIFFs in the "RangeExp_model_predictions" folder.
# edited 8/1/2023

# load packages
library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(moveVis)
library(tidyverse)

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

# ------------ load data -------------------
# 1.) load model workspace to extract latent state for N
# load model objects generated from:
# 'eval_model_insample.R'
# Latent states and parameter estimates are post burnin and thinning 

mod1 <- readRDS(file = "R:/Shriver_Lab/PJspread/evaluate_in_sample/model_objects_for_forecasting/base_model.rds") # base
med.lat.cover <- apply(mod1$Nlat, MARGIN = c(1,2), FUN = median)

# 2.) Forecasted cover for each pixel (TABULAR). Columns are pixel number, rows are year
pp <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv") # model version here

# 3.) Observed time series data (TABULAR)
tru.cover.dat <- read.csv("G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data/RAPtreecoverData.csv")

# 4.) Observed spatial tree cover data 1986-2021 (RASTER)
# load list (all RAP rasters)
RAP.list <- list.files("G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data", pattern = 'RAP_treeCover*', full.names = TRUE)
  # load
  tru.rasters <- lapply(lapply(RAP.list, FUN = raster::raster), FUN = raster::aggregate, fact = 30)
  # stack
  truthstack <- raster::stack(tru.rasters)

# convert RAP rasters to df
tree.df <- raster2df(truthstack) %>% 
  pivot_longer(cols = !c("cellnum", "x", "y"), names_to = "years", values_to = "cover") %>%
  separate(years, c("name", "year"), sep = "E.") %>% 
  dplyr::select(-"name") %>% 
  dplyr::mutate('year' = as.numeric(year)+1985) %>%
  pivot_wider(names_from = year, values_from = cover)

# Extract values for 1986 and xy data to join with forecasted values
start.xy <- tree.df %>% dplyr::select(c("cellnum", "x", "y", "1986")) 

# --------------- Recover forecasted XY data ----------------------
# Reformat forecasted cover data and join w/ start.xy dataframe to recover x,y data
t.pp <- t(pp) # row = pixel, column = year (transpose matrix)
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
dim(for.array)[3]
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

# --------------- Recover latent XY data ----------------
# Reformat forecasted cover data and join w/ start.xy dataframe to recover x,y data
t.lat <- t(med.lat.cover) # row = pixel, column = year (transpose matrix)

# add years as column names
oldnames <- t.lat %>% as.data.frame() %>% colnames()
tree.lat <- t.lat %>% as.data.frame() %>% rename_with(~ as.character(1986:2016)[which(oldnames == .x)], .cols = oldnames)

# join
lat.dat <- cbind(start.xy, tree.lat) %>% dplyr::select(-c("1986", "cellnum"))

# Create an array where each row is pixel, retain x,y, and year columns, 3rd dimension is year
cols <- c("x", "y", "cover")
lat.array <- array(NA, c(length(lat.dat$x),3, 31), dimnames = list(NULL,cols,NULL))

years <- 1986:2016
for(t in 1:31) {
  lat.array[,,t] <- lat.dat %>% dplyr::select(c(x, y, as.character(years[t]))) %>% as.matrix()
}

# Convert forecast array into a raster stack 
lat.rasters <- list()
for(t in 1:31) {
  lat.rasters[[t]] <- df2raster(as.data.frame(lat.array[,,t]))
}
lat.stack <- raster::stack(lat.rasters)

# -------- Plot true timeseries of cover --------------
animate(truthstack, pause=0.5, n = 1, main = "observed change in cover", col = brewer.pal(9,"YlGn")) # plays tree cover each year as a loop

cols <- colorRampPalette(brewer.pal(9,"YlGn"))
# levelplot(truthstack, col.regions=cols, names.attr = as.character(1986:2021)) # add names with names.attr

# subset and plot only 1990, 2000, 2010, 2020
subras <- c(5,15,25,36)

# plot in grid
levelplot(truthstack[[subras]], col.regions=colorRampPalette(brewer.pal(9,"YlGn")), names.attr = as.character(c(1990,2000,2010,2020)), scales = list(draw=F))

# ---------- plot latent timeseries of cover ----------
sublat <- c(5,15,25)
levelplot(lat.stack[[sublat]], col.regions=cols, names.attr = as.character(c(1990,2000,2010)))

# -------- Plot sequential cover Forecast -------------
# plot forecast
levelplot(forecast.stack, col.regions=colorRampPalette(brewer.pal(9,"YlGn")), names.attr = as.character(1986:2021))

# subset and plot only 1990, 2000, 2010, 2020
# this plot works best if you make plot window tall and skinny
levelplot(
    subset(raster::stack(truthstack[[subras]],forecast.stack[[subras]]), c(1,5,2,6,3,7,4,8)), #reorder rasters
    col.regions = colorRampPalette(brewer.pal(9,"YlGn")), 
    xlab = "Longitude", ylab = "Latitude", 
    main = "Observed    Predicted      ", 
    zlim = c(0,30),
    names.attr = as.character(c(1990,1990,2000,2000,2010,2010,2021,2021)), 
    scales = list(draw = F))

# #animate
raster::animate(forecast.stack, pause=0.5, n = 1, main = "tree cover forecast", col = brewer.pal(9,"YlGn"), zlim = c(0,25)) # plays tree cover each year as a loop

# --------- Plot forecasted change in cover -------------
# forecasted change in cover based of first year
predchangestack <- forecast.stack-truthstack[[1]]
obschangestack <- truthstack-truthstack[[1]]
latchangestack <- lat.stack-truthstack[[1]]
# plot forecast
cols2 <- c(rev(brewer.pal(9,"Blues")), "white", brewer.pal(9,"Reds"))

#
#
#
#***** try with out level plot!!! ******
par(mfrow = c(1,4))
plot(obschangestack[[4]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(obschangestack[[14]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(obschangestack[[24]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(obschangestack[[34]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))

par(mfrow = c(1,4))
plot(predchangestack[[4]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(predchangestack[[14]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(predchangestack[[24]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))
plot(predchangestack[[34]], col = rev(brewer.pal(7,"RdBu")), zlim = c(-15, 15))

#
#
#
#
breakpoints <- c(seq(-20,-1,1), 0, seq(1,20,1))

# subset and plot only 1990, 2000, 2010, 2020
levelplot(subset(raster::stack(obschangestack[[subras]],predchangestack[[subras]]), c(1,5,2,6,3,7,4,8)), col.regions=cols2, zlim = c(-10,10), main = "Observed    Predicted      ",names.attr = as.character(c(1990,1990,2000,2000,2010,2010,2020,2020)), xlab = "Longitude", ylab = "Latitude", scales = list(draw = F))

## #animate
raster::animate(predchangestack, pause=0.5, n = 1, col = rev(brewer.pal(7,"RdBu"))) # plays tree cover each year as a loop


# -------- plot Bias in cover forecast ------
# error based on latent state
#biastack <- (forecast.stack[[1:31]] - lat.stack[[1:31]]) # bias based on latent
biastack <- (forecast.stack[[1:36]] - truthstack[[1:36]]) # bias based on observed

## #animate
raster::animate(biastack, pause=0.5, n = 1, main = "forecast bias", col = rev(brewer.pal(7,"RdBu"))) # plays tree cover each year as a loop


subd <- sublat # for including latent
#subd <- subras # subset for w/o latent
#subd <- c(15,25,35) # only 3 years


cols3 <- colorRampPalette(rev(brewer.pal(9,"RdBu")))
breakp <- seq(-15,15,1)
levelplot(biastack[[subd]], col.regions=cols3, breaks = breakp, names.attr = as.character(c(1990,2000,2010)))

# # plot all maps of cover change together (without latent)
# # subset and plot only 1990, 2000, 2010, 2020 
# levelplot(subset(raster::stack(obschangestack[[subd]],predchangestack[[subd]],biastack[[subd]]), c(1,5,9,2,6,10,3,7,11,4,8,12)), col.regions=cols2, breaks = breakpoints,names.attr = as.character(c(1990,1990,1990,2000,2000,2000,2010,2010,2010,2020,2020,2020)), xlab = "Longitude", ylab = "Latitude", scales = list(draw = F))

# # plot all maps of cover change together (without latent) ONLY 2000,2010,2020
# (mapped_cover <- levelplot(subset(raster::stack(obschangestack[[subd]],predchangestack[[subd]],biastack[[subd]]), c(1,4,7,2,5,8,3,6,9)), col.regions=cols2, breaks = breakpoints,names.attr = as.character(c(2000,2000,2000,2010,2010,2010,2020,2020,2020)), xlab = "Longitude", ylab = "Latitude", scales = list(draw = F)))

# plot all maps of cover change together (with latent)
# subset and plot only 1990, 2000, 2010 
levelplot(subset(raster::stack(obschangestack[[sublat]],latchangestack[[sublat]],predchangestack[[sublat]],biastack[[sublat]]), c(1,4,7,10,2,5,8,11,3,6,9,12)), col.regions=cols2, breaks = breakpoints,main = "Obs.          Latent        Pred.       Error   ",names.attr = as.character(c(1990,1990,1990,1990,2000,2000,2000,2000,2010,2010,2010,2010)), xlab = "Longitude", ylab = "Latitude", scales = list(draw = F))

ggsave(filename = "G:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Figures/change_in_cover_mapped_base_NOLATENT.jpeg", plot = mapped_cover, dpi = 600, width = 8, height = 6, units = "in")
# ---------- Save forecasts and cover change in raster format
# predicted cover change
writeRaster(predchangestack, "RangeExp_model_predictions/NoObs_Base_mapped_pred_change.tif") # change based on model version here

# error
writeRaster(biastack, "RangeExp_model_predictions/NoObs_Base_mapped_pred_error.tif") # change based on model version here

# latent cover change
writeRaster(latchangestack, "RangeExp_model_predictions/Base_mapped_latent_change.tif") # change based on model version here

# #observed cover change
# writeRaster(obschangestack, "RangeExp_model_predictions/Mapped_obs_change.tif") 
