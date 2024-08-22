# load packages
library(maptools)
library(raster)
library(rgdal)
library(spatstat) # use pairdist to calc. distance matrix
library(tidyverse)
library(abind)

## --------- load data here ---------------
# Need tree cover and environmental covariate data converted from raster, to dataframes

# create raster to calc distance
r <- raster(nrows = 42, ncols = 56, xmn = -120.5686, xmx = -120.5535, ymn = 43.87138, ymx = 43.8827)

# tree cover data (raster converted to dataframe)
treedat <- read.csv("RAPtreecoverData.csv")

# clim data (raster converted to dataframe)
climdat <- read.csv("PrismClimateDataV3.csv") %>% dplyr::mutate(across(c("ppt", "tmean", "vpdmax"), scale))

# topo data (read in and standardize)
topodat <- read.csv("topographic_data.csv") %>% dplyr::mutate(across(c("heatload", "elev"), scale))

#--------- reformat data ------------
# matrix of tree cover data cols = pixel or cell number and rows = years
N <- treedat %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

# 3-dimensional array of environmental covariates by cellnum, year, and variable

ppt <- climdat %>% dplyr::select(c(cellnum, year, ppt)) %>% pivot_wider(names_from = cellnum, values_from = ppt) %>% dplyr::select(-year) %>% as.matrix()

tmean <- climdat %>% dplyr::select(c(cellnum, year, tmean)) %>% pivot_wider(names_from = cellnum, values_from = tmean) %>% dplyr::select(-year) %>% as.matrix()

vpdmax <- climdat %>% dplyr::select(c(cellnum, year, vpdmax)) %>% pivot_wider(names_from = cellnum, values_from = vpdmax) %>% dplyr::select(-year) %>% as.matrix()

heatload <- topodat %>% dplyr::select(c(cellnum, year, heatload)) %>% pivot_wider(names_from = cellnum, values_from = heatload) %>% dplyr::select(-year) %>% as.matrix()

elev <- topodat %>% dplyr::select(c(cellnum, year, elev)) %>% pivot_wider(names_from = cellnum, values_from = elev) %>% dplyr::select(-year) %>% as.matrix()

#stack
enviro.var <- abind(ppt, tmean, vpdmax, heatload, elev, along = 3)

#------ Build distance matrix -------
points <- rasterToPoints(r)

# specify window 
rwin <- as(extent(r), 'SpatialPolygons')
rowin <- as.owin(rwin)

# calculate pairwise distances, and output into a matrix 
D <- pairdist(as.ppp(X = points, W = rowin))*1000
# Dsc = D*1000 # rescale the distance matrix so values are close to 1-100
Dsq = D*D # square D for input into dispersal equation