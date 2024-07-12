## setwd to PJ_Photo/cover_spread/data in G drive
setwd("D:/.shortcut-targets-by-id/1FPlPAVacVgAROSPXMiiOGb2Takzm2241/PJ_Photo/cover_spread/Data")

library(tidyverse)
library(raster)
library(RColorBrewer)
library(rasterVis)
library(cowplot)
library(sf)

## functions
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

## load data

# load observed timeseries of cover data
N <- read.csv("RAPtreecoverData.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

N2 <- read.csv("out_of_sample/RAPtreecoverData_oos1.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

N3 <- read.csv("out_of_sample/RAPtreecoverData_oos2.csv") %>% dplyr::select(-c(x, y, X)) %>% pivot_wider(names_from = cellnum, values_from = cover) %>% dplyr::select(-year) %>% as.matrix()

# load rasters of cover data for 1986 for each landscape
N.rast <- raster("RAP_treeCover_1986.tif") %>% 
  raster2df(.) %>% dplyr::select(-cellnum) %>% df2raster(.)
N2.rast <- raster("out_of_sample/oos_raster_data/RAP_treeCover_oos1_1986.tif") %>% 
  raster2df(.) %>% dplyr::select(-cellnum) %>% df2raster(.)
N3.rast <- raster("out_of_sample/oos_raster_data/RAP_treeCover_oos2_1986.tif") %>% 
  raster2df(.) %>% dplyr::select(-cellnum) %>% df2raster(.)


## plot average cover over time
tiff("C:/Users/elise/OneDrive/Documents/UNR_Gr/chapter1/figures/study_areas_figure/landscapes_diagram.tif",width = 8,height=5.5,units="in", res=300)

# ggplot version
N.avg <- data.frame("avgcover" = rowMeans(N), "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
N2.avg <- data.frame("avgcover" = rowMeans(N2), "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)
N3.avg <- data.frame("avgcover" = rowMeans(N3), "year" = seq(1,36,1)) %>% dplyr::mutate(year = year+1985)

#linecols <- c("in sample" = "#4793AF", "oos (nearby)" = "#FFC470", "oos (far away)" = "#8B322C")
linecols <- c("in sample" = "#bec991", "oos (nearby)" = "#907350", "oos (distant)" = "#5f8971")
(b <- N.avg %>% ggplot(aes(x = year, y = avgcover)) + 
    geom_point(aes(col = "in sample"),size = 2) + 
    geom_line(aes(col = "in sample"),lwd = 1.1) + 
    geom_line(aes(col = "oos (nearby)"), data = N2.avg, lwd = 1.1) + 
    geom_point(aes(col = "oos (nearby)"), data = N2.avg, size = 2) +
    geom_line(aes(col = "oos (distant)"), data = N3.avg, lwd = 1.1) +
    geom_point(aes(col = "oos (distant)"), data = N3.avg, size = 2) +
    scale_color_manual(name = "landscape", values = linecols) +
    ylab("AVG. % TREE COVER") +
    xlab("YEAR") +
    theme_classic())
## plot mapped cover for 1986

cols <- colorRampPalette(brewer.pal(9,"YlGn"), bias = 1.75)(50)

# panel
a.1 <- levelplot(N.rast, col.regions = cols, scales = list(draw=F), margin = F, ylab = "", 
            xlab = "", colorkey = list(width = 0.75), at=seq(0,50), 
            main = "In sample")

a.1$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding',
     'main.key.padding') ] <- 1
a.1$aspect.fill <- TRUE


a.2 <- levelplot(N2.rast, col.regions = cols, scales = list(draw=F), margin = F, ylab = "", 
            xlab = "", colorkey = list(width = 0.75), at=seq(0,50), 
            main = "Out of sample (nearby)")

a.2$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding',
     'main.key.padding') ] <- 1
a.2$aspect.fill <- TRUE

a.3 <- levelplot(N3.rast, col.regions = cols, scales = list(draw=F), margin = F, ylab = "", 
            xlab = "", at=seq(0,50), colorkey = list(width = 0.75), 
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

(c <- ggplot() + theme_void())

b.c <- (plot_grid(b, c, ncol = 2, rel_widths = c(3.5,1)))

plot_grid(a,b.c, nrow = 2)
 
dev.off()

