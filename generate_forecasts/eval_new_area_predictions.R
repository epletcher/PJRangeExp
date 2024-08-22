# you will need to LOAD WORKSPACE FILES WITH OOS MODEL PREDICTIONS
load("FILEPATH/35y_insample_predictions_5y_average_initial.RData")
load("FILEPATH/35y_OOS_1_near_predictions_5y_average_initial.RData")
load("FILEPATH/35y_OOS_2_far_predictions_5y_average_initial.RData")

# load libraries
library(tidyverse)
library(cowplot)

# load functions
normalized_rmse <- function(rmsedat, N) {
  nrmse <- rmsedat/(max(N)-min(N))
  return(nrmse)
}
# --------------- PULL RMSE'S --------------------
# calculate RMSE all together
# calculate RMSE in 5-yr chunks
# choose a set of predictions for a specific study area

# ** Choose which landscape you would like to plot rmse's here for

# N
BASE = for.base.N$rmseTotOut
TOPO = for.topo.N$rmseTotOut
CLIM = for.clim.N$rmseTotOut
TOPOCLIM = for.topoclim.N$rmseTotOut

# # N2
# BASE = for.base.N2$rmseTotOut
# TOPO = for.topo.N2$rmseTotOut
# CLIM = for.clim.N2$rmseTotOut
# TOPOCLIM = for.topoclim.N2$rmseTotOut
# 
# # N3
# BASE = for.base.N3$rmseTotOut
# TOPO = for.topo.N3$rmseTotOut
# CLIM = for.clim.N3$rmseTotOut
# TOPOCLIM = for.topoclim.N3$rmseTotOut
#

# ---------- PLOT RMSE'S ACROSS MODELS TOGETHER ----------
## generate figure 4 (A, B, and C separately)
obs = N # ** study area here, N, N2 or N3

# base dataframe of rmse's
base <- BASE %>% 
  normalized_rmse(.,obs) %>% 
  as.data.frame() %>% 
  mutate('year' = seq(1,36,1)) %>% 
  pivot_longer(starts_with('V'), values_to = 'base', names_to = 'iteration') %>%
  filter(year != 1) %>%
  mutate(year_int = if_else(year<7,1,
                            if_else(year<12, 2,
                                    if_else(year<17, 3,
                                            if_else(year<22, 4,
                                                    if_else(year<27, 5,
                                                            if_else(year<32, 6, 7)))))))
  
  # for year_int, if any 'base'<-500 & >500 
  to.remove <- base %>% filter(base > 500) %>% 
      dplyr::select(year_int) %>% pull() %>% unique()

  # filter(year_int %in% 
  base <- base %>% filter(!year_int %in% to.remove)

# climate dataframe of rmse
clim <- CLIM %>% 
  normalized_rmse(.,obs) %>% 
  as.data.frame() %>% 
  mutate('year' = seq(1,36,1)) %>% 
  pivot_longer(starts_with('V'), values_to = 'clim', names_to = 'iteration') %>%
  filter(year != 1) %>%
  mutate(year_int = if_else(year<7,1,
                            if_else(year<12, 2,
                                    if_else(year<17, 3,
                                            if_else(year<22, 4,
                                                    if_else(year<27, 5,
                                                            if_else(year<32, 6, 7))))))) 
  
  # for year_int, if any 'clim'<-500 & >500 
  to.remove <- clim %>% filter(clim > 500) %>% 
      dplyr::select(year_int) %>% pull() %>% unique()

  # filter(year_int %in% 
  clim <- clim %>% filter(!year_int %in% to.remove)


# topo dataframe of rmse
topo <- TOPO %>% 
  normalized_rmse(.,obs) %>% 
  as.data.frame() %>% 
  mutate('year' = seq(1,36,1)) %>% 
  pivot_longer(starts_with('V'), values_to = 'topo', names_to = 'iteration') %>%
  filter(year != 1) %>%
  mutate(year_int = if_else(year<7,1,
                            if_else(year<12, 2,
                                    if_else(year<17, 3,
                                            if_else(year<22, 4,
                                                    if_else(year<27, 5,
                                                            if_else(year<32, 6, 7))))))) 

  # for year_int, if any 'topo'<-500 & >500 
  to.remove <- topo %>% filter(topo > 500) %>% 
  dplyr::select(year_int) %>% pull() %>% unique()

  # filter(year_int %in% 
  topo <- topo %>% filter(!year_int %in% to.remove)


# topoclim dataframe of rmse
topoclim <- TOPOCLIM %>% 
  normalized_rmse(.,obs) %>% 
  as.data.frame() %>% 
  mutate('year' = seq(1,36,1)) %>% 
  pivot_longer(starts_with('V'), values_to = 'topoclim', names_to = 'iteration') %>%
  filter(year != 1) %>%
  mutate(year_int = if_else(year<7,1,
                            if_else(year<12, 2,
                                    if_else(year<17, 3,
                                            if_else(year<22, 4,
                                                    if_else(year<27, 5,
                                                            if_else(year<32, 6, 7)))))))
  
  # for year_int, if any 'topoclim'<-500 & >500 
  to.remove <- topoclim %>% filter(topoclim > 500) %>% 
    dplyr::select(year_int) %>% pull() %>% unique()

  # filter(year_int %in% 
  topoclim <- topoclim %>% filter(!year_int %in% to.remove)

# combine rmse dataframes
dat <- base %>% full_join(., clim) %>% 
  full_join(., topo) %>%
  full_join(., topoclim) %>%
  pivot_longer(c(base, clim, topo, topoclim), names_to = 'model',
               values_to = 'rmse') %>%
  # pivot_longer(c(base, clim), names_to = 'model', values_to = 'rmse') %>%
  mutate(year_int = if_else(year_int==1,"2 - 6",
                            if_else(year_int==2,"7 - 11",
                                    if_else(year_int==3,"12 - 16",
                                            if_else(year_int==4,"17 - 21",
                                                    if_else(year_int==5,"22 - 26",
                                                            if_else(year_int==6,"27 - 31","32 - 36"))))))) %>%
  mutate(year_int = fct_relevel(year_int, 
                            "2 - 6", "7 - 11", "12 - 16", 
                            "17 - 21", "22 - 26", "27 - 31", 
                            "32 - 36"))
## plot

# colors
group.cols <- c("#F8766D","#00BFC4","#C77Cff","#7CAE00")
#group.cols <- c("#C77Cff","#7CAE00") # for plotting topos only

(rmse.plot <- dat %>% 
    #filter(model=='topo'|model=='topoclim') %>% # for plotting topos only
    ggplot(aes(x = year_int, y = rmse, col = model), col = cols) + 
    geom_boxplot(outlier.shape = NA) + 
    labs(y = 'NRMSE', x = "5-year interval") +
    scale_y_continuous(limits = c(0, 0.6)) +
    #scale_y_continuous(limits = c(0,30)) +
    scale_color_manual(values=group.cols) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key = element_rect(fill = "white"),
      panel.background = element_rect(linetype = "solid",fill = NA),
      panel.border = element_rect(linetype = "solid", fill = NA),
      panel.grid.major = element_line(colour = "lightgrey", linewidth = .4)
      )
    ) # NA's are the shorter length of iterations for 2 of the models

## **** Create table of all rmse values (in supplement)
rmse.35y <- dat %>% 
  group_by(model,year_int) %>%
  #group_by(model) %>%
  summarise(median = median(rmse, na.rm = T), lower_ci = quantile(rmse, na.rm = T, 0.05), upper_ci = quantile(rmse, na.rm = T, 0.95))

## save correct study areas to file
# ** update the file path for each study landscape
# write.csv(rmse.35y,"FILEPATH/rmses_for_model_studyarea_comparisons/in_sample_35y_nrmse_avg_allyears.csv")


# ------------------ Predicted median and credible intervals ------------
# Base and clim model onlu

# ** Choose which landscape you would like to calculate median and CI predictions for here:

# N
base.pred = for.base.N$predOut
topo.pred = for.topo.N$predOut
clim.pred = for.clim.N$predOut
topoclim.pred = for.topoclim.N$predOut

# # N2
# base.pred = for.base.N2$predOut
# topo.pred = for.topo.N2$predOut
# clim.pred = for.clim.N2$predOut
# topoclim.pred = for.topoclim.N2$predOut

# # N3
# base.pred = for.base.N3$predOut
# topo.pred = for.topo.N3$predOut
# clim.pred = for.clim.N3$predOut
# topoclim.pred = for.topoclim.N3$predOut

## median predictions and 90% credible intervals

# base
med.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = median)
low.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
up.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = quantile, 0.95)

# save median pixel-year predicted cover for mapping forecasts
# ** change file path depending on landscape
write.csv(med.pred.base, "FILEPATH/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")

# clim
med.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = median)
low.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
up.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = quantile, 0.95)

# ----------------- Plot predictions over time ------------------
## create a plotting function
plot_pix_gg <- function(pix, obs) {
  px <- as.character(pix)
  pxv <- paste("V",as.character(pix), sep = "")
  
  #base
  med.base <- med.pred.base %>% as.data.frame() 
  med.base$year <- row.names(med.base) 
  med <- med.base %>% mutate(year = as.numeric(year)+1985) %>% 
    pivot_longer(!year, values_to = "med.pred.base", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1986) # remove first year (true value, not part of the forecast)
  
  low.base <- low.pred.base %>% as.data.frame() # low
  low.base$year <- row.names(low.base)
  low <- low.base %>% mutate(year = as.numeric(year)+1985) %>% 
    pivot_longer(!year, values_to = "low.pred.base", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1) # remove first year (true value, not part of the forecast)
  
  up.base <- up.pred.base %>% as.data.frame() # high
  up.base$year <- row.names(up.base)
  up <- up.base %>% mutate(year = as.numeric(year)+1985) %>% 
    pivot_longer(!year, values_to = "up.pred.base", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1986) # remove first year (true value, not part of the forecast)
  
  #clim
  med.clim <- med.pred.clim %>% as.data.frame()
  med.clim$year <- row.names(med.clim)
  med.clim <- med.clim %>% mutate(year = as.numeric(year)+1985) %>%
    pivot_longer(!year, values_to = "med.pred.clim", names_to = "pixel") %>%
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1986) # remove first year (true value, not part of the forecast)

  low.clim <- low.pred.clim %>% as.data.frame() # low
  low.clim$year <- row.names(low.clim)
  low.clim <- low.clim %>% mutate(year = as.numeric(year)+1985) %>%
    pivot_longer(!year, values_to = "low.pred.clim", names_to = "pixel") %>%
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1986) # remove first year (true value, not part of the forecast)

  up.clim <- up.pred.clim %>% as.data.frame()  # high
  up.clim$year <- row.names(up.clim)
  up.clim <- up.clim %>% mutate(year = as.numeric(year)+1985) %>%
    pivot_longer(!year, values_to = "up.pred.clim", names_to = "pixel") %>%
    filter(pixel == pxv) %>% dplyr::select(-pixel) %>%
    filter(year != 1986) # remove first year (true value, not part of the forecast)

  # observed 
  Ndat <- obs %>% as.data.frame() %>% rownames_to_column("year") %>% 
    pivot_longer(!year, values_to = "obs", names_to = "pixel") %>% 
    filter(pixel == px) %>% 
    mutate(year = as.numeric(year)+1985) %>% 
    dplyr::select(-pixel)
  
  # combine predictions and observations into one dataframe
  plot.dat <- left_join(Ndat, med) %>% left_join(.,low) %>% left_join(.,up) %>% left_join(.,up.clim) %>% left_join(.,low.clim) %>% left_join(.,med.clim)
  
  ## plot
  
  # color palette
  cols = c('base' = '#F8766D', 'climate-only' = '#00BFC4')
 
  (plot.dat %>% 
    ggplot(aes(x = year, y = obs)) + 
    geom_point(aes(x = year, y = obs, shape = "observations"), fill = 'black') + geom_line(aes(x = year, y = obs), lty = 2) + 
    geom_ribbon(aes(ymin = low.pred.base, ymax = up.pred.base, fill = "base"), 
                alpha=0.4) + 
    geom_ribbon(aes(ymin = low.pred.clim, ymax = up.pred.clim, fill = "climate-only"), 
                alpha=0.3) + 
    geom_line(aes(x = year, y = med.pred.base, col = "base"), lwd = 1.5) + 
    geom_line(aes(x = year, y = med.pred.clim, col = "climate-only"), lwd = 1.5) +
    
    scale_y_continuous(limits = c(-5, 25)) + # good scale for N
    #scale_y_continuous(limits = c(-5, 35)) + # good scale for N2
    #scale_y_continuous(limits = c(-5, 80)) + # good scale for N3 
    scale_color_manual(name = '', values = cols) +
    scale_fill_manual(name = '', values = cols) +
    scale_shape_manual(name = '', values = 21, labels = 'observations') +
    
    labs(x = "YEAR", y = "TREE COVER (%)") + 
    
    theme(
      legend.justification = c("right", "top"),
      legend.position = c(.6, .99),
      legend.title=element_blank(),
      legend.spacing.y = unit(-.2, 'cm'),
      legend.box.margin = margin(4,4,4,4),
      legend.key=element_blank(),
      legend.key.size = unit(0.75, 'cm'),
      legend.text=element_text(size=13),
      text = element_text(size=18),
      panel.background = element_rect(linetype = "solid",fill = NA),
      panel.border = element_rect(linetype = "solid", fill = NA))
  )
}

## plot low, medium, and high cover all-together (Figure 5)

# ** Choose which landscape you would like to calculate median and CI predictions for here:

# N
(three.panel <- plot_grid(plot_pix_gg(707, N) + labs(x = '', y = ''), 
          plot_pix_gg(152, N) + theme(legend.position = "none") + labs(x = '', y = ''), 
          plot_pix_gg(120, N) + theme(legend.position = "none") + labs(x = '', y = ''), ncol = 3))

# # N2
# (three.panel <- plot_grid(plot_pix_gg(711, N2) + theme(legend.position = "none") + labs(x = ''), 
#                           plot_pix_gg(22, N2) + theme(legend.position = "none") + labs(x = '', y = ''), 
#                           plot_pix_gg(2050, N2) + theme(legend.position = "none") + labs(y = '', x = ''), ncol = 3))

# # N3
# (three.panel <- plot_grid(plot_pix_gg(1203, N3) + theme(legend.position = "none") + labs(x = '', y = ''), 
#           plot_pix_gg(1205, N3) + theme(legend.position = "none") + labs(y = ''), 
#           plot_pix_gg(1213, N3) + theme(legend.position = "none") + labs(x = '', y = ''), ncol = 3))