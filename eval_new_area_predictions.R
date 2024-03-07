# LOAD WORKSPACE FILES WITH OOS MODEL PREDICTIONS

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

# # N
# BASE = for.base.N$rmseTotOut
# TOPO = for.topo.N$rmseTotOut
# CLIM = for.clim.N$rmseTotOut
# TOPOCLIM = for.topoclim.N$rmseTotOut

# N2
BASE = for.base.N2$rmseTotOut
TOPO = for.topo.N2$rmseTotOut
CLIM = for.clim.N2$rmseTotOut
TOPOCLIM = for.topoclim.N2$rmseTotOut

# # N3
# BASE = for.base.N3$rmseTotOut
# TOPO = for.topo.N3$rmseTotOut
# CLIM = for.clim.N3$rmseTotOut
# TOPOCLIM = for.topoclim.N3$rmseTotOut

# ------------- PLOT RMSE FOR EACH MODEL -------------

# Plot RMSE across 5 year intervals
plot_rmse <- function(dat,N) {
  
  dat %>% normalized_rmse(.,N) %>% as.data.frame() %>% mutate('year' = seq(1,36,1)) %>% 
    pivot_longer(starts_with('V'), values_to = 'yearly_rmse', names_to = 'iteration') %>%
    filter(year != 1) %>%
    mutate(year_int = if_else(year<7,1,
                              if_else(year<12, 2,
                                      if_else(year<17, 3,
                                              if_else(year<22, 4,
                                                      if_else(year<27, 5,
                                                              if_else(year<32, 6, 7))))))) %>% 
    ggplot(aes(x = year_int, y = yearly_rmse, group = year_int)) + 
    geom_boxplot(col = 'aquamarine3', outlier.shape = NA) +
    labs(y = 'RMSE', x = "5 year chunks out") +
    theme_bw()
  
}

# models that work
plot_rmse(BASE, N2)
plot_rmse(CLIM, N2)
# models that don't (for all areas)
plot_rmse(TOPOCLIM, N2)
plot_rmse(TOPO, N2)

# ---------- PLOT RMSE'S ACROSS MODELS TOGETHER ----------
## Plot base and climate RMSE's together
obs = N3 # study area here

# base dataframe of rmse's
base <- BASE %>% normalized_rmse(.,obs) %>% as.data.frame() %>% 
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
clim <- CLIM %>% normalized_rmse(.,obs) %>% as.data.frame() %>% 
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
topo <- TOPO %>% normalized_rmse(.,obs) %>% as.data.frame() %>% 
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
topoclim <- TOPOCLIM %>% normalized_rmse(.,obs) %>% as.data.frame() %>% 
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
  mutate(year_int = as.factor(year_int))

## plot

# colors
group.cols <- c("#F8766D","#00BFC4","#C77Cff","#7CAE00")
#group.cols <- c("#C77Cff","#7CAE00") # for plotting topos only

(rmse.plot <- dat %>% 
    #filter(model=='topo'|model=='topoclim') %>% # for plotting topos only
    ggplot(aes(x = year_int, y = rmse, col = model), col = cols) + 
    geom_boxplot(outlier.shape = NA) + 
    labs(y = 'NRMSE', x = "5 year chunks out") +
    scale_y_continuous(limits = c(0, 0.6)) +
    scale_color_manual(values=group.cols) +
    theme_bw()) # NA's are the shorter length of iterations for 2 of the models

#saveRDS(rmse.plot, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/rmse_5yr_chunks_5y_avg_init_reletavized.rds")

# saveRDS(rmse.plot, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_1_near_predictions/rmse_5yr_chunks_5y_avg_init_reletavized_rm_NaNs.rds")

# saveRDS(rmse.plot, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/rmse_5yr_chunks_5y_avg_init_reletavized_rm_NaNs.rds")

# saveRDS(rmse.plot, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/rmse_5yr_chunks_5y_avg_init_reletavized_rm_negatives_topos_only.rds")

#ggsave("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/rmse_5yr_chunks_5y_avg_init_reletavized_rm_z.png", plot = last_plot(), width = 4, height = 4, dpi = 400)

## **** Create table of all rmse values
rmse.35y <- dat %>% 
  #group_by(model,year_int) %>%
  group_by(model) %>%
  summarise(median = median(rmse, na.rm = T), lower_ci = quantile(rmse, na.rm = T, 0.05), upper_ci = quantile(rmse, na.rm = T, 0.95))

## save correct study areas to file

# write.csv(rmse.35y,"R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/rmses_for_model_studyarea_comparisons/in_sample_35y_nrmse_avg_allyears.csv")

# write.csv(rmse.35y,"R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/rmses_for_model_studyarea_comparisons/in_sample_35y_nrmse_yearchunks.csv")

# ------------------ Predicted median and credible intervals ------------
# # N
# base.pred = for.base.N$predOut
# topo.pred = for.topo.N$predOut
# clim.pred = for.clim.N$predOut
# topoclim.pred = for.topoclim.N$predOut

# # N2
# base.pred = for.base.N2$predOut
# topo.pred = for.topo.N2$predOut
# clim.pred = for.clim.N2$predOut
# topoclim.pred = for.topoclim.N2$predOut

# N3
base.pred = for.base.N3$predOut
topo.pred = for.topo.N3$predOut
clim.pred = for.clim.N3$predOut
topoclim.pred = for.topoclim.N3$predOut

# ## save median predictions and calculate crediable intervals
med.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = median)
low.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
up.pred.base <- apply(base.pred, MARGIN = c(1,2), FUN = quantile, 0.95)

# # save mediam pixel-year predicted cover
# write.csv(med.pred.base, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_1_near_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")

# # save lower CI pixel-year predicted cover
# write.csv(low.pred.base, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_90_credible_lower_5y_avg_init.csv")
# 
# # save upper CI pixel-year predicted cover
# write.csv(up.pred.base, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_90_credible_upper_5y_avg_init.csv")

## LOAD predictions for the base model

# # base median
# med.pred.base <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_median_5y_avg_init.csv")
# 
# # base lower CI
# low.pred.base <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_90_credible_lower_5y_avg_init.csv")
# 
# # base upper CI
# up.pred.base <- read.csv("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/median_predictions_for_mapping/in_sample_35y_base_90_credible_upper_5y_avg_init.csv")

med.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = median)
low.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
up.pred.clim <- apply(clim.pred, MARGIN = c(1,2), FUN = quantile, 0.95)

# med.pred.topo <- apply(topo.pred, MARGIN = c(1,2), FUN = median)
# low.pred.topo <- apply(topo.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
# up.pred.topo <- apply(topo.pred, MARGIN = c(1,2), FUN = quantile, 0.95)
# 
# med.pred.topoclim <- apply(topoclim.pred, MARGIN = c(1,2), FUN = median)
# low.pred.topoclim <- apply(topoclim.pred, MARGIN = c(1,2), FUN = quantile, 0.05) # low
# up.pred.topoclim <- apply(topoclim.pred, MARGIN = c(1,2), FUN = quantile, 0.95)

# ------------------ Predicted vs. Observed plots -----------------------

# ** these plots had an error and are not currently used to generate any results **

# pred vs. observed change in cover over entire forecast period
obs = N

## use sweep to calc pred vs. latent (changes from 1 to accumulated 31)
plot(sweep(med.lat.base[2:31,], 2, med.lat.base[1,]),sweep(med.pred.base[2:31,], 2, med.pred.base[1,]), main = "BASE MODEL", xlab = "Latent CHANGE IN COVER (%)", ylab= "PREDICTED CHANGES IN COVER (%)", col = "aquamarine3", pch = 20, cex.lab = 1.4)
abline(0,1, col = "purple", lwd = 2.5)
# change in cover for single time step, 1 to 31
plot(med.lat.base[31,]-med.lat.base[1,],med.pred.base[31,]-med.pred.base[1,], main = "BASE MODEL", xlab = "Latent CHANGE IN COVER (%)", ylab= "PREDICTED CHANGES IN COVER (%)", col = "aquamarine3", pch = 20, cex.lab = 1.4)
abline(0,1, col = "purple", lwd = 2.5)

# #base model ## these need fixed!!
# plot(obs[2:36,]-obs[1,],med.pred.base[2:36,]-med.pred.base[1,], main = "BASE MODEL", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3", pch = 20, cex.lab = 1.4)
# abline(0,1, col = "purple", lwd = 2)
# # topo
# plot(obs[2:36,]-obs[1,],med.pred.topo[2:36,]-med.pred.topo[1,], main = "pred. vs. obs. change over forecast horizon", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3")
# abline(0,1, col = "purple", lwd = 2)
# #clim model
# plot(obs[2:36,]-obs[1,],med.pred.clim[2:36,]-med.pred.clim[1,], main = "pred. vs. obs. change over forecast horizon", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3")
# abline(0,1, col = "purple", lwd = 2)
# #topoclim model
# plot(obs[2:36,]-obs[1,],med.pred.topoclim[2:36,]-med.pred.topoclim[1,], main = "pred. vs. obs. change over forecast horizon", xlab = "obs. changes in cover (%)", ylab= "pred. changes in cover (%)", col = "aquamarine3")
# abline(0,1, col = "purple", lwd = 2)

# ----------------- Plot predictions over time ------------------
## ggplot version
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
    
    # scale_y_continuous(limits = c(-5, 45)) + # good scale for N2
    scale_y_continuous(limits = c(-5, 65)) + # good scale for N3
    # scale_y_continuous(limits = c(-5, 35)) + # good scale for N 
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

# test plot function for a given pixel
plot_pix_gg(711, N2)

## plot low, medium, and high cover all-together

# N
(three.panel <- plot_grid(plot_pix_gg(707, N) + labs(x = '', y = ''), 
          plot_pix_gg(152, N) + theme(legend.position = "none") + labs(x = '', y = ''), 
          plot_pix_gg(120, N) + theme(legend.position = "none") + labs(x = '', y = ''), ncol = 3))

# N=159 or 152, increasing trend pixel, N = 703 increasing trend pixel with high var, N = 707, low cover, increasing.

# N2
(three.panel <- plot_grid(plot_pix_gg(711, N2) + theme(legend.position = "none") + labs(x = ''), 
                          plot_pix_gg(22, N2) + theme(legend.position = "none") + labs(x = '', y = ''), 
                          plot_pix_gg(2050, N2) + theme(legend.position = "none") + labs(y = '', x = ''), ncol = 3))

# N2=22,100 low cover, modest increase; N2 = 160 moderate cover tiny increase; 
# N2 =220 moderate cover, no increase; N2 = 425, moderate, increasing

# N3
(three.panel <- plot_grid(plot_pix_gg(483, N3) + theme(legend.position = "none") + labs(x = '', y = ''), 
          plot_pix_gg(485, N3) + theme(legend.position = "none") + labs(y = ''), 
          plot_pix_gg(488, N3) + theme(legend.position = "none") + labs(x = '', y = ''), ncol = 3))

# N3 = 160 or 215 very low cover, minimal increase; N3 = 1035 low-moderate cover increasing; 
# N3 = 707 moderate cover, small increasing trend; N3 = 700 moderate cover, increasing cover trend;
# N3 = 332,333,336 nearly the highest cover, increasing at first, then levels off
# N3 = 219, 276 climate model did really well
# see N3 = 427-433 for a set of pixels in the same 'row' spanning the 'range margin'
# Also see N = 483,485,488 for a set of pixels in the same 'row' spanning the 'range margin'
# N3 = 2313, both models under-predict cover (climate usually over predicts)

# ggsave("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/single_pixel_low_mid_high_cover_5yr_avg_init.png", plot = last_plot(), dpi = 400)

# save as rds to put into a plot later
saveRDS(three.panel, "R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_OOS_2_far_predictions/single_pixel_low_mid_high_cover_5yr_avg_init.rds")
# ------- PLOT OBSERVED COVER OVER TIME FOR ANY PIXEL -----
# test out different pixels to plot based on 36-yr trends

(obs.plot <- N %>% as.data.frame() %>% rownames_to_column("year") %>% 
   pivot_longer(!year, values_to = "obs", names_to = "pixel") %>% 
   filter(pixel == "1530") %>% # any pixel between 1-2352
   mutate(year = as.numeric(year)+1985) %>% 
   ggplot(aes(x = year, y = obs)) + 
   geom_point(aes(x = year, y = obs)) + 
   geom_line(aes(x = year, y = obs), lty = 2) + 
   scale_y_continuous(limits = c(0, 20)) +
   labs(y = "% TREE COVER", x = "YEAR") +
   theme_bw())

ggsave("R:/Shriver_Lab/PJspread/figures/observed_margin_pix_1530.png", plot = last_plot(), width = 5, height = 4, dpi = 300)

#1506
#1530
#2226
# ------- PLOT LATENT COVER ON OBSERVED COVER OVER TIME ------
# test out ploting observed cover and latent cover (latent cover for base model)
med.lat <- apply(mod1$Nlat, MARGIN = c(1,2), FUN = median)
up.lat <- apply(mod1$Nlat, MARGIN = c(1,2), FUN = quantile, 0.95)
low.lat <- apply(mod1$Nlat, MARGIN = c(1,2), FUN = quantile, 0.05)

#

plot_lat_gg <- function(pix) {
  
  px <- as.character(pix)
  pxv <- paste("V",as.character(pix), sep = "")
  
  ## re-organize latent cover data
  
  #median values
  lat.med <- med.lat %>% as.data.frame() %>% 
    rownames_to_column("year") %>%
    pivot_longer(!year, values_to = "lat.med", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% 
    dplyr::select(!pixel) %>%
    mutate(year = as.numeric(year))
  
  # lat lower ci values
  lat.low <- low.lat %>% as.data.frame() %>% 
    rownames_to_column("year") %>%
    pivot_longer(!year, values_to = "lat.low", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% 
    dplyr::select(!pixel) %>%
    mutate(year = as.numeric(year))
  
  # lat upper ci values
  lat.up <- up.lat %>% as.data.frame() %>% 
    rownames_to_column("year") %>%
    pivot_longer(!year, values_to = "lat.up", names_to = "pixel") %>% 
    filter(pixel == pxv) %>% 
    dplyr::select(!pixel) %>%
    mutate(year = as.numeric(year))
  
  ## re-organize predicted cover data
  
  # # median predicted values
  # med.pred <- med.pred.clim %>% as.data.frame() %>% 
  #   rownames_to_column("year") %>%
  #   pivot_longer(!year, values_to = "med.pred", names_to = "pixel") %>% 
  #   filter(pixel == pxv) %>% 
  #   dplyr::select(!pixel) %>%
  #   mutate(year = as.numeric(year)) %>%
  #   filter(year !=1 )  # remove first year (true value, not part of the forecast)
  # 
  # # low ci prediction
  # low.pred <- low.pred.clim %>% as.data.frame() %>% 
  #   rownames_to_column("year") %>%
  #   pivot_longer(!year, values_to = "low.pred", names_to = "pixel") %>% 
  #   filter(pixel == pxv) %>% 
  #   dplyr::select(!pixel) %>%
  #   mutate(year = as.numeric(year)) %>%
  #   filter(year !=1 )  # remove first year (true value, not part of the forecast)
  # 
  # up.pred <- up.pred.clim %>% as.data.frame() %>% 
  #   rownames_to_column("year") %>%
  #   pivot_longer(!year, values_to = "up.pred", names_to = "pixel") %>% 
  #   filter(pixel == pxv) %>% 
  #   dplyr::select(!pixel) %>%
  #   mutate(year = as.numeric(year)) %>%
  #   filter(year !=1 )  # remove first year (true value, not part of the forecast)
  # 
  # observed 
  Ndat <- N %>% as.data.frame() %>% 
    rownames_to_column("year") %>% 
    pivot_longer(!year, values_to = "obs", names_to = "pixel") %>% 
    filter(pixel == px) %>% 
    mutate(year = as.numeric(year)) %>% 
    filter(year <= 31) %>%
    dplyr::select(!pixel)
  
  # combine predictions and observations into one dataframe
  plot.dat <- left_join(Ndat, lat.med) %>% 
    left_join(.,lat.low) %>% left_join(.,lat.up) %>% 
    #left_join(.,med.pred) %>% left_join(.,low.pred) %>% left_join(.,up.pred) %>% # not plotting predicted data
    mutate(year = year+1985)
  
  # plot
  (test <- plot.dat %>% 
      ggplot(aes(x = year, y = obs)) + 
      geom_point(aes(x = year, y = obs)) + 
      geom_line(aes(x = year, y = obs), lty = 2) + 
      geom_ribbon(aes(ymin = lat.low, ymax = lat.up), 
                  alpha=0.4, fill = "darkgrey") + 
      geom_line(aes(x = year, y = lat.med), col = "darkgrey", lwd = 1.1) +
      #geom_ribbon(aes(ymin = low.pred, ymax = up.pred), alpha=0.3, fill = "#00BFC4") +
      #geom_line(aes(x = year, y = med.pred), col = "#00BFC4") + 
      scale_y_continuous(limits = c(-2, 25)) + 
      labs(x = "YEAR", y = "TREE COVER (%)") + 
      theme(legend.position="none", text = element_text(size=5)) +
      theme_classic())
}

plot_lat_gg(159)


ggsave("R:/Shriver_Lab/PJspread/evaluate_out_of_sample/35y_insample_predictions/single_pixel_555_base_latent_w_obs.png", plot = last_plot(), width = 6, height = 4, dpi = 400)

# Good representative pixels
# N
# 171
# N2
# 608 - low cover, increasing