library(tidyverse)
library(coda)

# function to retrive parameter estimates for any model
retrieve_param <- function(model) {
  if(model=="topo") {model="topo_"}
  if(model=="base") {model="base_v"}
  if(model=="clim") {model="r_clim"}
  if(model %in% c("base_v","topo_","topoClim","r_clim","quad")) {
    
    files <- list.files("R:/Shriver_Lab/PJspread/sampleroutput", pattern=model, include.dirs = TRUE)
    
    mod <- list()
    
    for(c in 1:length(files)) {
      
      load(files[c])
    
      if(is.matrix(alphaOut)) {
        a_names <- dim(alphaOut)[2] 
        for(i in 1:dim(alphaOut)[2]) {
          a_names[i] <- paste("alpha",as.character(i-1),sep='')
        }
        colnames(alphaOut) <- a_names
      }
      
      if(is.matrix(betaOut)) {
        b_names <- dim(betaOut)[2]
        for(i in 1:dim(betaOut)[2]) {
          b_names[i] <- paste("beta",as.character(i-1),sep='')
        }
        colnames(betaOut) <- b_names
      }
      
      colnames(tauOut) <- "tau"
      colnames(sig.pOut) <- "sigma.p"
      
      chain <- cbind(alphaOut,betaOut,tauOut,sig.pOut)
      
      thin.iter <- seq(10,dim(chain)[1], by = 10) # every 10th iter
      thin.chain <- chain[thin.iter,] # thin chain
      
      mod[[c]] <- thin.chain
    }
    
    return(mod)
  }
  
}

mod1 <- retrieve_param("base") # model type here
mod2 <- retrieve_param("topo")
mod3 <- retrieve_param("clim")
mod4 <- retrieve_param("topoClim")

# ------ ggplot traceplot ----------

gg_traceplot <- function(mod) { # model object needs to be a list of matrices, w/ each matrix of [iter,param]
 
  iteration <- seq(1:dim(mod[[1]])[1])
  chain <- as.factor(rep(1, dim(mod[[1]])[1]))
  mod.df <- mod[[1]] %>% as.data.frame() %>% cbind(.,iteration) %>% cbind(.,chain)
  
  for(i in 2:length(mod)) {
    
    chain <- as.factor(rep(i, dim(mod[[i]])[1]))
    mod.df <- mod[[i]] %>% as.data.frame() %>% cbind(.,iteration) %>% cbind(.,chain) %>% rbind(mod.df)
    
  }
  
  mod.df %>% 
    pivot_longer(cols = -c(iteration,chain), names_to = "parameter", values_to = "value") %>% 
    filter(iteration > max(iteration)/2) %>% # filter mod.df to burnin = half of length of each chain
    ggplot(aes(x = iteration, y = value)) + 
    geom_line(aes(col = chain)) + 
    facet_wrap(~parameter, scales="free_y") + 
    theme_bw()
  
}

gg_traceplot(mod1)
ggsave("R:/Shriver_Lab/PJspread/figures/mcmc_base.png", dpi = 600)
gg_traceplot(mod2)
ggsave("R:/Shriver_Lab/PJspread/figures/mcmc_topo.png", dpi = 600)
gg_traceplot(mod3)
ggsave("R:/Shriver_Lab/PJspread/figures/mcmc_clim.png", dpi = 600)
gg_traceplot(mod4)
ggsave("R:/Shriver_Lab/PJspread/figures/mcmc_topoclim.png", dpi = 600)


# ###------- base traceplot --------
# traceplot <- function(mod) {
#   
#   if(dim(mod[[1]])[2]<=4) {par(mfrow = c(2,2))}
#   if(dim(mod[[1]])[2]>4 & dim(mod[[1]])[2]<10) {par(mfrow = c(3,3))}
#   if(dim(mod[[1]])[2]>=10) {par(mfrow = c(3,4))}
#   
#   for(i in 1:dim(mod[[1]])[2]) {
#     
#     plot(mod[[1]][,i], type = "l", ylab = colnames(mod[[1]])[i], xlab = "iteration", col = "purple")
#     lines(mod[[2]][,i], col = "aquamarine4")
#   }
#   
# }
# 
# traceplot(mod1)
# traceplot(mod2)
# traceplot(mod3)
# traceplot(mod4)
