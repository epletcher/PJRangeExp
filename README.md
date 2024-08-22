# Pinyon-Juniper woodland range expansion model

## There are four folders:
* model_data_prepping: a script for taking tree cover and environmental covariates for the in sample study area and reformatting them as inputs for the mcmc_samplers.
* mcmc_samplers: contains multiple scripts. Uses reformatted data generated from model_data_prepping to fit spatiotemporal models. The output are workspace files with fit parameter estimates.
* generate_forecasts: contains multiple scripts, all for taking parameter estimates and forecasting changes in tree cover for a 5 year interval in sample, and 35 y for all 3 study landscapes.
* additional_figures: contains scripts for mapping forecasts in raster format and plotting supplemental figures.

### Each Sampler*.R file is a different version of an MCMC sampler written to model Pinyon-Juniper spread across a landscape in Central OR. Samplers vary by the complexity of the process model. All model versions use the same observation and process error.

Sampler.R - Basic spatiotemporal model with no covariates

Sampler_Topo.R - Basic model + two spatially variable topographic covariates: heatload and elevation 

Sampler_climate.R - Basic model + two temporally variable climatic covariates: temp. and precip.

Sampler_TopoClim.R - Basic model + climatic covariates + topographic covariates

Sampler_quad.R - Basic model, with an added quadratic term to allow for negative and postive density dependence

### SamplerFunctions.R file contains functions for updating each version of the process model and for sampling posterior distributions of error terms and the latent state. As well as functions for 
