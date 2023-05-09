# Pinyon-Juniper woodland spread model

### Each Sampler*.R file is a different version of an MCMC sampler written to model Pinyon-Juniper spread across a landscape in Central OR. Samplers vary by the complexity of the process model. All model versions use the same observation and process error.

Sampler.R - Basic spatiotemporal model with no covariates

Sampler_Topo.R - Basic model + two spatially variable topographic covariates: heatload and elevation 

Sampler_climate.R - Basic model + two temporally variable climatic covariates: temp. and precip.

Sampler_TopoClim.R - Basic model + climatic covariates + topographic covariates

Sampler_quad.R - Basic model, with an added quadratic term to allow for negative and postive density dependence

### SamplerFunctions.R file contains functions for updating each version of the process model and for sampling posterior distributions of error terms and the latent state. As well as functions for 
