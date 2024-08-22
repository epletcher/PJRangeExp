# MCMC Samplers
### Each Sampler*.R file is a different version of an MCMC sampler written to model Pinyon-Juniper spread across a landscape in Central OR. Samplers vary by the complexity of the process model. All model versions use the same observation and process error.

## Contents:
* **Sampler.R** - Basic spatiotemporal model with no covariates
* **Sampler_Topo.R** - Basic model + two spatially variable topographic covariates: heatload and elevation 
* **Sampler_climate.R** - Basic model + two temporally variable climatic covariates: temp. and precip.
* **Sampler_TopoClim.R** - Basic model + climatic covariates + topographic covariates
* **SamplerFunctions.R** - file contains functions for updating each version of the process model and for sampling posterior  distributions of error terms and the latent state.
