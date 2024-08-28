# Pinyon-Juniper woodland range expansion model

### Modeling and figure generating scripts for *Intrinsic population processes alone can accurately forecast range shifts in a dominant dryland tree species* by Elise Pletcher and Robert K. Shriver 

## Folders:
* **model_data_prepping**: a script for taking tree cover and environmental covariates for the in sample study area and reformatting them as inputs for the mcmc_samplers.
* **mcmc_samplers**: contains multiple scripts. Uses reformatted data generated from model_data_prepping to fit spatiotemporal models. The output are workspace files with fit parameter estimates.
* **generate_forecasts**: contains multiple scripts, all for taking parameter estimates and forecasting changes in tree cover for a 5 year interval in sample, and 35 y for all 3 study landscapes.
* **additional_figures**: contains scripts for mapping forecasts in raster format and plotting supplemental figures.
