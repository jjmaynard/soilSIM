# soilSIM: Multi-Source Soil Property Simulation and Bayesian Data Fusion

A data-source-agnostic framework for simulating soil properties and
fusing estimates from multiple soil data products. Provides
percentile-triplet and arbitrary-percentile distribution fitting,
correlated Monte Carlo simulation (including compositional texture
handling via isometric log-ratio transforms), Gaussian-process
depth-trend modeling, descriptive statistics/correlation analysis, and
Bayesian updating/fusion in both scalar and raster-native forms. Ships
with SSURGO and SOLUS100 data-source adapters wired into the fusion
pipeline as prior and likelihood sources respectively; the simulation
and fusion core is designed to accommodate additional gridded soil data
sources (e.g. HWSD, SoilGrids) as further adapters.

## Overview

The package is split into a generic simulation/fusion **core** -
percentile fitting, Monte Carlo simulation, GP depth modeling, and
Bayesian fusion - that operates on any low/representative/high triplet
or prior/likelihood pair regardless of where the numbers came from, and
a set of source-specific **adapters** that acquire, process, and shape a
given soil data product's output into that generic form. SSURGO
(`ssurgo-acquisition.R`, `ssurgo-processing.R`, `ssurgo-simulation.R`)
and SOLUS100 (`solus-simulation.R`) are the two adapters shipped today;
new data sources are intended to plug in as additional adapters
alongside them, reusing the same core rather than duplicating it.

## Data acquisition & processing

The SSURGO adapter:
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)/
[`download_and_prepare_ssurgo()`](https://jjmaynard.github.io/soilSIM/reference/download_and_prepare_ssurgo.md)
acquire and cache tabular SSURGO data for an area of interest;
[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)
cleans and standardizes it into modeling-ready horizon/component tables;
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md)/[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)
fill missing property values using pedologically informed recovery
strategies.

## Statistics & diagnostics

[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)
characterizes processed data (correlations, distribution fitting,
outliers, summary statistics);
[`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md)/[`generate_validation_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_validation_report.md)
provide end-to-end quality assurance across the full simulation
workflow.

## Distribution fitting (core)

[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md)/[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md)
fit and sample closed-form distributions (Normal, Beta, metalog,
triangular, linear-CDF) from a low/representative/high triplet;
[`simulate_from_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_from_percentiles.md)
generalizes this to an arbitrary number of percentiles via
linear-CDF/spline/KDE methods;
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)
provide the isometric log-ratio transforms used for compositional
texture data.

## Monte Carlo simulation (core)

[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)
is the master pipeline for correlated Monte Carlo simulation of soil
properties from percentile inputs;
[`simulate_correlated_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_properties.md)
is its Cholesky-copula core simulator.

## GP depth modeling & multivariate adjustment (core)

[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md)/[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
fit and apply Gaussian-process depth-trend models;
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)
combines Monte Carlo output with GP depth trends across a whole dataset
while preserving cross-property correlation structure.

## Profile, component & depth simulation

[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)/[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
simulate component composition and per-cokey properties;
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)/[`simulate_profile_depths_by_mukey()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths_by_mukey.md)
simulate horizon depths and thicknesses for whole soil profiles.

## AWS / Van Genuchten modeling

[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)/[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)/[`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)
estimate available water storage from ROSETTA pedotransfer parameters.

## Bayesian updating (scalar, core)

[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
and family provide closed-form same-family fusion of a prior and
observed data;
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
provides a fully general grid-KDE fusion for arbitrary distributions;
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)
fuses compositional texture data jointly in ILR space.

## Multi-source raster fusion pipeline

[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)/[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
(in `raster-fusion.R`) are the generic raster-native fusion core,
combining prior and likelihood `SpatRaster` inputs cell-by-cell over an
area of interest; the SSURGO adapter
([`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md))
supplies the prior side and the SOLUS100 adapter
([`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md))
supplies the likelihood side, with
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)
providing disk caching for both. Designed so additional data-source
adapters can supply either side without changes to the fusion core.

## Utilities

[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)/[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)/[`read_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/read_soil_data.md)/
[`write_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/write_soil_data.md)
and a shared logging/configuration framework used throughout the
package.

## See also

Useful links:

- <https://jjmaynard.github.io/soilSIM/>

- <https://github.com/jjmaynard/soilSIM>

- Report bugs at <https://github.com/jjmaynard/soilSIM/issues>

## Author

**Maintainer**: Jonathan Maynard <jjmaynard9@gmail.com>
