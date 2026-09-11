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
(`adapter-ssurgo-acquire.R`, `adapter-ssurgo-process.R`,
`adapter-ssurgo-simulate.R`) and SOLUS100 (`adapter-solus.R`) are the
two adapters shipped today; new data sources are intended to plug in as
additional adapters alongside them, reusing the same core rather than
duplicating it.

## Data acquisition & processing

The SSURGO adapter:
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)/
[`fetch_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_data.md)
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
[`diagnose_workflow()`](https://jjmaynard.github.io/soilSIM/reference/diagnose_workflow.md)/[`generate_validation_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_validation_report.md)
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

[`simulate_monte_carlo()`](https://jjmaynard.github.io/soilSIM/reference/simulate_monte_carlo.md)
is the master pipeline for correlated Monte Carlo simulation of soil
properties from percentile inputs;
[`simulate_correlated_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_properties.md)
is its Cholesky-copula core simulator.

## GP depth modeling & multivariate adjustment (core)

[`fit_depth_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_depth_gp_models.md)/[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
fit and apply Gaussian-process depth-trend models;
[`apply_depth_gp_to_simulation()`](https://jjmaynard.github.io/soilSIM/reference/apply_depth_gp_to_simulation.md)
combines Monte Carlo output with GP depth trends across a whole dataset
while preserving cross-property correlation structure.

## Profile, component & depth simulation

[`simulate_component_composition()`](https://jjmaynard.github.io/soilSIM/reference/simulate_component_composition.md)/[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
simulate component composition and per-cokey properties;
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)/[`simulate_profile_depths()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths.md)
simulate horizon depths and thicknesses for whole soil profiles.

## AWS / Van Genuchten modeling

[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)/[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)/[`compute_aws()`](https://jjmaynard.github.io/soilSIM/reference/compute_aws.md)
estimate available water storage from ROSETTA pedotransfer parameters.

## Bayesian updating (scalar, core)

[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md)
and family provide closed-form same-family fusion of a prior and
observed data;
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
provides a fully general grid-KDE fusion for arbitrary distributions;
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)
fuses compositional texture data jointly in ILR space.

## Multi-source raster fusion pipeline

[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)/[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
(in `core-fusion.R`) are the generic raster-native fusion core,
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
