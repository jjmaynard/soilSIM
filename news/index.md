# Changelog

## soilSIM 0.1.0

Initial release.

- Data-source-agnostic simulation/fusion core: percentile-triplet and
  arbitrary-percentile distribution fitting
  ([`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md),
  [`simulate_from_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_from_percentiles.md)),
  correlated Monte Carlo simulation with compositional (ILR) texture
  handling
  ([`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)),
  Gaussian-process depth-trend modeling
  ([`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md),
  [`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)),
  and Bayesian updating/fusion in both scalar
  ([`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
  and family,
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md))
  and raster-native
  ([`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md),
  [`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md))
  forms.
- SSURGO data-source adapter: tabular acquisition
  ([`download_and_prepare_ssurgo()`](https://jjmaynard.github.io/soilSIM/reference/download_and_prepare_ssurgo.md)),
  cleaning
  ([`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)),
  infilling
  ([`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md)),
  and simulation-ready raster/percentile extraction
  ([`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)).
- SOLUS100 data-source adapter: raster percentile fetch
  ([`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md))
  supplying the likelihood side of the raster fusion pipeline.
- Profile, component, and depth simulation: component composition
  ([`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md),
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md))
  and horizon depth/thickness simulation
  ([`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md),
  [`simulate_profile_depths_by_mukey()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths_by_mukey.md)).
- Available water storage modeling via Van Genuchten / ROSETTA
  ([`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md),
  [`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md),
  [`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)).
- Statistics and workflow validation/diagnostics
  ([`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md),
  [`validate_complete_workflow()`](https://jjmaynard.github.io/soilSIM/reference/validate_complete_workflow.md),
  [`generate_validation_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_validation_report.md)).
- Four vignettes covering the tabular Monte Carlo pipeline,
  profile/depth simulation, available water storage, and multi-source
  raster fusion (SSURGO x SOLUS100), all built against real SSURGO data
  for a Sierra Nevada foothills and Salinas Valley area of interest.
