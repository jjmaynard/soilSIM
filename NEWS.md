# soilSIM 0.1.0

Initial release.

* Data-source-agnostic simulation/fusion core: percentile-triplet and arbitrary-percentile
  distribution fitting (`fit_percentile_triplet()`, `simulate_from_percentiles()`), correlated
  Monte Carlo simulation with compositional (ILR) texture handling
  (`generate_monte_carlo_realizations()`), Gaussian-process depth-trend modeling
  (`build_stratified_gp_models()`, `integrate_monte_carlo_with_gp()`), and Bayesian
  updating/fusion in both scalar (`bayes_update_normal_normal()` and family,
  `bayesian_update()`) and raster-native (`fuse_property_adaptive()`, `fuse_texture_group()`)
  forms.
* SSURGO data-source adapter: tabular acquisition (`download_and_prepare_ssurgo()`), cleaning
  (`process_ssurgo_data()`), infilling (`process_soil_properties_comprehensive()`), and
  simulation-ready raster/percentile extraction (`simulate_ssurgo_mapunit_draws()`).
* SOLUS100 data-source adapter: raster percentile fetch (`fetch_solus_percentiles()`) supplying
  the likelihood side of the raster fusion pipeline.
* Profile, component, and depth simulation: component composition
  (`sim_component_comp()`, `simulate_cokey_generalized()`) and horizon depth/thickness
  simulation (`simulate_and_perturb_soil_profiles()`, `simulate_profile_depths_by_mukey()`).
* Available water storage modeling via Van Genuchten / ROSETTA
  (`van_genuchten()`, `simulate_vg_aws()`, `calculate_aws_df()`).
* Statistics and workflow validation/diagnostics (`analyze_soil_statistics()`,
  `validate_complete_workflow()`, `generate_validation_report()`).
* Four vignettes covering the tabular Monte Carlo pipeline, profile/depth simulation, available
  water storage, and multi-source raster fusion (SSURGO x SOLUS100), all built against real SSURGO
  data for a Sierra Nevada foothills and Salinas Valley area of interest.
