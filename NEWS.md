# soilSIM 0.1.0.9000 (development version)

* Multi-property raster fusion. `simulate_cokey_generalized()`, `simulate_ssurgo_mapunit_draws()`,
  `fetch_ssurgo_percentiles()`, and `extract_mukey_joint_ensemble()` gain a `requested_properties`
  argument that restricts the SSURGO Monte Carlo simulation to the properties actually needed - the
  per-cokey depth-trend GP fit is the pipeline's dominant cost and scales with the property count.
* `run_stage1_fusion()` now simulates **only** its own property (plus the full sand/silt/clay draw
  for texture members) by default. The prior is statistically equivalent to - not bit-identical to -
  the previous all-property simulation (the pipeline was already unseeded), and a single-property
  map can show marginally more `NA` where a map unit's components carry no data for that one
  property. Restriction is honoured only under the default `vertical_correlation_method =
  "joint_copula"`; `"gp_quantile_retrofit"` falls back to a full simulation with a warning.
* New `run_stage1_fusion_multi()`: fuses many properties over an AOI (and multiple depth windows)
  from a **single** SSURGO simulation, seeding the same per-property disk caches
  `run_stage1_fusion()` reads. Replaces an N-properties x M-windows loop of `run_stage1_fusion()`
  calls (each of which re-ran the full simulation); every leaf now shares one draw set, so their
  `NA` masks are mutually consistent.
* New optional `seed` argument on `simulate_ssurgo_mapunit_draws()`, `fetch_ssurgo_percentiles()`,
  `extract_mukey_joint_ensemble()`, `run_stage1_fusion()`, `run_stage1_fusion_group()`, and
  `run_stage1_fusion_multi()` for opt-in reproducibility (default `NULL` = unchanged stochastic
  behavior). Seeds the sequential RNG stream and the parallel depth-trend `future.seed`;
  determinism is conditional on the live SSURGO/SOLUS data being unchanged.

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
