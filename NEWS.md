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
* `simulate_ssurgo_mapunit_draws()`'s raw SSURGO tabular-download cache is now keyed by AOI alone
  instead of AOI + depth window (`download_ssurgo_tabular()` never depth-filters, so the cached
  content was already depth-window-independent) - one AOI's tabular download is now shared across
  every depth window/call requested for it, including between `run_stage1_fusion_multi()`'s
  wide-span simulation and a later single-window `run_stage1_fusion()` call.
* **Correctness fix**: `simulate_cokey_generalized()`'s simulated `soc` property is now a genuine
  soil-organic-carbon estimate (`om * OM_TO_SOC_FACTOR`, the Van Bemmelen factor) instead of raw
  SSURGO organic matter passed through unconverted. Previously every `"soc"` fusion result (via
  `run_stage1_fusion()`/`run_stage1_fusion_multi()` with `solus_variable = "soc"`) fused SSURGO
  organic matter % (prior) against SOLUS100's real organic carbon % (likelihood) - two different
  quantities differing by roughly a factor of 1.724. `"soc"` percentile/posterior values from this
  release will differ from earlier releases; the change should be read as a bug fix, not a
  regression. No change was needed to the KSSL reference correlation matrix (correlation is
  invariant under a linear rescale of one variable).
* New `fetch_solus_low_pred_high_multi()` / `fetch_solus_percentiles_multi()`: fetch every needed
  SOLUS100 variable's low/prediction/high rasters for one depth window in a **single**
  `soilDB::fetchSOLUS()` call (`variables`, `depth_slices`, and `output_type` are all requested as
  vectors at once), instead of one call per `(variable, output_type)`. `run_stage1_fusion_multi()`
  now fetches SOLUS via this batched path - one call per window covering every standalone
  property and texture-group member, down from up to `3 * n_variables` calls per window - falling
  back to the existing per-variable scalar fetch only if the batched request itself errors. The
  scalar `fetch_solus_low_pred_high()` / `fetch_solus_percentiles()` (used by single-property
  `run_stage1_fusion()`) are unchanged.
* SSURGO chemistry-property coverage extended to 5 new properties: `caco3`, `ec`, `ecec`,
  `gypsum`, `sar` (SSURGO chorizon column names live-confirmed against a real gSSURGO query).
  `download_ssurgo_tabular()`'s default `properties`, `create_ssurgo_property_lookup_working()`,
  `simulate_cokey_generalized()`'s `param_order`, `property_to_sim_column()`,
  `normalize_requested_properties()`, `SSURGO_SIM_PROPERTY_COLUMNS`, and `infill_soil_data()` all
  now recognize them, so they can be fetched via `fetch_ssurgo_percentiles()` and fused via
  `run_stage1_fusion()`/`run_stage1_fusion_multi()` like any of the original 9 properties. They
  have no entry in the static KSSL reference correlation matrices (fit years before these
  properties existed in this pipeline) - `simulate_ssurgo_mapunit_draws()` now builds its
  per-genhz correlation matrices via `build_kssl_fallback_matrix()` (sized to the full property
  vocabulary), so these 5 simulate **uncorrelated** with everything else (and each other) rather
  than erroring; real correlations can be added later via a KSSL data re-fit with no code changes.
  `simulate_cokey_generalized()` also gained a defensive fallback
  (`extend_corr_matrix_with_identity()`) so a direct caller supplying an older correlation-matrix
  list that predates this change degrades to the same independence default instead of erroring.
  SSURGO coverage for these 5 properties is real but sparser than the original 9 in spot-checks
  (e.g. `ecec` was populated in only ~25% of sampled horizons for one AOI) - expect more `NA` in
  their fused output.
* New optional `restriction_depth` argument on `remarginalized_awc()`: bedrock/restriction-depth-
  aware AWC truncation. SOLUS100 predicts every property at every depth even for shallow soils, so
  without this, a pixel with bedrock at 40 cm was credited with water storage in depth windows
  entirely or partially below bedrock - a real, physically-wrong AWC inflation for shallow soils.
  `NULL` (default) preserves prior output exactly; when supplied (typically
  `fetch_solus_restriction_depth(aoi_vect)`, new below), each window's thickness contribution
  becomes the per-pixel effective thickness `pmax(0, pmin(bottom, restriction_depth) - top)` - a
  window straddling the restriction gets partial credit, one entirely below it contributes exactly
  0 (not `NA`, which would incorrectly blank the whole pixel's AWC rather than just that window).
* New `fetch_solus_site_level()` / `fetch_solus_restriction_depth()`: fetch a SOLUS100 site-level
  (depth-independent) variable - `anylithicdpt` (depth to bedrock) / `resdept` (depth to
  restriction) - in one `fetchSOLUS(depth_slices = "all", ...)` call.
  `fetch_solus_restriction_depth()` defaults to `anylithicdpt`, confirmed (per SOLUS100's own
  published training methodology) to be trained only on lithic/paralithic bedrock restriction
  records - hard-bedrock-specific by construction, needing no separate SSURGO-side classifier;
  `resdept` (trained on every restriction record, any kind) is offered as a broader, more
  conservative non-default alternative. Both are right-censored at 201 cm in SOLUS100's training
  data (new `SOLUS_RESTRICTION_CENSOR_CM` constant) - values at/above that are recoded to `Inf`
  (no restriction detected), not treated as a literal depth.

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
