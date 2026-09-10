# soilSIM 0.1.0.9000 (development version)

## Naming and package-structure re-architecture (in progress)

* **`R/` reorganized** into a generic core + data-source adapters layout
  (`core-*.R`, `adapter-ssurgo-*.R`, `adapter-solus.R`, `model-aws.R`, plus
  `geometry.R` / `logging.R` / `config.R` / `io.R` / `property-names.R` split out
  of the old `utils.R`, and `diagnostics.R` / `diagnostics-checks.R` out of
  `validation-diagnostics.R`). No user-visible change from the move itself.
* **S3 result objects.** `simulate_monte_carlo()`, `analyze_soil_statistics()` and
  `fit_depth_gp_models()` now return classed lists (`soilSIM_simulation`,
  `soilSIM_statistics`, `soilSIM_gp_models`) with `print()` / `summary()` methods.
  `$`-access is unchanged.
* **Tier-1 entry points renamed** (old names soft-deprecated, forwarding shims kept
  until 0.4.0):
  * `generate_monte_carlo_realizations()` -> `simulate_monte_carlo()`
  * `build_stratified_gp_models()` -> `fit_depth_gp_models()`
  * `adjust_multivariate_depthwise_GP()` -> `adjust_simulation_depthwise()`
  * `integrate_monte_carlo_with_gp()` -> `apply_depth_gp_to_simulation()`
  * `run_stage1_fusion()` / `_group()` / `_multi()` -> `run_fusion()` / `run_fusion_group()` / `run_fusion_multiproperty()`
  * `download_and_prepare_ssurgo()` -> `fetch_ssurgo_data()`
  * `infill_soil_data()` -> `infill_ssurgo_data()`
  * `sim_component_comp()` -> `simulate_component_composition()`
  * `calculate_aws_df()` -> `compute_aws()`
  * `validate_complete_workflow()` -> `diagnose_workflow()`
* **Systematic renames** (also soft-deprecated with shims):
  * Migration-scar suffixes dropped: `*_working` / `*_working_compatible` / `*_ssurgo_compatible`
    (e.g. `create_ssurgo_property_lookup_working()` -> `build_ssurgo_property_lookup()`,
    `process_horizon_data_working_compatible()` -> `process_ssurgo_horizons()`).
  * Config accessors: `get_default_configuration()` -> `default_config()`,
    `get_monte_carlo_defaults()` -> `default_monte_carlo_config()`,
    `get_predefined_properties()` -> `predefined_properties()`, and siblings.
  * `calculate_*()` -> `compute_*()` (`calculate_mode()`, `calculate_confidence_intervals()`, ...).
  * Bayesian verbs unified under `fuse_`/`update_`: `bayes_fuse()` -> `fuse_distribution()`,
    `bayes_update_normal_normal()` -> `fuse_normal_normal()`, `bayesian_update()` -> `update_prior()`.
  * `_multi` -> `_multiproperty`; noun-first names flipped verb-first
    (`mukey_draws_lookup()` -> `lookup_mukey_draws()`).
  * QA reporters moved from `validate_*`/`assess_*` to `diagnose_*`
    (`validate_monte_carlo_quality()` -> `diagnose_simulation()`, etc.); `validate_*` now means an
    assertion that returns invisibly or throws.

## Data infilling consolidation

* **One property-value cleaner.** `clean_property_data()` gains an `outlier_policy` argument
  (`"soil_aware"` default, `"aggressive_iqr"`, `"none"`) and is now called by both
  `process_ssurgo_data()` and `infill_soil_property()`.
  `clean_property_data_ssurgo_compatible()` is **deprecated** - a shim forwarding with
  `outlier_policy = "aggressive_iqr"`.
* **`process_ssurgo_data()` now cleans with the soil-aware outlier policy** (physically-impossible
  values only, conservative `IQR x 5` otherwise) instead of a generic `IQR x 3` on every property.
  The old policy nulled genuine distribution-tail values that the infilling then re-estimated,
  narrowing the very distributions the Monte Carlo is meant to characterise. This changes
  `process_ssurgo_data()` output for AOIs with legitimately wide-spread properties.
* **Correctness fix**: `simulate_cokey_generalized()` no longer silently disables texture
  simulation for horizons whose `hzname` doesn't map to a master horizon (e.g. the generic
  `H1`/`H2`/`H3` designations some SSURGO components use). An `NA` generalized-horizon value used
  as a list key is write-only, which fed a `NULL` correlation matrix into the texture draw and
  produced whole map units of `NA` in the raster-fusion output. Such horizons now fall back to the
  pooled genhz-agnostic correlation matrices.
* `infill_ssurgo_data()` (the raster-fusion infiller) now mirrors
  `process_soil_properties_comprehensive()`'s phase structure and derives `wthirdbar`/`wfifteenbar`
  via the Saxton-Rawls pedotransfer function (`water_retention_method = "saxton_rawls"`, default)
  rather than the generic recovery hierarchy.
* `learn_property_ranges()` no longer errors on a tibble whose property has an `_r` column but no
  `_l`/`_h` columns.
* New internal constants: `DEFAULT_MAX_DEPTH_CM` (250), `RANGE_FALLBACK_HALFWIDTH` (2),
  `SAXTON_RAWLS_RANGE_FRAC` (0.15), replacing repeated bare literals.

## Other changes

* Multi-property raster fusion. `simulate_cokey_generalized()`, `simulate_ssurgo_mapunit_draws()`,
  `fetch_ssurgo_percentiles()`, and `extract_mukey_joint_ensemble()` gain a `requested_properties`
  argument that restricts the SSURGO Monte Carlo simulation to the properties actually needed - the
  per-cokey depth-trend GP fit is the pipeline's dominant cost and scales with the property count.
* `run_fusion()` now simulates **only** its own property (plus the full sand/silt/clay draw
  for texture members) by default. The prior is statistically equivalent to - not bit-identical to -
  the previous all-property simulation (the pipeline was already unseeded), and a single-property
  map can show marginally more `NA` where a map unit's components carry no data for that one
  property. Restriction is honoured only under the default `vertical_correlation_method =
  "joint_copula"`; `"gp_quantile_retrofit"` falls back to a full simulation with a warning.
* New `run_fusion_multiproperty()`: fuses many properties over an AOI (and multiple depth windows)
  from a **single** SSURGO simulation, seeding the same per-property disk caches
  `run_fusion()` reads. Replaces an N-properties x M-windows loop of `run_fusion()`
  calls (each of which re-ran the full simulation); every leaf now shares one draw set, so their
  `NA` masks are mutually consistent.
* New optional `seed` argument on `simulate_ssurgo_mapunit_draws()`, `fetch_ssurgo_percentiles()`,
  `extract_mukey_joint_ensemble()`, `run_fusion()`, `run_fusion_group()`, and
  `run_fusion_multiproperty()` for opt-in reproducibility (default `NULL` = unchanged stochastic
  behavior). Seeds the sequential RNG stream and the parallel depth-trend `future.seed`;
  determinism is conditional on the live SSURGO/SOLUS data being unchanged.
* `simulate_ssurgo_mapunit_draws()`'s raw SSURGO tabular-download cache is now keyed by AOI alone
  instead of AOI + depth window (`download_ssurgo_tabular()` never depth-filters, so the cached
  content was already depth-window-independent) - one AOI's tabular download is now shared across
  every depth window/call requested for it, including between `run_fusion_multiproperty()`'s
  wide-span simulation and a later single-window `run_fusion()` call.
* **Correctness fix**: `simulate_cokey_generalized()`'s simulated `soc` property is now a genuine
  soil-organic-carbon estimate (`om * OM_TO_SOC_FACTOR`, the Van Bemmelen factor) instead of raw
  SSURGO organic matter passed through unconverted. Previously every `"soc"` fusion result (via
  `run_fusion()`/`run_fusion_multiproperty()` with `solus_variable = "soc"`) fused SSURGO
  organic matter % (prior) against SOLUS100's real organic carbon % (likelihood) - two different
  quantities differing by roughly a factor of 1.724. `"soc"` percentile/posterior values from this
  release will differ from earlier releases; the change should be read as a bug fix, not a
  regression. No change was needed to the KSSL reference correlation matrix (correlation is
  invariant under a linear rescale of one variable).
* New `fetch_solus_low_pred_high_multi()` / `fetch_solus_percentiles_multi()`: fetch every needed
  SOLUS100 variable's low/prediction/high rasters for one depth window in a **single**
  `soilDB::fetchSOLUS()` call (`variables`, `depth_slices`, and `output_type` are all requested as
  vectors at once), instead of one call per `(variable, output_type)`. `run_fusion_multiproperty()`
  now fetches SOLUS via this batched path - one call per window covering every standalone
  property and texture-group member, down from up to `3 * n_variables` calls per window - falling
  back to the existing per-variable scalar fetch only if the batched request itself errors. The
  scalar `fetch_solus_low_pred_high()` / `fetch_solus_percentiles()` (used by single-property
  `run_fusion()`) are unchanged.
* SSURGO chemistry-property coverage extended to 5 new properties: `caco3`, `ec`, `ecec`,
  `gypsum`, `sar` (SSURGO chorizon column names live-confirmed against a real gSSURGO query).
  `download_ssurgo_tabular()`'s default `properties`, `create_ssurgo_property_lookup_working()`,
  `simulate_cokey_generalized()`'s `param_order`, `property_to_sim_column()`,
  `normalize_requested_properties()`, `SSURGO_SIM_PROPERTY_COLUMNS`, and `infill_ssurgo_data()` all
  now recognize them, so they can be fetched via `fetch_ssurgo_percentiles()` and fused via
  `run_fusion()`/`run_fusion_multiproperty()` like any of the original 9 properties. They
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
  (`simulate_monte_carlo()`), Gaussian-process depth-trend modeling
  (`fit_depth_gp_models()`, `apply_depth_gp_to_simulation()`), and Bayesian
  updating/fusion in both scalar (`bayes_update_normal_normal()` and family,
  `bayesian_update()`) and raster-native (`fuse_property_adaptive()`, `fuse_texture_group()`)
  forms.
* SSURGO data-source adapter: tabular acquisition (`fetch_ssurgo_data()`), cleaning
  (`process_ssurgo_data()`), infilling (`process_soil_properties_comprehensive()`), and
  simulation-ready raster/percentile extraction (`simulate_ssurgo_mapunit_draws()`).
* SOLUS100 data-source adapter: raster percentile fetch (`fetch_solus_percentiles()`) supplying
  the likelihood side of the raster fusion pipeline.
* Profile, component, and depth simulation: component composition
  (`simulate_component_composition()`, `simulate_cokey_generalized()`) and horizon depth/thickness
  simulation (`simulate_and_perturb_soil_profiles()`, `simulate_profile_depths_by_mukey()`).
* Available water storage modeling via Van Genuchten / ROSETTA
  (`van_genuchten()`, `simulate_vg_aws()`, `compute_aws()`).
* Statistics and workflow validation/diagnostics (`analyze_soil_statistics()`,
  `diagnose_workflow()`, `generate_validation_report()`).
* Vignettes covering the tabular Monte Carlo pipeline, profile/depth simulation, available water
  storage, multi-source raster fusion (SSURGO x SOLUS100), per-pixel fused ensembles, and
  function-by-function tours of each subsystem, all built against real SSURGO data for a Sierra
  Nevada foothills and Salinas Valley area of interest.
