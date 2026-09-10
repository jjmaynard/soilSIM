# Changelog

## soilSIM 0.1.0.9000 (development version)

### Data infilling consolidation

- **One property-value cleaner.**
  [`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md)
  gains an `outlier_policy` argument (`"soil_aware"` default,
  `"aggressive_iqr"`, `"none"`) and is now called by both
  [`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)
  and
  [`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md).
  [`clean_property_data_ssurgo_compatible()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data_ssurgo_compatible.md)
  is **deprecated** - a shim forwarding with
  `outlier_policy = "aggressive_iqr"`.
- **[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)
  now cleans with the soil-aware outlier policy** (physically-impossible
  values only, conservative `IQR x 5` otherwise) instead of a generic
  `IQR x 3` on every property. The old policy nulled genuine
  distribution-tail values that the infilling then re-estimated,
  narrowing the very distributions the Monte Carlo is meant to
  characterise. This changes
  [`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md)
  output for AOIs with legitimately wide-spread properties.
- **Correctness fix**:
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
  no longer silently disables texture simulation for horizons whose
  `hzname` doesn’t map to a master horizon (e.g. the generic
  `H1`/`H2`/`H3` designations some SSURGO components use). An `NA`
  generalized-horizon value used as a list key is write-only, which fed
  a `NULL` correlation matrix into the texture draw and produced whole
  map units of `NA` in the raster-fusion output. Such horizons now fall
  back to the pooled genhz-agnostic correlation matrices.
- [`infill_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_data.md)
  (the raster-fusion infiller) now mirrors
  [`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md)’s
  phase structure and derives `wthirdbar`/`wfifteenbar` via the
  Saxton-Rawls pedotransfer function
  (`water_retention_method = "saxton_rawls"`, default) rather than the
  generic recovery hierarchy.
- [`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md)
  no longer errors on a tibble whose property has an `_r` column but no
  `_l`/`_h` columns.
- New internal constants: `DEFAULT_MAX_DEPTH_CM` (250),
  `RANGE_FALLBACK_HALFWIDTH` (2), `SAXTON_RAWLS_RANGE_FRAC` (0.15),
  replacing repeated bare literals.

### Other changes

- Multi-property raster fusion.
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md),
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
  and
  [`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
  gain a `requested_properties` argument that restricts the SSURGO Monte
  Carlo simulation to the properties actually needed - the per-cokey
  depth-trend GP fit is the pipeline’s dominant cost and scales with the
  property count.
- [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
  now simulates **only** its own property (plus the full sand/silt/clay
  draw for texture members) by default. The prior is statistically
  equivalent to - not bit-identical to - the previous all-property
  simulation (the pipeline was already unseeded), and a single-property
  map can show marginally more `NA` where a map unit’s components carry
  no data for that one property. Restriction is honoured only under the
  default `vertical_correlation_method = "joint_copula"`;
  `"gp_quantile_retrofit"` falls back to a full simulation with a
  warning.
- New
  [`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md):
  fuses many properties over an AOI (and multiple depth windows) from a
  **single** SSURGO simulation, seeding the same per-property disk
  caches
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
  reads. Replaces an N-properties x M-windows loop of
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
  calls (each of which re-ran the full simulation); every leaf now
  shares one draw set, so their `NA` masks are mutually consistent.
- New optional `seed` argument on
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
  [`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md),
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md),
  [`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md),
  and
  [`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)
  for opt-in reproducibility (default `NULL` = unchanged stochastic
  behavior). Seeds the sequential RNG stream and the parallel
  depth-trend `future.seed`; determinism is conditional on the live
  SSURGO/SOLUS data being unchanged.
- [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)’s
  raw SSURGO tabular-download cache is now keyed by AOI alone instead of
  AOI + depth window
  ([`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
  never depth-filters, so the cached content was already
  depth-window-independent) - one AOI’s tabular download is now shared
  across every depth window/call requested for it, including between
  [`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)’s
  wide-span simulation and a later single-window
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
  call.
- **Correctness fix**:
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)’s
  simulated `soc` property is now a genuine soil-organic-carbon estimate
  (`om * OM_TO_SOC_FACTOR`, the Van Bemmelen factor) instead of raw
  SSURGO organic matter passed through unconverted. Previously every
  `"soc"` fusion result (via
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)
  with `solus_variable = "soc"`) fused SSURGO organic matter % (prior)
  against SOLUS100’s real organic carbon % (likelihood) - two different
  quantities differing by roughly a factor of 1.724. `"soc"`
  percentile/posterior values from this release will differ from earlier
  releases; the change should be read as a bug fix, not a regression. No
  change was needed to the KSSL reference correlation matrix
  (correlation is invariant under a linear rescale of one variable).
- New
  [`fetch_solus_low_pred_high_multi()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high_multi.md)
  /
  [`fetch_solus_percentiles_multi()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles_multi.md):
  fetch every needed SOLUS100 variable’s low/prediction/high rasters for
  one depth window in a **single**
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
  call (`variables`, `depth_slices`, and `output_type` are all requested
  as vectors at once), instead of one call per
  `(variable, output_type)`.
  [`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)
  now fetches SOLUS via this batched path - one call per window covering
  every standalone property and texture-group member, down from up to
  `3 * n_variables` calls per window - falling back to the existing
  per-variable scalar fetch only if the batched request itself errors.
  The scalar
  [`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
  /
  [`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)
  (used by single-property
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md))
  are unchanged.
- SSURGO chemistry-property coverage extended to 5 new properties:
  `caco3`, `ec`, `ecec`, `gypsum`, `sar` (SSURGO chorizon column names
  live-confirmed against a real gSSURGO query).
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)’s
  default `properties`,
  [`create_ssurgo_property_lookup_working()`](https://jjmaynard.github.io/soilSIM/reference/create_ssurgo_property_lookup_working.md),
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)’s
  `param_order`,
  [`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md),
  [`normalize_requested_properties()`](https://jjmaynard.github.io/soilSIM/reference/normalize_requested_properties.md),
  `SSURGO_SIM_PROPERTY_COLUMNS`, and
  [`infill_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_data.md)
  all now recognize them, so they can be fetched via
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)
  and fused via
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_multi()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_multi.md)
  like any of the original 9 properties. They have no entry in the
  static KSSL reference correlation matrices (fit years before these
  properties existed in this pipeline) -
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  now builds its per-genhz correlation matrices via
  [`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)
  (sized to the full property vocabulary), so these 5 simulate
  **uncorrelated** with everything else (and each other) rather than
  erroring; real correlations can be added later via a KSSL data re-fit
  with no code changes.
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
  also gained a defensive fallback
  ([`extend_corr_matrix_with_identity()`](https://jjmaynard.github.io/soilSIM/reference/extend_corr_matrix_with_identity.md))
  so a direct caller supplying an older correlation-matrix list that
  predates this change degrades to the same independence default instead
  of erroring. SSURGO coverage for these 5 properties is real but
  sparser than the original 9 in spot-checks (e.g. `ecec` was populated
  in only ~25% of sampled horizons for one AOI) - expect more `NA` in
  their fused output.
- New optional `restriction_depth` argument on
  [`remarginalized_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalized_awc.md):
  bedrock/restriction-depth- aware AWC truncation. SOLUS100 predicts
  every property at every depth even for shallow soils, so without this,
  a pixel with bedrock at 40 cm was credited with water storage in depth
  windows entirely or partially below bedrock - a real, physically-wrong
  AWC inflation for shallow soils. `NULL` (default) preserves prior
  output exactly; when supplied (typically
  `fetch_solus_restriction_depth(aoi_vect)`, new below), each window’s
  thickness contribution becomes the per-pixel effective thickness
  `pmax(0, pmin(bottom, restriction_depth) - top)` - a window straddling
  the restriction gets partial credit, one entirely below it contributes
  exactly 0 (not `NA`, which would incorrectly blank the whole pixel’s
  AWC rather than just that window).
- New
  [`fetch_solus_site_level()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_site_level.md)
  /
  [`fetch_solus_restriction_depth()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_restriction_depth.md):
  fetch a SOLUS100 site-level (depth-independent) variable -
  `anylithicdpt` (depth to bedrock) / `resdept` (depth to restriction) -
  in one `fetchSOLUS(depth_slices = "all", ...)` call.
  [`fetch_solus_restriction_depth()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_restriction_depth.md)
  defaults to `anylithicdpt`, confirmed (per SOLUS100’s own published
  training methodology) to be trained only on lithic/paralithic bedrock
  restriction records - hard-bedrock-specific by construction, needing
  no separate SSURGO-side classifier; `resdept` (trained on every
  restriction record, any kind) is offered as a broader, more
  conservative non-default alternative. Both are right-censored at 201
  cm in SOLUS100’s training data (new `SOLUS_RESTRICTION_CENSOR_CM`
  constant) - values at/above that are recoded to `Inf` (no restriction
  detected), not treated as a literal depth.

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
- Vignettes covering the tabular Monte Carlo pipeline, profile/depth
  simulation, available water storage, multi-source raster fusion
  (SSURGO x SOLUS100), per-pixel fused ensembles, and
  function-by-function tours of each subsystem, all built against real
  SSURGO data for a Sierra Nevada foothills and Salinas Valley area of
  interest.
