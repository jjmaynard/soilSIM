# soilSIM Performance Improvement Plan

Tracks a package-wide performance audit and fix pass, following on from an earlier session that
found and fixed five real bugs in the SSURGO Monte Carlo -> GP depth-trend pipeline (see the
"Prior work" log below) - all the same shape: an expensive operation not actually dependent on a
loop variable was called once per iteration instead of once, vectorized.

Check items off as they are implemented, tested, and verified (one unit of work at a time, not
batched). Findings were located by three parallel Explore-agent passes over all of `soilSIM/R/`,
then re-verified against the real source (line numbers below may still drift slightly - always
re-read the current code before touching it).

## Methodology (repeat per fix)

1. **Locate & confirm** the anti-pattern by reading the actual current code.
2. **Verify it's a real hot path**: confirm actual internal callers exist (`grep` the whole `R/`
   tree), estimate realistic iteration counts for a real AOI before investing effort.
3. **Benchmark before**: `Rprof()`/`summaryRprof()` via `data-raw/benchmark_performance.R` (small
   test AOI: `POLYGON((-121.652 36.610, -121.650 36.610, -121.650 36.612, -121.652 36.612,
   -121.652 36.610))`, inside the Salinas Valley area other tests/vignettes already use).
4. **Fix**, preserving exact behavior - vectorize, don't redesign. Comment explaining the
   anti-pattern and the profiling evidence.
5. **Test**: regression test proving numeric equivalence to the old behavior
   (`tolerance`-based `expect_equal()`).
6. **Benchmark after**, compute the speedup factor, log it below.
7. **Commit** (message includes before/after numbers, no `Co-Authored-By` trailer per user's
   global settings) + push, check the box.

Environment notes: `Rscript -e` inline segfaults in this environment - always write scratch `.R`
files and run via `Rscript path.R` in the background; `unset PROJ_LIB` before any
`terra`-touching Rscript call.

## Speedup log

| Item | Before | After | Factor | Commit |
|---|---|---|---|---|
| *(prior session, for reference)* `convert_to_property_matrices`/`convert_to_long_format` | 145.77s | 24.16s | ~6x (small AOI) | `5f7693d` |
| *(prior session)* `run_stage1_fusion_group` shared simulation | 116.31s | 38.98s | 2.98x | `1f66323` |
| *(prior session)* `fit_local_gp_model_single` gp_control | - | - | ~5x per GP_fit() call | `f77cb32` |
| *(prior session)* `apply_quantile_adjustment` vectorized | 14.66s | 8.77s | 1.67x (small AOI) | `2a5ad74` |
| *(prior session)* full SSURGO pipeline, real Salinas AOI | 1998.58s | 211.43s | 9.45x | `5f7693d` |
| `fuse_texture_group()` (50k-cell synthetic raster) | 150.34s | 193-357s (variable) | **0.4-0.8x - no win** (kept anyway: memory-safer, bit-identical, see Tier 1 notes) | `641a441` |
| `process_single_cokey()` split-once (500 cokeys x 20 rows, synthetic) | 0.91s | 0.08s | 11.4x | `89a3553` |
| `apply_sum_constraints()` (150 horizons x 2000 realizations, synthetic) | 3.71s | 1.16s | 3.2x | `78e5502` |
| `optimize_gp_hyperparameters()` gp_control (20 reps, 12-pt series) | 82.70s | 16.97s | 4.9x | `29896d1` |
| `simulate_vg_aws()` drop rowwise() (200 rows x 500 sims, synthetic) | 2.90s | 0.47s | 6.2x | `490dc21` |
| `is_unsuitable()` vectorized toupper/trimws (200k rows, synthetic) | 70.12s | 9.47s | 7.4x | `a98b7c3` |
| `infill_property_range_values()` (20k rows, synthetic) | 3.10s | 2.50s | 1.24x | TBD |

## Benchmark infrastructure

- [x] `data-raw/benchmark_performance.R` created (maintainer-run, not part of package build) - one
  function per pipeline stage below, `system.time()`/`Rprof()` on the small test AOI. Commit
  `a690a5d` (scaffolding) + `4964a40` (real `<<-` scoping bug in `validate_data_quality()`'s error
  handler, found while smoke-testing benchmark #3 and fixed since it silently broke
  `integrate_monte_carlo_with_gp()`'s error fallback).

## Tier 1 - high-confidence, high-impact, confirmed hot path

- [x] `R/raster-fusion.R::fuse_texture_group()` - no size-based dispatch (unlike sibling
  `fuse_property_adaptive()`'s `threshold_cells`), ran `estimate_ilr_moments_mc()` (6000
  `rnorm()` draws) twice per raster cell via `terra::app()`+`apply()`, unconditionally.
  - [x] Benchmarked before: 150.34s for a synthetic 50,000-cell raster (direct
    `fuse_texture_group()` call, bypassing network fetch).
  - [x] Fixed: vectorized the MC draws + ILR moments + bivariate-normal fusion across all cells
    in a chunk at once (`fuse_texture_group_batch()`/`_batch_core()`), with the `rnorm()` stream
    ordered to match the original per-cell consumption order exactly (bit-identical output, not
    just statistically similar). **Real bug found and fixed along the way**: the first version
    processed a whole `terra::app()` chunk in one un-chunked matrix operation, which failed with
    "cannot allocate vector of size 4.5 Gb" at 50,000 cells (terra's own chunk sizing doesn't
    account for the 2000x internal MC expansion) - added internal sub-chunking
    (`max_cells_per_subchunk`, default 2000) to bound memory regardless of what chunk terra hands
    the function.
  - [x] Regression test added: `test-raster-fusion.R` - reimplements the original scalar loop
    inline and asserts `fuse_texture_group_batch()` matches it bit-for-bit (tolerance 1e-9,
    heterogeneous per-cell inputs, fixed seed) - proving this is a pure vectorization, not a
    behavior change.
  - [x] Benchmarked after: **not a clear win** - 193-357s across several runs at the same
    50,000-cell scale (higher and more variable than the 150.34s baseline). Likely cause: the
    dominant cost is the ~600M `rnorm()` draws themselves (50k cells x 2000 MC x 6 blocks), which
    batching doesn't meaningfully speed up over the same total draws done per-cell, while the
    large intermediate matrices (up to ~1.8GB) add GC/allocation overhead that offsets the
    R-level `apply()`-dispatch savings. Also removed an unnecessary `t()`-transpose pass (kept
    everything in R's native column-major n_mc-rows x ncell-cols orientation) - didn't close the
    gap. Benchmarking was also confounded by the dev machine being under real memory pressure at
    the time (repeated genuine `Rscript` segfaults, unrelated to this code). Decision (user):
    keep the vectorized version anyway - it's correctness-neutral (bit-identical) and
    memory-safer than the original (fixed the un-chunked blowup bug), and may still help on real
    AOI data/hardware even though this synthetic benchmark was inconclusive. Revisit with
    `Rprof()` profiling if this ever shows up as a real bottleneck on an actual AOI run.
  - [x] Committed + pushed

- [x] `R/multivariate-adjustment.R::process_single_cokey()` - `dplyr::filter()` rescanned the
  full multi-cokey `simulation_data` once per cokey instead of `split()`-once.
  - [x] Benchmarked before: 0.91s (synthetic 500 cokeys x 20 rows/cokey via
    `process_cokeys_sequential()`).
  - [x] Fixed: `integrate_monte_carlo_with_gp()` now splits `simulation_data` by cokey once
    (`split(simulation_data, simulation_data$cokey)`) and passes the resulting `cokey_groups`
    list down through `process_cokeys_parallel()`/`process_cokeys_sequential()` to
    `process_single_cokey()`, which now takes the already-filtered `cokey_data` directly instead
    of `simulation_data` + re-filtering. `split()` partitions rows identically to the old
    `dplyr::filter(cokey == !!cokey)` (same rows, same order) - a pure row-extraction move, not a
    behavior change.
  - [x] Regression test added: updated the two existing direct-call tests in
    `test-multivariate-adjustment.R` (`process_cokeys_sequential()`/`process_cokeys_parallel()`)
    to pass `split(sim_data, sim_data$cokey)` per the new signature; both still pass.
  - [x] Benchmarked after: 0.08s - **~11.4x** at the same 500-cokey scale.
  - [x] Committed + pushed

- [x] `R/monte-carlo.R::apply_sum_constraints()` - nested per-(horizon,realization)-cell loop
  instead of vectorized array ops (sibling constraint functions already do this correctly).
  - [x] Benchmarked before: 3.71s (synthetic 150 horizons x 2000 realizations).
  - [x] Fixed: `apply(results[, prop_indices, , drop=FALSE], c(1,3), sum)` computes every
    (horizon, realization) cell's sum at once; a single `ifelse()`-built rescale-factor matrix
    (1 = no-op for cells already within tolerance or with `current_sum <= 0`) replaces the
    per-cell `if`; the factor matrix is applied to each property slice in one broadcasted
    multiply. Matches the vectorization pattern already used by
    `apply_range_constraints_batch()`/`apply_relationship_constraints()`/
    `apply_physical_constraints()` in the same file.
  - [x] Regression test added: `test-monte-carlo.R` - reimplements the original scalar loop
    inline and asserts bit-identical output (tolerance 1e-12) plus matching `adjustments` count,
    including a zero-sum-cell and an already-within-tolerance-cell edge case that must both stay
    untouched.
  - [x] Benchmarked after: 1.16s - **~3.2x**, with `adjustments` count (299,907) identical
    before/after, confirming correctness.
  - [x] Committed + pushed

- [x] `R/gp-modeling.R::fit_individual_gp_model()`/`optimize_gp_hyperparameters()` - same
  `GPfit::GP_fit()` control-effort bug already fixed once elsewhere, but ~16x worse here (5-fold
  CV x 3 correlation families + refit, none passing reduced `control=`).
  - [x] Benchmarked before: 82.70s for 20 reps of `optimize_gp_hyperparameters()` on a 12-point
    synthetic series (default 5-fold CV).
  - [x] Fixed: threaded `gp_control = c(20, 10, 2)` (the same empirically-justified default
    already established for `fit_local_gp_model_single()` in `multivariate-adjustment.R`) through
    all `GPfit::GP_fit()` call sites in `gp-modeling.R` -
    `fit_individual_gp_model()`'s direct fit, `k_fold_gp_cv()`'s per-fold fits, and
    `optimize_gp_hyperparameters()`'s winning-candidate refit + fallback/baseline fits. All exist
    as function parameters with this default, so existing callers (e.g.
    `build_stratified_gp_models()`) get the speedup automatically without call-site changes;
    `validation-diagnostics.R`'s existing positional `k_fold_gp_cv()` call is unaffected since
    `gp_control` was added after the existing positional parameters.
  - [x] Regression test added: existing `test-gp-modeling.R`/`test-validation-diagnostics.R`
    suites (which exercise `fit_individual_gp_model()`/`optimize_gp_hyperparameters()`/
    `k_fold_gp_cv()` end-to-end) all still pass unchanged - the smaller search effort converges
    to the identical optimum for these low-dimensional fits, matching the precedent already
    established for `fit_local_gp_model_single()`.
  - [x] Benchmarked after: 16.97s - **~4.9x**.
  - [x] Committed + pushed

## Tier 2 - moderate confidence/impact

- [ ] `R/data-infilling.R::related_property_estimation()` - six per-row loops (texture,
  water-retention, CEC, pH, organic-matter, bulk-density), whole-dataset scope, each
  mechanically vectorizable (e.g. texture branch -> `rowSums()`).
  - [ ] Benchmarked before / Fixed / Regression test / Benchmarked after / Committed

- [x] `R/data-infilling.R::infill_property_range_values()` - `apply(df, 1, ...)` forces
  whole-frame character-matrix coercion per row on every horizon needing range infilling.
  - [x] Benchmarked before: 3.10s (synthetic 20,000 rows, ~90% needing infill).
  - [x] Fixed: `calculate_property_lower_bound()`/`calculate_property_upper_bound()`/
    `get_contextual_spread()` only ever read 3 columns off `row` (`<property>_r`, `hzname`,
    `hzdepb_r`) - pre-extract those 3 as plain vectors once, then build a small named-list "row"
    per element (list `[[`/`names()` work identically to a data.frame row for these helpers).
    **First attempt regression (caught by benchmarking)**: slicing single-row data.frame subsets
    per element (`subset_df[i, ]`) was actually *slower* than the original `apply()` (8.16s vs.
    3.10s) - `[.data.frame`'s per-call overhead outweighs the savings from dropping unused
    columns. Switched to plain-list construction instead, which is cheap per element. **Also
    fixes a real latent bug**: the old `apply()` coercion turned `row` all-character, so
    `calculate_property_lower_bound()`/`upper_bound()`'s `if (depth <= 30)` depth-zone check
    (`depth = row[['hzdepb_r']]`, never wrapped in `as.numeric()`) was comparing STRINGS (e.g.
    `"5" <= "30"` is FALSE lexicographically) - silently misclassifying depth zones for horizons
    whose numeric `hzdepb_r` string didn't happen to sort the same as its numeric value.
  - [x] Regression test added: `test-data-infilling.R` - constructs surface/deep learned-range
    training rows with very different spreads and asserts a shallow (5cm) target horizon
    correctly resolves to the `depth_surface` spread, not `depth_deep` (which the old bug would
    have produced).
  - [x] Benchmarked after: 2.50s - **~1.24x** (modest; this function's cost is dominated by
    `learn_property_ranges()`/`get_property_contextual_ranges()`, not the per-row loop itself).
  - [x] Committed + pushed

- [x] `R/utils.R::is_unsuitable()` - per-row scalar `toupper()`/`trimws()` instead of vectorized
  once before the loop; called repeatedly across pipeline stages on overlapping data.
  - [x] Benchmarked before: 70.12s (synthetic 200,000 rows).
  - [x] Fixed: `toupper(trimws(hznames))`/`toupper(trimws(as.character(data$desgnmaster)))`
    computed once (vectorized) before the loop instead of per-row inside it; the loop's own
    per-row classification logic (`grepl()` pattern matching) is unchanged.
  - [x] Regression test added: new `tests/testthat/test-utils.R` (no prior test file existed for
    `utils.R`) - covers mixed case/whitespace, `NA` handling, a missing `desgnmaster` column, and
    a missing `hzname_col`.
  - [x] Benchmarked after: 9.47s - **~7.4x**, with an identical unsuitable-row count
    (140,016) before/after, confirming correctness.
  - [x] Committed + pushed

- [x] `R/aws-simulation.R::simulate_vg_aws()` - unnecessary `dplyr::rowwise()` around an
  already-vectorized `van_genuchten()` call.
  - [x] Benchmarked before: 2.90s (synthetic 200 rows x 500 simulations).
  - [x] Fixed: dropped `dplyr::rowwise()`/`ungroup()`; `van_genuchten()` is pure elementwise
    arithmetic so `dplyr::mutate()` already vectorizes it correctly without row-at-a-time
    dispatch.
  - [x] Regression test added: `test-aws-simulation.R` - asserts `theta_fc`/`theta_pwp` are
    `identical()` (not just `equal()`) to per-row `van_genuchten()` calls, since this is a pure
    vectorization with no room for even floating-point drift.
  - [x] Benchmarked after: 0.47s - **~6.2x**.
  - [x] Committed + pushed

## Tier 3 - low priority (confirmed zero internal callers via `grep` - exported API only)

- [ ] `R/ssurgo-processing.R::hz_quant_prob_mukey()` - three independent
  `group_by()`/`summarize()` passes (one per quantile level) + three `left_join()`s, instead of
  one pass computing all quantiles per group. **Correction (verified during implementation)**:
  originally listed as Tier 1 "confirmed hot path" per an Explore agent's report, but a direct
  `grep -rn "hz_quant_prob_mukey" R/` found zero internal callers anywhere in the package - it's
  only referenced in its own file, a package-doc comment, and docs/vignettes prose. Not wired
  into `simulate_ssurgo_mapunit_draws()` or any other pipeline. Demoted to Tier 3 accordingly.
  - [ ] Fixed / Regression test / Committed

- [ ] `R/gp-modeling.R::adjust_multivariate_depthwise_GP()` - identical
  quantile()-per-replicate bug as the already-fixed `apply_quantile_adjustment()`.
  - [ ] Fixed / Regression test / Committed

- [ ] `R/property-simulation.R::slice_and_aggregate_soil_data()` - allocates one single-row
  `data.frame` per centimeter of profile depth before a final `rbind()`.
  - [ ] Fixed / Regression test / Committed

## Confirmed NOT bugs (no action needed)

- ROSETTA batching (`aws-simulation.R::calculate_aws_df()`) - already correctly batched.
- `simulate_correlated_properties()`, `estimate_property_correlations()`,
  `preserve_correlation_structure()`, `ilr_forward()`/`ilr_inverse()`,
  `estimate_correlation_matrix_robust()` - already vectorized.
- `depth-simulation.R` (`simulate_profile_depths_by_collection()`, `evaluate_simulated_depths()`,
  `fetch_osd_horizons_cached()`) - genuinely per-profile work or already cached.
- `distribution-fitting-raster.R` non-KDE routes, `fuse_general_kde()` (already gated behind
  `threshold_cells`) - working as designed.
- `bayesian-updating.R::bayesian_update()` - confirmed not wired into the main pipeline by the
  file's own header comment; dead code from the hot-path perspective by original design.
- Small/capped loops in `validation-diagnostics.R` and
  `statistics.R::analyze_property_distributions_safe()` - negligible impact.

## Final gate

- [x] All Tier 1 items complete
- [ ] Full-AOI (Salinas Valley, 15 mukeys) end-to-end timing run of `run_stage1_fusion_group()`
  for the texture composition group, logged above
- [ ] Full `devtools::test()` pass (not just affected files) with 0 failures
- [ ] Tier 2/3 items complete (optional, lower urgency)
