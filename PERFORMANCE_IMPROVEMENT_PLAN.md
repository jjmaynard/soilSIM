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
| `infill_property_range_values()` (20k rows, synthetic) | 3.10s | 2.50s | 1.24x | `84b22b8` |
| `related_property_estimation()` (30k rows, organic-matter branch) | 7.82s | 0.04s | 195x | `b3a82bb` |
| `hz_quant_prob_mukey()` single-pass (10k rows, 500 mukeys) | ~2.9s | ~2.9s | ~1x (no meaningful change - see Tier 3 notes) | `725d8a8` |
| `adjust_multivariate_depthwise_GP()` (8 depths x 2000 replicates, synthetic) | 7.47s | 0.15s | 50x | `0b5f81c` |
| `slice_and_aggregate_soil_data()` (50 horizons, 1500cm total, synthetic) | 0.58s | 0.03s | 19x | `4a1ab9c` |

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

- [x] `R/data-infilling.R::related_property_estimation()` - six per-row loops (texture,
  water-retention, CEC, pH, organic-matter, bulk-density), whole-dataset scope, each
  mechanically vectorizable (e.g. texture branch -> `rowSums()`).
  - [x] Benchmarked before: 7.82s (synthetic 30,000 rows, organic-matter branch).
  - [x] Fixed: all six branches rewritten as vectorized array operations (`rowSums()` for the
    texture sum-constraint, `ifelse()`-chains replacing the per-row `if`/`else if` cascades for
    the other five) operating on the whole missing-row set at once; `mark_estimated_vec()`
    writes `infill_method` for a subset of indices in one call instead of once per row.
  - [x] Regression test added: extends `test-data-infilling.R` with a 40-row synthetic dataset
    (randomized missingness across `hzname`/`hzdept_r`/`claytotal_r`/`sandtotal_r`/`om_r`) run
    through all 6 branches, reimplementing the original per-row loop inline for comparison -
    exact match (tolerance 1e-12) for every branch. The existing single-row tests for each
    branch also still pass unchanged.
  - [x] Benchmarked after: 0.04s - **~195x** (the organic-matter branch's 3 conditional
    adjustments compounded the per-row overhead the most of the six).
  - [x] Committed + pushed

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

- [x] `R/ssurgo-processing.R::hz_quant_prob_mukey()` - three independent
  `group_by()`/`summarize()` passes (one per quantile level) + three `left_join()`s, instead of
  one pass computing all quantiles per group. **Correction (verified during implementation)**:
  originally listed as Tier 1 "confirmed hot path" per an Explore agent's report, but a direct
  `grep -rn "hz_quant_prob_mukey" R/` found zero internal callers anywhere in the package - it's
  only referenced in its own file, a package-doc comment, and docs/vignettes prose. Not wired
  into `simulate_ssurgo_mapunit_draws()` or any other pipeline. Demoted to Tier 3 accordingly.
  - [x] Fixed: replaced the three `group_by()`/`summarize()` passes + three `left_join()`s with
    one grouped pass (`across()`'s multi-function form computing all 3 quantiles per group at
    once) + a vectorized PIW90 computation, matching the plan's suggested direction.
  - [x] Regression test added: `test-ssurgo-processing.R` - reimplements the original three-pass
    version inline and asserts identical output (tolerance 1e-12) across multiple mukeys/depths/
    properties. Existing single-row and error-path tests for this function all still pass
    unchanged.
  - [x] Benchmarked: **inconclusive/no meaningful change** (2 repeated before/after run pairs on
    a synthetic 10,000-row/500-mukey dataset landed within ~10% of each other with no consistent
    direction - `dplyr`'s per-call grouping/hashing overhead dominates at this scale regardless
    of 1 vs. 3 summarize() passes). Kept anyway - the code is genuinely simpler (one join instead
    of three, less intermediate object churn) and correctness is verified; not a wasted change,
    just not a measurable wall-clock win in this synthetic benchmark.
  - [x] Committed + pushed

- [x] `R/gp-modeling.R::adjust_multivariate_depthwise_GP()` - identical
  quantile()-per-replicate bug as the already-fixed `apply_quantile_adjustment()`.
  - [x] Fixed: `quantile(curr_values, probs = reference_quantiles, na.rm = TRUE)` computes all
    `n_sims` replicates' quantiles in one call, replacing the `for (j in 1:n_sims)` loop that
    called `quantile()` once per replicate with a single-probability `q`.
  - [x] Regression test added: `test-gp-modeling.R` - reimplements the original per-replicate
    loop inline and asserts identical output (tolerance 1e-10) across a synthetic 2-property,
    4-depth, 30-replicate dataset. The existing dimension/validation test for this function also
    still passes unchanged.
  - [x] Benchmarked: 7.47s -> 0.15s (synthetic 8 depths x 2,000 replicates) - **~50x**.
  - [x] Committed + pushed

- [x] `R/property-simulation.R::slice_and_aggregate_soil_data()` - allocates one single-row
  `data.frame` per centimeter of profile depth before a final `rbind()`.
  - [x] Fixed: replaced the per-centimeter `as.data.frame()` + `rbind()` loop with a vectorized
    row-index expansion (`row_idx <- rep(seq_len(nrow(df)), times = n_expand)`) and a vectorized
    per-row depth sequence via `sequence()`'s "concatenated per-group `seq_len()`" trick, then a
    single `df[row_idx, data_columns]` subset instead of one tiny data.frame per depth.
  - [x] Regression test added: `test-property-simulation.R` - reimplements the original
    per-centimeter loop inline and cross-checks the new vectorized version's per-depth-range
    aggregates match it, across multiple horizons/properties with uneven depth ranges. The
    existing single-property aggregation test also still passes unchanged.
  - [x] Benchmarked: 0.58s -> 0.03s (synthetic 50 horizons, 1500cm total profile depth) -
    **~19x**.
  - [x] Committed + pushed

## Confirmed NOT bugs (no action needed)

- ROSETTA batching (`aws-simulation.R::calculate_aws_df()`) - already correctly batched.
- `simulate_correlated_properties()`, `estimate_property_correlations()`,
  `preserve_correlation_structure()`, `ilr_forward()`/`ilr_inverse()`,
  `estimate_correlation_matrix_robust()` - already vectorized.
- `depth-simulation.R` (`simulate_profile_depths_by_collection()`, `evaluate_simulated_depths()`,
  `fetch_osd_horizons_cached()`) - genuinely per-profile work or already cached.
- `distribution-fitting-raster.R` non-KDE routes - working as designed.
- `bayesian-updating.R::bayesian_update()` - confirmed not wired into the main pipeline by the
  file's own header comment; dead code from the hot-path perspective by original design.
- Small/capped loops in `validation-diagnostics.R` and
  `statistics.R::analyze_property_distributions_safe()` - negligible impact.

## Tier 4 - KDE fusion + new row-loop candidates (found by comparing against a sibling Python
project's own performance report, `soil-id-algorithm-api/docs/PERFORMANCE_PROFILING_REPORT.md`)

**Correction to the "Confirmed NOT bugs" entry above**: `fuse_general_kde()` was previously waved
off as "already gated behind `threshold_cells` - working as designed" with no benchmark behind
that call. It was not fine - see below.

- [x] `R/raster-fusion.R::fuse_general_kde()` - per-cell `bayesian_update()` call (itself calling
  `stats::density()` twice) inside `terra::app()` + `apply(row_mat, 1, ...)`, for every cell at or
  below `threshold_cells` (default 80,000).
  - [x] Benchmarked before (synthetic percentile rasters, bypassing network fetch): 5.88s / 29.58s
    / 109.94s / 209.13s at 506 / 2,024 / 10,000 / 20,022 actual cells respectively (~linear,
    ~10.4ms/cell) - extrapolating to the default `threshold_cells = 80,000` implies ~14 minutes in
    the worst case.
  - [x] `Rprof()` profiling (10,000-cell case) found `bayesian_update()` = 77% of total time, of
    which `stats::density()` alone = 65% (`dnorm` 26%, `fft` 16%) - `simulate_from_percentiles()`
    (the piece originally assumed to be the main target) was only 12%, of which
    `extract_percentile_pairs()` dispatch overhead alone was 9% (more than the actual
    `approxfun()`-based sampling itself). This reprioritized the fix: `density()`'s evaluation
    grid length (driven by `bayesian_update()`'s `grid_resolution` parameter, default `0.01`) was
    the real lever, not per-cell R dispatch.
  - [x] Accuracy-vs-speed sweep across 4 representative percentile scenarios (narrow/wide/skewed
    distributions) comparing `grid_resolution` against a `0.01` reference: at `0.1`, max mean error
    0.018%, max variance error 0.12% (analytical, sampling-noise-free comparison) - both far inside
    the ~5-10% variability typical of field-measured soil properties (the same bar the sibling
    Python report's cubic-spline-to-linear-interpolation swap was validated against). Coarser than
    `0.1` degrades fast (variance error reaches double digits by `0.25`, 470% by `2.0`).
  - [x] Fixed: added `FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION <- 0.1` (`R/raster-fusion.R`) and
    changed `fuse_general_kde()`'s `NULL`-`grid_resolution` fallback to use it instead of silently
    deferring to `bayesian_update()`'s own standalone default of `0.01`. `bayesian_update()`'s own
    default is deliberately left unchanged - only `fuse_general_kde()`/`fuse_adaptive()`'s
    caller-facing default shifted, so any other direct caller of `bayesian_update()` is unaffected.
  - [x] Regression test added: `test-raster-fusion.R` - compares mean posterior mu/sigma (averaged
    over 8 seeds, to isolate the systematic `grid_resolution` effect from `bayesian_update()`'s own
    two-layer Monte Carlo sampling noise) between the new default and the `0.01` reference, across
    the same 3 representative scenarios, within a tolerance looser than the analytical bound
    (2%/10%) to stay robust while still catching a real regression.
  - [x] Benchmarked after: 3.65s / 10.35s / 49.17s / 84.92s at the same 4 scales - **1.6x-2.9x**
    (largest at the 2,024-cell scale).
  - [x] Full `devtools::test()`: 0 failures (only pre-existing live-network skips/unrelated
    warnings). Full `devtools::check()`: 0 errors, 0 warnings, 1 pre-existing unrelated NOTE
    (untracked `pkgdown/` directory, predates this work).
  - [x] Committed (`c1b64fb`); not yet pushed.

- [x] `simulate_from_percentiles()`/`extract_percentile_pairs()` batching across cells - originally
  planned as the primary fix before profiling reprioritized it as secondary; turned out bigger
  than the ~10-12% estimate.
  - [x] Added `sim_linear_cdf_batch()` (`R/percentile-sampling.R`) - vectorized equivalent of
    calling `sim_linear_cdf()` once per row via `stats::approxfun()`, using `findInterval()` +
    matrix arithmetic across all rows sharing the same `probs` knots at once. Validated against
    `stats::approxfun()` at fixed evaluation points (bit-identical, max abs diff = 0) and
    distributionally equivalent to `sim_linear_cdf()` via KS test (both checks now in
    `test-percentile-sampling.R`).
  - [x] Wired into `fuse_general_kde()`: cells with no NA/`-1`-sentinel percentile values (the
    common case) now go through a batched fast path (`sim_linear_cdf_batch()` for both prior and
    likelihood, once per chunk); any cell with a missing value falls back to the original
    per-cell `simulate_from_percentiles()` path unchanged, preserving the existing
    NA-degrades-to-NA-output contract exactly (confirmed by the existing dedicated test still
    passing). `bayesian_update()`'s `density()` call still runs per cell - no vectorized form
    exists in base R.
  - [x] Benchmarked after (on top of the `grid_resolution` fix): 1.27s / 5.18s / 23.37s / 44.92s
    at the same 4 scales - a further **1.9x-2.9x**, for a **cumulative 4.6x-5.7x** from the
    original baseline (5.88s/29.58s/109.94s/209.13s). At the default `threshold_cells = 80,000`,
    this takes the original ~14-minute worst case down to roughly 3 minutes.
  - [x] Full `devtools::test()`: 0 failures. Full `devtools::check()`: 0 errors, 0 warnings, 1
    pre-existing unrelated NOTE.
  - [x] Committed (`a00a993`); not yet pushed.

- [x] `R/property-simulation.R::simulate_cokey_generalized()` - **assumption corrected by
  profiling before fixing**: the per-row Cholesky/multivariate draw itself
  (`simulate_correlated_triangular()`) turned out to only be 8.4% of total wall-clock. A full
  vectorization was not attempted, since `simulate_correlated_triangular()` takes one (a,b,c)
  triplet per property and draws `n` realizations from it, and both the triplets **and** `n`
  (`sim_comppct`) vary per row - true multi-row batching would need a substantial redesign of
  that function's contract. `Rprof()` profiling (2,000 synthetic rows) found the *real* cost
  elsewhere:
  - [x] `calculate_mode()` (called twice per row, in the texture step) = **22% of total
    wall-clock** via `table()`/`factor()`/`unique()` - a full factor/hash-table build just to
    count occurrences of up to `sim_comppct` (~15-60) continuous values.
  - [x] `ensure_positive_definite_matrix(txt_corr)` = **18% of total wall-clock** via
    `eigen()`/`isSymmetric.matrix()`, recomputed identically on every row despite being a
    deterministic function of `genhz_val` alone (typically a handful of distinct values per
    cokey, not one per row).
  - [x] Fixed `calculate_mode()`: replaced `table(x)`/`factor(x)` with
    `tabulate(match(x, sort(unique(x))))` - verified identical output (including `table()`'s
    implicit smallest-value tie-break) across ties, negatives, singletons, and random floats.
    Only caller is `simulate_cokey_generalized()`; has its own pre-existing unit test.
  - [x] Fixed the PD-matrix recomputation: added `pd_txt_corr_cache` (keyed by `genhz_val`)
    inside `simulate_cokey_generalized()`, computed once per distinct genhz value instead of
    once per row. Confirmed `list[[NA_character_]]` always returns `NULL` even after
    "storing" under that key (an R quirk) - the one row-level edge case (unparseable/NA genhz)
    simply never benefits from the cache and recomputes every time, same correct result, no bug.
  - [x] Regression test added: `test-property-simulation.R` - reimplements the ORIGINAL
    uncached-per-row version inline (as it existed before this fix) and asserts bit-identical
    output at the same seed, using a 2-row fixture that intentionally shares one `genhz` value
    (the case the caching actually exercises - the existing single-row fixtures never did).
  - [x] Benchmarked after: 4.74s -> 2.94s (2,000 synthetic rows) - **~1.6x**, from two low-risk
    fixes with zero change to the actual statistical/RNG behavior.
  - [x] Full `devtools::test()`: 0 failures. Full `devtools::check()`: 0 errors, 0 warnings, 1
    pre-existing unrelated NOTE.
  - [x] Committed (`d81248d`); not yet pushed.

- [x] `R/multivariate-adjustment.R::merge_adjusted_data()` - per-row `which()` full-table scan
  "join," reached per-cokey via `apply_gp_depth_trends()`.
  - [x] Benchmarked before: 3.22s (10,000 rows, 10 depths x 1,000 realizations - one cokey's
    scale).
  - [x] Fixed: replaced the per-row `which()` scan with a single vectorized key match
    (`paste(hzdept_r, simulation_number)` composite key, `match()` + a uniqueness check via
    `table()` to preserve the original's "exactly one match" contract - zero or duplicate
    matches in `result_data` are silently skipped, not an error). Verified empirically that R's
    vectorized `[<-` assignment applies duplicate indices in order with the last one winning
    (`x[c(2,2,3)] <- c(10,20,30)` -> `x[2] == 20`), matching the original sequential loop's
    last-write-wins behavior for `adjusted_data` rows that share a key.
  - [x] Regression test added: `test-multivariate-adjustment.R` - reimplements the original
    per-row loop inline and asserts bit-identical output, plus explicit spot-checks of every
    edge case (duplicate key within `result_data` -> never matched; duplicate key within
    `adjusted_data` -> last wins; `NA` value in `adjusted_data` -> skipped).
  - [x] Benchmarked after: 3.22s -> 0.07s - **~46x**.
  - [x] Full `devtools::test()`: 0 failures. Full `devtools::check()`: 0 errors, 0 warnings, 1
    pre-existing unrelated NOTE.
  - [x] Committed (`6af405f`); not yet pushed.

- [x] `R/multivariate-adjustment.R::apply_cross_property_constraints()` - per-row texture
  sum/rescale via `data[i, texture_props]` row-slicing. **Note**: `correct_distribution_shapes()`
  (this function's only caller) has zero internal call sites anywhere in `R/` (confirmed via
  `grep -rn "correct_distribution_shapes(" R/`) - exported-API-only, same status as the earlier
  Tier 3 `hz_quant_prob_mukey()` item. Fixed for correctness/consistency regardless.
  - [x] Benchmarked before: **31.61s** (50,000 synthetic rows) - much more expensive than its
    trivial per-row body suggested.
  - [x] Fixed: `rowSums()`-based vectorization over the whole texture-column matrix at once,
    mirroring the already-fixed `related_property_estimation()` texture branch exactly.
    `texture_sum > 0 & !is.na(texture_sum)` preserves the original's exact `&&`-based NA
    handling (`NA & FALSE` resolves to `FALSE` for both `&`/`&&`).
  - [x] Regression test added: `test-multivariate-adjustment.R` - reimplements the original
    per-row loop inline and asserts bit-identical output across NA-sum, zero-sum, and
    already-100 rows.
  - [x] Benchmarked after: 31.61s -> 0.02s - **~1,580x**.
  - [x] Full `devtools::test()`: 0 failures. Full `devtools::check()`: 0 errors, 0 warnings, 1
    pre-existing unrelated NOTE.
  - [ ] Committed + pushed

- [ ] `R/monte-carlo.R::check_property_data_availability()` - nested per-row/per-property NA check
  with repeated column re-indexing. Benchmarked before: 0.42s (20,000 rows x 5 properties) - real
  caller (`monte-carlo.R:1840`) but already cheap; lowest priority of the 5. Not yet fixed.

## Final gate

- [x] All Tier 1 items complete
- [ ] Full-AOI (Salinas Valley, 15 mukeys) end-to-end timing run of `run_stage1_fusion_group()`
  for the texture composition group, logged above
- [x] Full `devtools::test()` pass (not just affected files): 2 failures, both pre-existing and
  unrelated to this audit - `run_stage1_fusion()`'s two live-network tests
  (`test-raster-fusion.R:277`/`:305`) fail in `property_to_sim_column()` (a file never touched by
  any fix here) because their fixtures use synthetic property IDs (`test_clay`/`test_clay2`)
  never registered in that mapping. Confirmed via call-stack trace that neither failure passes
  through any file modified in this plan.
- [x] Tier 2/3 items complete (optional, lower urgency)
