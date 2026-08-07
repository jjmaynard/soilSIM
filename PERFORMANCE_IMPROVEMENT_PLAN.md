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
| `fuse_texture_group()` (50k-cell synthetic raster) | 150.34s | 193-357s (variable) | **0.4-0.8x - no win** (kept anyway: memory-safer, bit-identical, see Tier 1 notes) | TBD |

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

- [ ] `R/multivariate-adjustment.R::process_single_cokey()` - `dplyr::filter()` rescans full
  multi-cokey `simulation_data` once per cokey instead of `split()`-once.
  - [ ] Benchmarked before
  - [ ] Fixed
  - [ ] Regression test added
  - [ ] Benchmarked after, logged above
  - [ ] Committed + pushed

- [ ] `R/monte-carlo.R::apply_sum_constraints()` - nested per-(horizon,realization)-cell loop
  instead of vectorized array ops (sibling constraint functions already do this correctly).
  - [ ] Benchmarked before
  - [ ] Fixed
  - [ ] Regression test added
  - [ ] Benchmarked after, logged above
  - [ ] Committed + pushed

- [ ] `R/gp-modeling.R::fit_individual_gp_model()`/`optimize_gp_hyperparameters()` - same
  `GPfit::GP_fit()` control-effort bug already fixed once elsewhere, but ~16x worse here (5-fold
  CV x 3 correlation families + refit, none passing reduced `control=`).
  - [ ] Benchmarked before
  - [ ] Fixed (apply established `gp_control = c(20, 10, 2)` default to all call sites in file)
  - [ ] Regression test added
  - [ ] Benchmarked after, logged above
  - [ ] Committed + pushed

## Tier 2 - moderate confidence/impact

- [ ] `R/data-infilling.R::related_property_estimation()` - six per-row loops (texture,
  water-retention, CEC, pH, organic-matter, bulk-density), whole-dataset scope, each
  mechanically vectorizable (e.g. texture branch -> `rowSums()`).
  - [ ] Benchmarked before / Fixed / Regression test / Benchmarked after / Committed

- [ ] `R/data-infilling.R::infill_property_range_values()` - `apply(df, 1, ...)` forces
  whole-frame character-matrix coercion per row on every horizon needing range infilling.
  - [ ] Benchmarked before / Fixed / Regression test / Benchmarked after / Committed

- [ ] `R/utils.R::is_unsuitable()` - per-row scalar `toupper()`/`trimws()` instead of vectorized
  once before the loop; called repeatedly across pipeline stages on overlapping data.
  - [ ] Benchmarked before / Fixed / Regression test / Benchmarked after / Committed

- [ ] `R/aws-simulation.R::simulate_vg_aws()` - unnecessary `dplyr::rowwise()` around an
  already-vectorized `van_genuchten()` call.
  - [ ] Benchmarked before / Fixed / Regression test / Benchmarked after / Committed

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

- [ ] All Tier 1 items complete
- [ ] Full-AOI (Salinas Valley, 15 mukeys) end-to-end timing run of `run_stage1_fusion_group()`
  for the texture composition group, logged above
- [ ] Full `devtools::test()` pass (not just affected files) with 0 failures
- [ ] Tier 2/3 items complete (optional, lower urgency)
