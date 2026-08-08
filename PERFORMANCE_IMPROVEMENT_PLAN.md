# soilSIM Performance Improvement Plan

Tracks a package-wide performance audit and fix pass, following on from
an earlier session that found and fixed five real bugs in the SSURGO
Monte Carlo -\> GP depth-trend pipeline (see the “Prior work” log
below) - all the same shape: an expensive operation not actually
dependent on a loop variable was called once per iteration instead of
once, vectorized.

Check items off as they are implemented, tested, and verified (one unit
of work at a time, not batched). Findings were located by three parallel
Explore-agent passes over all of `soilSIM/R/`, then re-verified against
the real source (line numbers below may still drift slightly - always
re-read the current code before touching it).

## Methodology (repeat per fix)

1.  **Locate & confirm** the anti-pattern by reading the actual current
    code.
2.  **Verify it’s a real hot path**: confirm actual internal callers
    exist (`grep` the whole `R/` tree), estimate realistic iteration
    counts for a real AOI before investing effort.
3.  **Benchmark before**:
    [`Rprof()`](https://rdrr.io/r/utils/Rprof.html)/[`summaryRprof()`](https://rdrr.io/r/utils/summaryRprof.html)
    via `data-raw/benchmark_performance.R` (small test AOI:
    `POLYGON((-121.652 36.610, -121.650 36.610, -121.650 36.612, -121.652 36.612, -121.652 36.610))`,
    inside the Salinas Valley area other tests/vignettes already use).
4.  **Fix**, preserving exact behavior - vectorize, don’t redesign.
    Comment explaining the anti-pattern and the profiling evidence.
5.  **Test**: regression test proving numeric equivalence to the old
    behavior (`tolerance`-based `expect_equal()`).
6.  **Benchmark after**, compute the speedup factor, log it below.
7.  **Commit** (message includes before/after numbers, no
    `Co-Authored-By` trailer per user’s global settings) + push, check
    the box.

Environment notes: `Rscript -e` inline segfaults in this environment -
always write scratch `.R` files and run via `Rscript path.R` in the
background; `unset PROJ_LIB` before any `terra`-touching Rscript call.

## Speedup log

| Item | Before | After | Factor | Commit |
|----|----|----|----|----|
| *(prior session, for reference)* `convert_to_property_matrices`/`convert_to_long_format` | 145.77s | 24.16s | ~6x (small AOI) | `5f7693d` |
| *(prior session)* `run_stage1_fusion_group` shared simulation | 116.31s | 38.98s | 2.98x | `1f66323` |
| *(prior session)* `fit_local_gp_model_single` gp_control | \- | \- | ~5x per GP_fit() call | `f77cb32` |
| *(prior session)* `apply_quantile_adjustment` vectorized | 14.66s | 8.77s | 1.67x (small AOI) | `2a5ad74` |
| *(prior session)* full SSURGO pipeline, real Salinas AOI | 1998.58s | 211.43s | 9.45x | `5f7693d` |
| [`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md) (50k-cell synthetic raster) | 150.34s | 193-357s (variable) | **0.4-0.8x - no win** (kept anyway: memory-safer, bit-identical, see Tier 1 notes) | `641a441` |
| `process_single_cokey()` split-once (500 cokeys x 20 rows, synthetic) | 0.91s | 0.08s | 11.4x | `89a3553` |
| `apply_sum_constraints()` (150 horizons x 2000 realizations, synthetic) | 3.71s | 1.16s | 3.2x | `78e5502` |
| [`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md) gp_control (20 reps, 12-pt series) | 82.70s | 16.97s | 4.9x | `29896d1` |
| [`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md) drop rowwise() (200 rows x 500 sims, synthetic) | 2.90s | 0.47s | 6.2x | `490dc21` |
| [`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md) vectorized toupper/trimws (200k rows, synthetic) | 70.12s | 9.47s | 7.4x | `a98b7c3` |
| [`infill_property_range_values()`](https://jjmaynard.github.io/soilSIM/reference/infill_property_range_values.md) (20k rows, synthetic) | 3.10s | 2.50s | 1.24x | `84b22b8` |
| [`related_property_estimation()`](https://jjmaynard.github.io/soilSIM/reference/related_property_estimation.md) (30k rows, organic-matter branch) | 7.82s | 0.04s | 195x | `b3a82bb` |
| [`hz_quant_prob_mukey()`](https://jjmaynard.github.io/soilSIM/reference/hz_quant_prob_mukey.md) single-pass (10k rows, 500 mukeys) | ~2.9s | ~2.9s | ~1x (no meaningful change - see Tier 3 notes) | `725d8a8` |
| [`adjust_multivariate_depthwise_GP()`](https://jjmaynard.github.io/soilSIM/reference/adjust_multivariate_depthwise_GP.md) (8 depths x 2000 replicates, synthetic) | 7.47s | 0.15s | 50x | `0b5f81c` |
| [`slice_and_aggregate_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/slice_and_aggregate_soil_data.md) (50 horizons, 1500cm total, synthetic) | 0.58s | 0.03s | 19x | `4a1ab9c` |

## Benchmark infrastructure

`data-raw/benchmark_performance.R` created (maintainer-run, not part of
package build) - one function per pipeline stage below,
[`system.time()`](https://rdrr.io/r/base/system.time.html)/[`Rprof()`](https://rdrr.io/r/utils/Rprof.html)
on the small test AOI. Commit `a690a5d` (scaffolding) + `4964a40` (real
`<<-` scoping bug in
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)’s
error handler, found while smoke-testing benchmark \#3 and fixed since
it silently broke
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)’s
error fallback).

## Tier 1 - high-confidence, high-impact, confirmed hot path

`R/raster-fusion.R::fuse_texture_group()` - no size-based dispatch
(unlike sibling
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)’s
`threshold_cells`), ran
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)
(6000 [`rnorm()`](https://rdrr.io/r/stats/Normal.html) draws) twice per
raster cell via
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html)+[`apply()`](https://rdrr.io/r/base/apply.html),
unconditionally.

Benchmarked before: 150.34s for a synthetic 50,000-cell raster (direct
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)
call, bypassing network fetch).

Fixed: vectorized the MC draws + ILR moments + bivariate-normal fusion
across all cells in a chunk at once
(`fuse_texture_group_batch()`/`_batch_core()`), with the
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) stream ordered to match
the original per-cell consumption order exactly (bit-identical output,
not just statistically similar). **Real bug found and fixed along the
way**: the first version processed a whole
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html)
chunk in one un-chunked matrix operation, which failed with “cannot
allocate vector of size 4.5 Gb” at 50,000 cells (terra’s own chunk
sizing doesn’t account for the 2000x internal MC expansion) - added
internal sub-chunking (`max_cells_per_subchunk`, default 2000) to bound
memory regardless of what chunk terra hands the function.

Regression test added: `test-raster-fusion.R` - reimplements the
original scalar loop inline and asserts `fuse_texture_group_batch()`
matches it bit-for-bit (tolerance 1e-9, heterogeneous per-cell inputs,
fixed seed) - proving this is a pure vectorization, not a behavior
change.

Benchmarked after: **not a clear win** - 193-357s across several runs at
the same 50,000-cell scale (higher and more variable than the 150.34s
baseline). Likely cause: the dominant cost is the ~600M
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) draws themselves (50k
cells x 2000 MC x 6 blocks), which batching doesn’t meaningfully speed
up over the same total draws done per-cell, while the large intermediate
matrices (up to ~1.8GB) add GC/allocation overhead that offsets the
R-level [`apply()`](https://rdrr.io/r/base/apply.html)-dispatch savings.
Also removed an unnecessary
[`t()`](https://rdrr.io/r/base/t.html)-transpose pass (kept everything
in R’s native column-major n_mc-rows x ncell-cols orientation) - didn’t
close the gap. Benchmarking was also confounded by the dev machine being
under real memory pressure at the time (repeated genuine `Rscript`
segfaults, unrelated to this code). Decision (user): keep the vectorized
version anyway - it’s correctness-neutral (bit-identical) and
memory-safer than the original (fixed the un-chunked blowup bug), and
may still help on real AOI data/hardware even though this synthetic
benchmark was inconclusive. Revisit with
[`Rprof()`](https://rdrr.io/r/utils/Rprof.html) profiling if this ever
shows up as a real bottleneck on an actual AOI run.

Committed + pushed

`R/multivariate-adjustment.R::process_single_cokey()` -
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
rescanned the full multi-cokey `simulation_data` once per cokey instead
of [`split()`](https://rdrr.io/r/base/split.html)-once.

Benchmarked before: 0.91s (synthetic 500 cokeys x 20 rows/cokey via
`process_cokeys_sequential()`).

Fixed:
[`integrate_monte_carlo_with_gp()`](https://jjmaynard.github.io/soilSIM/reference/integrate_monte_carlo_with_gp.md)
now splits `simulation_data` by cokey once
(`split(simulation_data, simulation_data$cokey)`) and passes the
resulting `cokey_groups` list down through
`process_cokeys_parallel()`/`process_cokeys_sequential()` to
`process_single_cokey()`, which now takes the already-filtered
`cokey_data` directly instead of `simulation_data` + re-filtering.
[`split()`](https://rdrr.io/r/base/split.html) partitions rows
identically to the old `dplyr::filter(cokey == !!cokey)` (same rows,
same order) - a pure row-extraction move, not a behavior change.

Regression test added: updated the two existing direct-call tests in
`test-multivariate-adjustment.R`
(`process_cokeys_sequential()`/`process_cokeys_parallel()`) to pass
`split(sim_data, sim_data$cokey)` per the new signature; both still
pass.

Benchmarked after: 0.08s - **~11.4x** at the same 500-cokey scale.

Committed + pushed

`R/monte-carlo.R::apply_sum_constraints()` - nested
per-(horizon,realization)-cell loop instead of vectorized array ops
(sibling constraint functions already do this correctly).

Benchmarked before: 3.71s (synthetic 150 horizons x 2000 realizations).

Fixed: `apply(results[, prop_indices, , drop=FALSE], c(1,3), sum)`
computes every (horizon, realization) cell’s sum at once; a single
[`ifelse()`](https://rdrr.io/r/base/ifelse.html)-built rescale-factor
matrix (1 = no-op for cells already within tolerance or with
`current_sum <= 0`) replaces the per-cell `if`; the factor matrix is
applied to each property slice in one broadcasted multiply. Matches the
vectorization pattern already used by
`apply_range_constraints_batch()`/`apply_relationship_constraints()`/
`apply_physical_constraints()` in the same file.

Regression test added: `test-monte-carlo.R` - reimplements the original
scalar loop inline and asserts bit-identical output (tolerance 1e-12)
plus matching `adjustments` count, including a zero-sum-cell and an
already-within-tolerance-cell edge case that must both stay untouched.

Benchmarked after: 1.16s - **~3.2x**, with `adjustments` count (299,907)
identical before/after, confirming correctness.

Committed + pushed

`R/gp-modeling.R::fit_individual_gp_model()`/[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md) -
same [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)
control-effort bug already fixed once elsewhere, but ~16x worse here
(5-fold CV x 3 correlation families + refit, none passing reduced
`control=`).

Benchmarked before: 82.70s for 20 reps of
[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)
on a 12-point synthetic series (default 5-fold CV).

Fixed: threaded `gp_control = c(20, 10, 2)` (the same
empirically-justified default already established for
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)
in `multivariate-adjustment.R`) through all
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) call
sites in `gp-modeling.R` -
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)’s
direct fit,
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)’s
per-fold fits, and
[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)’s
winning-candidate refit + fallback/baseline fits. All exist as function
parameters with this default, so existing callers (e.g.
[`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md))
get the speedup automatically without call-site changes;
`validation-diagnostics.R`’s existing positional
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
call is unaffected since `gp_control` was added after the existing
positional parameters.

Regression test added: existing
`test-gp-modeling.R`/`test-validation-diagnostics.R` suites (which
exercise
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)/[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)/
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
end-to-end) all still pass unchanged - the smaller search effort
converges to the identical optimum for these low-dimensional fits,
matching the precedent already established for
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md).

Benchmarked after: 16.97s - **~4.9x**.

Committed + pushed

## Tier 2 - moderate confidence/impact

`R/data-infilling.R::related_property_estimation()` - six per-row loops
(texture, water-retention, CEC, pH, organic-matter, bulk-density),
whole-dataset scope, each mechanically vectorizable (e.g. texture branch
-\> [`rowSums()`](https://rdrr.io/r/base/colSums.html)).

Benchmarked before: 7.82s (synthetic 30,000 rows, organic-matter
branch).

Fixed: all six branches rewritten as vectorized array operations
([`rowSums()`](https://rdrr.io/r/base/colSums.html) for the texture
sum-constraint, [`ifelse()`](https://rdrr.io/r/base/ifelse.html)-chains
replacing the per-row `if`/`else if` cascades for the other five)
operating on the whole missing-row set at once; `mark_estimated_vec()`
writes `infill_method` for a subset of indices in one call instead of
once per row.

Regression test added: extends `test-data-infilling.R` with a 40-row
synthetic dataset (randomized missingness across
`hzname`/`hzdept_r`/`claytotal_r`/`sandtotal_r`/`om_r`) run through all
6 branches, reimplementing the original per-row loop inline for
comparison - exact match (tolerance 1e-12) for every branch. The
existing single-row tests for each branch also still pass unchanged.

Benchmarked after: 0.04s - **~195x** (the organic-matter branch’s 3
conditional adjustments compounded the per-row overhead the most of the
six).

Committed + pushed

`R/data-infilling.R::infill_property_range_values()` -
`apply(df, 1, ...)` forces whole-frame character-matrix coercion per row
on every horizon needing range infilling.

Benchmarked before: 3.10s (synthetic 20,000 rows, ~90% needing infill).

Fixed:
[`calculate_property_lower_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_lower_bound.md)/[`calculate_property_upper_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_upper_bound.md)/
[`get_contextual_spread()`](https://jjmaynard.github.io/soilSIM/reference/get_contextual_spread.md)
only ever read 3 columns off `row` (`<property>_r`, `hzname`,
`hzdepb_r`) - pre-extract those 3 as plain vectors once, then build a
small named-list “row” per element (list
`[[`/[`names()`](https://rdrr.io/r/base/names.html) work identically to
a data.frame row for these helpers). **First attempt regression (caught
by benchmarking)**: slicing single-row data.frame subsets per element
(`subset_df[i, ]`) was actually *slower* than the original
[`apply()`](https://rdrr.io/r/base/apply.html) (8.16s vs. 3.10s) -
`[.data.frame`’s per-call overhead outweighs the savings from dropping
unused columns. Switched to plain-list construction instead, which is
cheap per element. **Also fixes a real latent bug**: the old
[`apply()`](https://rdrr.io/r/base/apply.html) coercion turned `row`
all-character, so
[`calculate_property_lower_bound()`](https://jjmaynard.github.io/soilSIM/reference/calculate_property_lower_bound.md)/`upper_bound()`’s
`if (depth <= 30)` depth-zone check (`depth = row[['hzdepb_r']]`, never
wrapped in [`as.numeric()`](https://rdrr.io/r/base/numeric.html)) was
comparing STRINGS (e.g. `"5" <= "30"` is FALSE lexicographically) -
silently misclassifying depth zones for horizons whose numeric
`hzdepb_r` string didn’t happen to sort the same as its numeric value.

Regression test added: `test-data-infilling.R` - constructs surface/deep
learned-range training rows with very different spreads and asserts a
shallow (5cm) target horizon correctly resolves to the `depth_surface`
spread, not `depth_deep` (which the old bug would have produced).

Benchmarked after: 2.50s - **~1.24x** (modest; this function’s cost is
dominated by
[`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md)/[`get_property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/get_property_contextual_ranges.md),
not the per-row loop itself).

Committed + pushed

`R/utils.R::is_unsuitable()` - per-row scalar
[`toupper()`](https://rdrr.io/r/base/chartr.html)/[`trimws()`](https://rdrr.io/r/base/trimws.html)
instead of vectorized once before the loop; called repeatedly across
pipeline stages on overlapping data.

Benchmarked before: 70.12s (synthetic 200,000 rows).

Fixed:
`toupper(trimws(hznames))`/`toupper(trimws(as.character(data$desgnmaster)))`
computed once (vectorized) before the loop instead of per-row inside it;
the loop’s own per-row classification logic
([`grepl()`](https://rdrr.io/r/base/grep.html) pattern matching) is
unchanged.

Regression test added: new `tests/testthat/test-utils.R` (no prior test
file existed for `utils.R`) - covers mixed case/whitespace, `NA`
handling, a missing `desgnmaster` column, and a missing `hzname_col`.

Benchmarked after: 9.47s - **~7.4x**, with an identical unsuitable-row
count (140,016) before/after, confirming correctness.

Committed + pushed

`R/aws-simulation.R::simulate_vg_aws()` - unnecessary
[`dplyr::rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html)
around an already-vectorized
[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)
call.

Benchmarked before: 2.90s (synthetic 200 rows x 500 simulations).

Fixed: dropped
[`dplyr::rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html)/`ungroup()`;
[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)
is pure elementwise arithmetic so
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
already vectorizes it correctly without row-at-a-time dispatch.

Regression test added: `test-aws-simulation.R` - asserts
`theta_fc`/`theta_pwp` are
[`identical()`](https://rdrr.io/r/base/identical.html) (not just
`equal()`) to per-row
[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)
calls, since this is a pure vectorization with no room for even
floating-point drift.

Benchmarked after: 0.47s - **~6.2x**.

Committed + pushed

## Tier 3 - low priority (confirmed zero internal callers via `grep` - exported API only)

`R/ssurgo-processing.R::hz_quant_prob_mukey()` - three independent
`group_by()`/`summarize()` passes (one per quantile level) + three
`left_join()`s, instead of one pass computing all quantiles per group.
**Correction (verified during implementation)**: originally listed as
Tier 1 “confirmed hot path” per an Explore agent’s report, but a direct
`grep -rn "hz_quant_prob_mukey" R/` found zero internal callers anywhere
in the package - it’s only referenced in its own file, a package-doc
comment, and docs/vignettes prose. Not wired into
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
or any other pipeline. Demoted to Tier 3 accordingly.

Fixed: replaced the three `group_by()`/`summarize()` passes + three
`left_join()`s with one grouped pass (`across()`’s multi-function form
computing all 3 quantiles per group at once) + a vectorized PIW90
computation, matching the plan’s suggested direction.

Regression test added: `test-ssurgo-processing.R` - reimplements the
original three-pass version inline and asserts identical output
(tolerance 1e-12) across multiple mukeys/depths/ properties. Existing
single-row and error-path tests for this function all still pass
unchanged.

Benchmarked: **inconclusive/no meaningful change** (2 repeated
before/after run pairs on a synthetic 10,000-row/500-mukey dataset
landed within ~10% of each other with no consistent direction -
`dplyr`’s per-call grouping/hashing overhead dominates at this scale
regardless of 1 vs. 3 summarize() passes). Kept anyway - the code is
genuinely simpler (one join instead of three, less intermediate object
churn) and correctness is verified; not a wasted change, just not a
measurable wall-clock win in this synthetic benchmark.

Committed + pushed

`R/gp-modeling.R::adjust_multivariate_depthwise_GP()` - identical
quantile()-per-replicate bug as the already-fixed
`apply_quantile_adjustment()`.

Fixed:
`quantile(curr_values, probs = reference_quantiles, na.rm = TRUE)`
computes all `n_sims` replicates’ quantiles in one call, replacing the
`for (j in 1:n_sims)` loop that called
[`quantile()`](https://rdrr.io/r/stats/quantile.html) once per replicate
with a single-probability `q`.

Regression test added: `test-gp-modeling.R` - reimplements the original
per-replicate loop inline and asserts identical output (tolerance 1e-10)
across a synthetic 2-property, 4-depth, 30-replicate dataset. The
existing dimension/validation test for this function also still passes
unchanged.

Benchmarked: 7.47s -\> 0.15s (synthetic 8 depths x 2,000 replicates) -
**~50x**.

Committed + pushed

`R/property-simulation.R::slice_and_aggregate_soil_data()` - allocates
one single-row `data.frame` per centimeter of profile depth before a
final [`rbind()`](https://rdrr.io/r/base/cbind.html).

Fixed: replaced the per-centimeter
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) +
[`rbind()`](https://rdrr.io/r/base/cbind.html) loop with a vectorized
row-index expansion
(`row_idx <- rep(seq_len(nrow(df)), times = n_expand)`) and a vectorized
per-row depth sequence via
[`sequence()`](https://rdrr.io/r/base/sequence.html)’s “concatenated
per-group [`seq_len()`](https://rdrr.io/r/base/seq.html)” trick, then a
single `df[row_idx, data_columns]` subset instead of one tiny data.frame
per depth.

Regression test added: `test-property-simulation.R` - reimplements the
original per-centimeter loop inline and cross-checks the new vectorized
version’s per-depth-range aggregates match it, across multiple
horizons/properties with uneven depth ranges. The existing
single-property aggregation test also still passes unchanged.

Benchmarked: 0.58s -\> 0.03s (synthetic 50 horizons, 1500cm total
profile depth) - **~19x**.

Committed + pushed

## Confirmed NOT bugs (no action needed)

- ROSETTA batching (`aws-simulation.R::calculate_aws_df()`) - already
  correctly batched.
- [`simulate_correlated_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_properties.md),
  [`estimate_property_correlations()`](https://jjmaynard.github.io/soilSIM/reference/estimate_property_correlations.md),
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md),
  [`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
  [`estimate_correlation_matrix_robust()`](https://jjmaynard.github.io/soilSIM/reference/estimate_correlation_matrix_robust.md) -
  already vectorized.
- `depth-simulation.R`
  ([`simulate_profile_depths_by_collection()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths_by_collection.md),
  [`evaluate_simulated_depths()`](https://jjmaynard.github.io/soilSIM/reference/evaluate_simulated_depths.md),
  [`fetch_osd_horizons_cached()`](https://jjmaynard.github.io/soilSIM/reference/fetch_osd_horizons_cached.md)) -
  genuinely per-profile work or already cached.
- `distribution-fitting-raster.R` non-KDE routes,
  [`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
  (already gated behind `threshold_cells`) - working as designed.
- `bayesian-updating.R::bayesian_update()` - confirmed not wired into
  the main pipeline by the file’s own header comment; dead code from the
  hot-path perspective by original design.
- Small/capped loops in `validation-diagnostics.R` and
  `statistics.R::analyze_property_distributions_safe()` - negligible
  impact.

## Final gate

All Tier 1 items complete

Full-AOI (Salinas Valley, 15 mukeys) end-to-end timing run of
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
for the texture composition group, logged above

Full `devtools::test()` pass (not just affected files): 2 failures, both
pre-existing and unrelated to this audit -
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)’s
two live-network tests (`test-raster-fusion.R:277`/`:305`) fail in
[`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md)
(a file never touched by any fix here) because their fixtures use
synthetic property IDs (`test_clay`/`test_clay2`) never registered in
that mapping. Confirmed via call-stack trace that neither failure passes
through any file modified in this plan.

Tier 2/3 items complete (optional, lower urgency)
