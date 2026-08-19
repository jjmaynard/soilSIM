# Vertical Correlation Improvement Plan

Tracks the redesign of soilSIM's depth-to-depth (vertical) correlation modeling, replacing the
current GP-mean-trend + quantile-retrofit approach with a joint depth x property Gaussian-copula
sampler. Source plan (approved): `C:\Users\jmaynard\.claude\plans\evaluate-these-files-from-quirky-pretzel.md`.

Follows the same working conventions as `PERFORMANCE_IMPROVEMENT_PLAN.md`: one unit of work at a
time, test before moving to the next item, scratch `.R` files + background `Rscript` runs (inline
`Rscript -e` segfaults in this environment), re-read current code before touching it (line numbers
below may drift).

## Context (why this work exists)

A sibling Python project (`soil-id-algorithm-api`) documents a known flaw in its vertical
autocorrelation approach: independent per-property AR(1) applied *after* within-horizon Cholesky
correlation distorts cross-property correlation by up to ~50% when per-property rho values differ
a lot (`CORRELATION_PRESERVATION_ANALYSIS.md` in that repo). soilSIM already does something more
sophisticated - a per-property Gaussian Process (GPfit) depth-trend fit, cross-validated over
exponential/Matern kernels - but has the *same architectural flaw*: depth correlation and property
correlation are two separate sequential stages, reconciled only by a retrofit trick (one "primary
property"'s rank copied to every depth, nudged by the GP mean ratio, remapped to the local
marginal). The code itself flags this as approximate: it warns when post-hoc correlation drift
exceeds 0.1 (`gp-modeling.R:770`, pre-Phase-0) and only spot-checked the first 5 depths.

"Currently gaussian" (as asked) refers to the `GPfit` kernel family (exponential/Matern - a
Gaussian Process), not a separate legacy noise model - confirmed by grep, no AR(1)/gaussian-noise
alternative exists in `soilSIM/R`.

Goal: replace the quantile-retrofit adjustment with a **Gaussian-copula joint sampling step** that
draws depth correlation and property correlation simultaneously from a Kronecker-structured
covariance (depth kernel (x) property-correlation matrix), so both are satisfied by construction -
no retrofit, no drift warning - while reusing the valuable parts of the current design (GP-fitted
mean trends, existing property-correlation estimation, KSSL reference matrices).

## Status tracker

| Phase | Description | Status |
|---|---|---|
| 0 | Instrumentation - joint correlation validator | [x] Done |
| 1 | Depth kernel + length-scale extraction | [x] Done |
| 1b | Thread `bound_sd` (OSD boundary distinctness) through pipeline | [x] Done |
| 1c | Discontinuity-gated depth kernel | [x] Done |
| 1d | Gating regression/unit tests | [x] Done (combined with 1c) |
| 1e (optional) | Add `chorizon.hzdistinctness` to SDA fetch | [ ] Deferred |
| 2 | Joint copula sampler | [x] Done |
| 3 | `preserve_correlation_structure_joint()` | [x] Done |
| 4 | Config-flag wiring (`vertical_correlation_method`) | [x] Done |
| 5 | Validation on real data (Amador County AOI, via vignette) | [x] Done |
| 6 | Regression tests (old-path bit-identical + new-path fidelity) | [x] Done |
| 7 | Benchmark (cost/quality tradeoff) | [x] Done |
| 8 | Decouple core method from discontinuity gating (`vertical_correlation_gating` flag) | [x] Done |
| 9 | Multi-AOI validation (Salinas Valley + discontinuity profile) | [ ] Not started |
| 10 | Full-AOI performance benchmark | [ ] Not started |
| 11 | Thread `config`/`gp_models` through the NRCS GP path | [ ] Not started |
| 12 | Empirical check on GP-recentering for strictly-positive properties | [ ] Not started |
| 13 | Execute the default flip | [ ] Not started |

---

## Phase 0 - Instrumentation: joint correlation validator

**Goal**: a validation helper that reports achieved depth-lag correlation *and* achieved property
correlation against targets, at all depths (not just first 5) - the acceptance test every later
phase is judged against.

- [x] Added `validate_joint_correlation_structure()` to `R/gp-modeling.R` (after
  `validate_correlation_preservation()`, ~line 866). Takes `simulated_list` (named list of
  `[n_depths x n_sims]` matrices, same shape used throughout this pipeline), optional
  `target_property_corr` (k x k) and `target_depth_corr` (n_depths x n_depths). Returns:
  - `achieved_property_correlation` - list of one k x k matrix per depth (all depths, not 5)
  - `property_correlation_max_diff` - numeric vector, one per depth
  - `achieved_depth_correlation` - named list of one n_depths x n_depths matrix per property (new -
    `validate_correlation_preservation()` never checked this axis at all)
  - `depth_correlation_max_diff` - named numeric vector, one per property
  - `overall_property_max_diff`, `overall_depth_max_diff` - worst-case summaries
- [x] Exported in `NAMESPACE` (alphabetized, between `validate_integration_results` and
  `validate_monte_carlo_config`).
- [x] Regression tests added to `tests/testthat/test-gp-modeling.R`:
  - Constructs synthetic data via the exact separable/Kronecker-MVN sampling identity the Phase 2
    joint sampler will use (`Z = L_depth %*% Eps %*% t(L_prop)`, target correlation matrices known
    by construction) at `n_sims = 4000`; asserts achieved max diff from both targets is `< 0.1`
    (the same threshold `adjust_multivariate_depthwise_GP()` already uses for its own drift
    warning), and that a deliberately wrong property-correlation target is caught (`> 0.3` diff).
  - Confirms all-`NA` diffs (not an error) when no targets are supplied, and correct list
    shapes/names in both cases.
- [x] Full `tests/testthat/test-gp-modeling.R` run (`devtools::load_all()` +
  `testthat::test_file()`, via scratch `.R` file per this repo's `Rscript -e` segfault
  workaround): **26/26 test blocks pass, 0 failures, 0 errors** (24 pre-existing blocks + 2 new).

**Notes for future reference**: this function is deliberately *not* wired into any pipeline call
site yet - it is a standalone diagnostic, used starting in Phase 5 to A/B the old
(`gp_quantile_retrofit`) and new (`joint_copula`) methods on real Salinas Valley AOI data. Its
`0.1` pass/fail threshold is inherited from the existing (pre-Phase-0) drift-warning threshold in
`adjust_multivariate_depthwise_GP()` purely for continuity/comparability - not re-derived here, and
worth revisiting once Phase 5 has real achieved-distortion numbers to calibrate against.

---

## Phase 1 - Depth kernel + length-scale extraction

**Goal**: build `R_depth` (n_depths x n_depths correlation matrix) via an exponential/Matern kernel
over real depth distance, with its length-scale derived from the already-fitted `GPfit` models
(not a new estimation step).

- [x] `extract_depth_length_scale(gp_model_list, depth_scaling = NULL, target_corr = exp(-1))`
  added to `R/gp-modeling.R` (new "3b. VERTICAL-CORRELATION REDESIGN: DEPTH CORRELATION KERNEL"
  section, placed after the Phase 0 validator). Reuses `GPfit::corr_matrix()` (GPfit's *own*
  exported correlation-evaluation function - confirmed via `getNamespaceExports("GPfit")`) plus
  `stats::uniroot()` to solve for the scaled-`[0,1]`-units distance at which the already-fitted
  model's correlation function crosses `exp(-1)` (the standard "range" definition), then rescales
  by the model's own stored `depth_scaling$range` to real depth units (cm). Works generically
  across GPfit's correlation families (exponential, Matern) without hand-deriving a closed-form
  inverse per family, since it evaluates GPfit's actual fitted correlation function numerically
  rather than assuming a specific parametric form.
- [x] `build_depth_correlation_kernel(depths, length_scale, kernel = c("exponential","matern"), nu = 1.5)`
  added alongside it. Builds an `n_depths x n_depths` matrix from real Euclidean depth *distances*
  (handles uneven horizon spacing), via `exp(-dist/length_scale)` (exponential) or a standard
  Matern-nu correlation form (parameterized directly by `length_scale`, independent of GPfit's
  internal beta/power form). Falls back to the identity matrix (no depth correlation) on
  non-finite/non-positive `length_scale`, matching this package's established graceful-fallback
  pattern. Repairs floating-point asymmetry and confirms positive-definiteness via the existing
  `ensure_positive_definite_matrix()` (`R/distributions.R:667`) rather than a new PD-repair
  routine.
- [x] Both exported in `NAMESPACE` (alphabetized) and documented (`devtools::document()` generated
  `extract_depth_length_scale.Rd`/`build_depth_correlation_kernel.Rd` cleanly, no warnings).
- [x] Unit tests added to `tests/testthat/test-gp-modeling.R`:
  - `extract_depth_length_scale()`: fits a real `GPfit` model on 5 synthetic depths (0-150cm),
    confirms a positive, finite, real-units length-scale roughly commensurate with the depth
    range; confirms the scaled-`[0,1]`-units result (no `depth_scaling` supplied) rescales back to
    the same real-units value (`scaled * depth_scaling$range == real_units`, tolerance `1e-8`);
    confirms `NA_real_` (not an error) on `NULL`/malformed input.
  - `build_depth_correlation_kernel()`: for both `"exponential"` and `"matern"`, confirms
    unit-diagonal, symmetric, positive-semi-definite output at uneven depth spacing
    `c(0, 10, 30, 100)`, and confirms monotonic non-increasing correlation with distance from a
    common anchor depth; confirms identity-matrix fallback for `NA`/`0`/negative `length_scale`;
    confirms the single-depth degenerate case returns `matrix(1,1,1)`.
- [x] **Bug found and fixed during testing**: `stats::dist()` labels matrix dimnames `"1","2",...`
  by default (not `NULL`), which `ensure_positive_definite_matrix()` then preserves through its
  eigen-repair round-trip - `expect_equal(diag(R), rep(1,4))` failed on a names mismatch (values
  were correct; only `names(diag(R))` differed). Fixed by explicitly stripping dimnames off the
  distance matrix (`dimnames(dist_mat) <- NULL`) right after `stats::dist()`/`as.matrix()`, before
  any downstream use - a real latent gotcha worth remembering for any other code in this package
  that builds matrices via `stats::dist()`.
- [x] Full `tests/testthat/test-gp-modeling.R` run (scratch `.R` file, per this repo's
  `Rscript -e` segfault workaround): **102/102 assertions pass, 0 failures, 0 errors** (after the
  dimnames fix above; first run caught the bug with 3 failures, all in the new
  `build_depth_correlation_kernel()` test, confirming the test itself was doing real work).

**Notes for future reference**: `extract_depth_length_scale()` intentionally does NOT hand-roll
GPfit's correlation formula (which differs between the exponential/power family and Matern family,
and could silently drift out of sync with a future GPfit version) - it evaluates GPfit's own
exported `corr_matrix()` numerically and solves for the range via root-finding instead, so it stays
correct automatically as long as GPfit's `GP_fit()`/`corr_matrix()` contract doesn't change.
`build_depth_correlation_kernel()`'s Matern implementation is a *separate*, standard
geostatistical Matern form (parameterized directly by real-units `length_scale`) - deliberately not
reusing GPfit's internal Matern parameterization, since by the time this function is called the
length-scale is already in real depth-distance units, not GPfit's internal scaled/beta form.

## Phase 1b - Thread `bound_sd` through the property-simulation pipeline

**Goal**: expose the already-computed OSD-derived `bound_sd` boundary-distinctness column (today
computed in `depth-simulation.R::query_osd_distinctness()` and consumed only by that file's
`aqp::perturb()` call) on the `sim_cokey`/`cokey_data` frames that flow into
`simulate_cokey_generalized()` and `apply_gp_depth_trends()`.

- [x] Added `attach_osd_boundary_distinctness(hz_data)` to `R/depth-simulation.R` (right after
  `query_osd_distinctness()`). Requires `compname`/`hzname`/`genhz` columns; calls
  `query_osd_distinctness()` and reuses the EXACT SAME
  `dplyr::group_by(id, genhz) |> dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE))`
  aggregation `simulate_and_perturb_soil_profiles()` already used internally
  (`depth-simulation.R:722-725`), then left-joins the result onto `hz_data` by
  `(compname, genhz)`. On OSD lookup failure (no network, `soilDB::fetchOSD()` error, etc.),
  degrades to an all-`NA` `bound_sd` column via `tryCatch()` rather than erroring - callers
  (including the Phase 1c kernel gating, once built) already treat missing distinctness data as
  "no gating, fall back to the plain kernel."
- [x] Wired into the main production entry point,
  `R/ssurgo-simulation.R::simulate_ssurgo_mapunit_draws()`, immediately after
  `hz_data$genhz <- classify_genhz(hz_data$hzname)` and before `sim_component_comp()` - this is
  the function real AOI runs actually call (`simulate_cokey_generalized()`'s production caller;
  the `depth-simulation.R` `SoilProfileCollection`/`aqp::perturb()` path that originally computed
  `bound_sd` is a separate, parallel horizon-depth-perturbation workflow, not currently wired into
  this property-value simulation pipeline at all - confirmed via
  `grep -rn simulate_and_perturb_soil_profiles` showing zero callers inside
  `ssurgo-simulation.R`/`property-simulation.R`).
- [x] `R/property-simulation.R::simulate_cokey_generalized()` now copies `row$bound_sd` onto its
  output `sim_data$bound_sd` whenever the input `row` has a `bound_sd` column - guarded
  (`if ("bound_sd" %in% names(row))`), so every pre-existing caller whose `sim_cokey` has no
  `bound_sd` column is completely unaffected (confirmed by the full test suite - see below).
  `@return` roxygen doc updated to document the new optional column.
- [x] Unit tests added:
  - `tests/testthat/test-depth-simulation.R`: `attach_osd_boundary_distinctness()` tested fully
    offline via `testthat::local_mocked_bindings(query_osd_distinctness = ...)` (the same
    mocking pattern `test-gp-modeling.R` already uses for `predict_gp_depth_trends()`) - confirms
    correct per-`(compname, genhz)` mean aggregation (two OSD rows for one genhz collapse to their
    mean), confirms graceful all-`NA` degradation on a simulated OSD failure, confirms the
    required-columns error message.
  - `tests/testthat/test-property-simulation.R`: `simulate_cokey_generalized()` tested both
    without a `bound_sd` column (output has none - existing contract preserved) and with one
    (output broadcasts the value onto every simulated realization for that horizon row).
- [x] Full test-file runs (`test-depth-simulation.R`: 36/36 pass; `test-property-simulation.R`:
  48/48 pass) plus a full-package `testthat::test_dir()` run (scratch `.R` file, per this repo's
  `Rscript -e` segfault workaround): **1513/1513 assertions pass, 0 failures, 0 errors** - all 15
  skips are pre-existing live-network/CRAN-gated tests unrelated to this change (confirmed
  `simulate_ssurgo_mapunit_draws()`'s own test is skipped for the same pre-existing "requires live
  SDA service" reason it always was, not a new failure).
- [x] `devtools::document()` regenerated `attach_osd_boundary_distinctness.Rd` and
  `simulate_cokey_generalized.Rd` (updated `@return`) cleanly, no warnings.
- [x] `NAMESPACE`: `attach_osd_boundary_distinctness` exported (alphabetized, between
  `assess_workflow_quality` and `backup_data`).

**Notes for future reference**: the OSD-based `bound_sd` this phase threads through is
**genhz-generic** (one value per soil series x generalized-horizon-group, from the Official Series
Description), not cokey-specific - this is the same limitation flagged in the approved plan's open
questions and left to the optional Phase 1e (`chorizon.hzdistinctness` SDA fetch) as a future
refinement, not a blocker here. Also worth noting: `simulate_and_perturb_soil_profiles()`
(`depth-simulation.R`) and `simulate_ssurgo_mapunit_draws()` (`ssurgo-simulation.R`) turned out to
be two genuinely separate pipelines in this codebase (horizon-DEPTH perturbation vs.
horizon-PROPERTY-VALUE Monte Carlo simulation) that don't currently call each other - `bound_sd`
had to be independently computed for the property-simulation side via
`attach_osd_boundary_distinctness()` rather than simply "unlocking" an existing pass-through,
since no pass-through existed between the two pipelines.

## Phase 1c - Discontinuity-gated depth kernel

**Goal**: extend the Phase 1 kernel builder to down-weight/zero the correlation between adjacent
depths whose intervening boundary is Abrupt/Clear (via `bound_sd`), falling back to the plain
kernel when no distinctness data is available.

- [x] `build_depth_correlation_kernel()` (`R/gp-modeling.R`) extended with three new parameters:
  `boundary_distinctness` (optional, same length as `depths` - `bound_sd` value for "the boundary
  immediately above this depth" in sorted-ascending order), `distinctness_range` (default
  `c(min=1, max=10)`, chosen to span `aqp::hzDistinctnessCodeToOffset()`'s actual `abrupt`
  (`1.0`) to `diffuse` (`10.0`) values for the codes `infill_missing_distinctness()` ever assigns
  - confirmed via a live `hzDistinctnessCodeToOffset(c("abrupt","clear","gradual","diffuse"))`
  call: `1.0, 2.5, 7.5, 10.0`), and `min_gate_weight` (default `0.05`, a floor so even the
  sharpest boundary leaves a small residual correlation rather than exactly zero).
- [x] **Compounds across multiple intervening boundaries, not just adjacent pairs**: implemented
  via cumulative log-sums of per-boundary gate weights (`log_cum <- c(0, cumsum(log(gate_weights)))`,
  gate between sorted depths i/j = `exp(-|log_cum[i] - log_cum[j]|)` = the product of every
  boundary weight strictly between them) - two abrupt boundaries in a row suppress the endpoints'
  correlation more than one abrupt boundary does (locked in by a dedicated test).
- [x] **Order-independent**: `depths` need not be pre-sorted - internally sorts via `order(depths)`
  for the gating pass, applies the gate matrix in sorted-index space, then maps back
  (`R[ord, ord] <- R_sorted`) to `depths`' original order. Verified by a shuffled-order regression
  test comparing against the same gating computed on already-sorted input.
- [x] Missing/NA `bound_sd` for a specific boundary degrades to gate weight `1` (no suppression
  for that boundary only) - `NULL` or all-`NA` `boundary_distinctness` skips gating entirely,
  reproducing the ungated Phase 1 kernel exactly (not merely approximately).
- [x] Re-symmetrizes and repairs via the same `ensure_positive_definite_matrix()` used throughout
  this package - no new PD-repair logic introduced.
- [x] Wrong-length `boundary_distinctness` errors immediately rather than silently misaligning
  with `depths`.

## Phase 1d - Gating regression/unit tests

Combined into the same unit of work as Phase 1c (one cohesive test block covering the whole
gating feature, added to the existing `build_depth_correlation_kernel()` test file section rather
than a separate test file).

- [x] `tests/testthat/test-gp-modeling.R`: one comprehensive test block covers -
  - An injected Abrupt boundary (`bound_sd = 1`) suppresses every depth-pair straddling it to
    `< 50%` of its plain-kernel value, while same-side pairs stay within `0.05` of the plain
    kernel (not just "near zero everywhere").
  - All-Diffuse (`bound_sd = 10`, the `distinctness_range` max -> gate weight exactly `1`)
    reproduces the plain kernel to `1e-8` tolerance - confirms gating is a strict generalization.
  - `NULL` and all-`NA` `boundary_distinctness` both reproduce the plain kernel exactly.
  - Two consecutive abrupt boundaries suppress the endpoints more than one abrupt boundary alone
    (compounding, not just per-adjacent-pair gating).
  - Wrong-length `boundary_distinctness` errors with a message naming the contract violated.
  - A shuffled-order `depths` input produces the same gating result as the already-sorted input,
    after un-shuffling for comparison (confirms the internal sort/un-sort round-trips correctly).
  - Output remains symmetric, unit-diagonal, and positive-semi-definite throughout (same checks as
    the Phase 1 plain-kernel test).
- [x] Full `test-gp-modeling.R` run: **118/118 assertions pass, 0 failures, 0 errors** (all new
  gating tests passed on the first run - no bugs found needing a fix this time, unlike Phase 1's
  `stats::dist()` dimnames gotcha).
- [x] Full-package `testthat::test_dir()` run (scratch `.R` file, per this repo's `Rscript -e`
  segfault workaround): **1529/1529 assertions pass, 0 failures, 0 errors** (16 more than Phase
  1b's 1513, matching the 16 new gating-test expectations; same 15 pre-existing live-network/CRAN
  skips, none new).
- [x] `devtools::document()` regenerated `build_depth_correlation_kernel.Rd` cleanly for the
  extended signature/parameters, no warnings.

**Notes for future reference**: the compounding-via-cumulative-log-sum trick means gating cost is
`O(n)` to build `log_cum` plus one `outer()` call - no explicit loop over depth pairs, consistent
with this package's vectorization conventions (see `PERFORMANCE_IMPROVEMENT_PLAN.md`). The
`distinctness_range`/`min_gate_weight` defaults are reasonable starting points grounded in the
actual `aqp` offset values this pipeline produces, but - like `VERTICAL_RHO`/`GLOBAL_CORRELATION_MATRIX`
in the sibling Python project - they are not yet empirically calibrated against real profile data;
Phase 5's real-AOI validation is the natural place to check whether the default suppression
strength (a `bound_sd=1` "abrupt" boundary gates to a `~0.05` floor after just two boundaries'
worth of compounding) matches actual observed cross-boundary correlation in KSSL/SSURGO profiles,
and to adjust `distinctness_range`/`min_gate_weight` if not.

## Phase 1e (optional, deferred) - `chorizon.hzdistinctness` SDA fetch

Not required for first ship of the gated kernel (Phase 1c uses OSD-derived, genhz-generic
`bound_sd`). Would add cokey-specific boundary distinctness via
`R/ssurgo-acquisition.R::execute_ssurgo_query_working()`'s `SELECT` list.

- [ ] Deferred.

## Phase 2 - Joint copula sampler

**Goal**: `sample_joint_depth_property_copula(R_depth, R_prop, n_sims)` +
`apply_copula_to_marginals()`, vectorized across realizations.

- [x] Added a new "1b. VERTICAL-CORRELATION REDESIGN: JOINT DEPTH x PROPERTY COPULA (Phase 2)"
  section to `R/multivariate-adjustment.R`, placed immediately before
  `preserve_correlation_structure()` (the function Phase 3 will provide a drop-in alternative to).
- [x] `sample_joint_depth_property_copula(R_depth, R_prop, n_sims, seed = NULL)`: draws the
  standard separable/Kronecker-MVN identity `Z = L_depth %*% Eps %*% t(L_prop)`
  (`L_depth <- t(chol(R_depth))`, `L_prop <- t(chol(R_prop))`), which achieves
  `Cov(vec(Z)) = R_prop (x) R_depth` without ever materializing the full Kronecker product matrix.
  Vectorized across all `n_sims` realizations via two array-reshape + single-matrix-multiply
  passes (no explicit per-realization loop): left-multiply by `L_depth` over the whole
  `n_depths x (k*n_sims)` flattened array at once, then `aperm()` to bring the property dimension
  first and right-multiply by `L_prop` over the whole `k x (n_depths*n_sims)` flattened array at
  once. Returns a `c(n_depths, k, n_sims)` standard-normal array.
- [x] `apply_copula_to_marginals(Z, property_matrices, gp_predictions = NULL)`: probability-
  integral-transforms each depth/property/realization cell (`pnorm()`) and inverse-transforms
  (`quantile()`) against that depth/property's own already-simulated marginal
  (`property_matrices[[prop]][depth, ]`) - a standard Gaussian copula. When `gp_predictions` is
  supplied for a property, that property's marginal is first re-centered onto the GP-predicted
  mean via a pure location shift (`curr_values + (gp_mean[i] - mean(curr_values))`) - shape/spread
  preserved exactly - **replacing** `preserve_correlation_structure()`'s sequential `gp_ratio`
  nudge with a single direct shift, per the approved plan's design.
- [x] Both exported in `NAMESPACE` (alphabetized) and documented cleanly via
  `devtools::document()` (`sample_joint_depth_property_copula.Rd`,
  `apply_copula_to_marginals.Rd`).
- [x] Unit tests added to `tests/testthat/test-multivariate-adjustment.R`:
  - **The key acceptance test**: draws a joint sample from known `R_depth`/`R_prop` targets, then
    feeds it straight into Phase 0's `validate_joint_correlation_structure()` - `overall_property_max_diff`
    and `overall_depth_max_diff` both `< 0.1` at `n_sims = 4000`, confirming the joint sampler
    actually achieves both targets simultaneously (this is the same acceptance test Phase 0 was
    built for, now exercised against the real Phase 2 sampler it was designed to judge, not just
    synthetic Kronecker-constructed data as in Phase 0's own test).
  - `n_sims < 1` errors rather than silently misbehaving.
  - `apply_copula_to_marginals()` with an identity copula (`R_depth = R_prop = diag()`) reproduces
    each depth's original marginal mean/sd within tolerance - confirms the function does not
    distort marginals on its own; any correlation comes entirely from `Z`.
  - GP-mean recentering: confirms the per-depth mean shifts to the supplied `gp_predictions`
    target while spread (sd) stays close to the original - the location-shift-only contract.
  - Dimension-contract validation errors (`property_matrices` with no names; `Z`'s property
    dimension mismatched against `length(property_matrices)`).
- [x] Full `test-multivariate-adjustment.R` run: **79/79 assertions pass, 0 failures, 0 errors**
  (all new tests passed on the first run). Full-package `testthat::test_dir()` run (scratch `.R`
  file, per this repo's `Rscript -e` segfault workaround): **1550/1550 assertions pass, 0
  failures, 0 errors** (21 more than Phase 1c/1d's 1529, matching the new test count; same
  pre-existing live-network/CRAN skips, none new).

**Notes for future reference**: this phase deliberately stops at the sampler + marginal-mapping
primitives - it does NOT yet wire into `preserve_correlation_structure()`'s callers
(`apply_gp_depth_trends()`/`process_single_cokey()`), which is Phase 3's job (a drop-in
`preserve_correlation_structure_joint()` matching the existing function's exact signature) and
Phase 4's job (gating behind a config flag so this doesn't change production behavior until
Phase 5's real-AOI validation passes). The `apply_copula_to_marginals()` GP-recentering design
(pure location shift) is simpler than the old `gp_ratio`-based multiplicative nudge - worth
double-checking during Phase 5 validation whether a purely additive shift is the right choice for
all property types (e.g. bulk density, which is strictly positive and might warrant a
multiplicative/log-scale shift instead - not yet decided either way, flagging for that phase).

## Phase 3 - `preserve_correlation_structure_joint()`

**Goal**: drop-in alternative matching `preserve_correlation_structure()`'s exact signature.

- [x] Added `preserve_correlation_structure_joint()` to `R/multivariate-adjustment.R`, immediately
  after `preserve_correlation_structure()`. Required parameters match
  `preserve_correlation_structure()` exactly, in the same order
  (`property_matrices, gp_predictions, depths, primary_property, verbose`), so it is callable
  anywhere the original is called without changing the call site - `primary_property` is accepted
  but intentionally unused (the joint method has no need for a single reference property; every
  property's correlation to every other is modeled directly via `R_prop`).
- [x] Three new OPTIONAL trailing parameters (default `NULL`/sensible values, so they never break
  a drop-in call): `gp_models`, `boundary_distinctness`, `kernel`. This is how Phase 4's config-flag
  wiring will thread the fitted GP models and OSD boundary-distinctness data through, once that
  phase updates the actual caller (`apply_gp_depth_trends()`).
- [x] **Property correlation (`R_prop`)**: estimated empirically from the surface (first) depth's
  already-drawn values - the identical computation `adjust_multivariate_depthwise_GP()`'s own
  verification step already performs (`gp-modeling.R`'s `achieved_property_correlation` check) -
  repaired via the same `ensure_positive_definite_matrix()` used throughout this package. Falls
  back to the identity matrix on a single property or an estimation failure.
- [x] **Depth correlation (`R_depth`)**: length-scale reused from `extract_depth_length_scale()`
  applied to every property in `gp_models` with a usable fit, averaged across them (Phase 1's
  "blend of everyone's fitted length-scales" design). When `gp_models` is `NULL`/empty (e.g. called
  standalone, or before Phase 4 wires the real models through), falls back to a length-scale
  spanning the full depth range (`diff(range(depths))`) - documented as a conservative default that
  keeps the function usable and testable without requiring fitted GP models.
- [x] Joint draw via `sample_joint_depth_property_copula()` + marginal mapping via
  `apply_copula_to_marginals()` (both Phase 2), each wrapped in its own `tryCatch()` that degrades
  gracefully to the original, unadjusted `property_matrices` (with a `WARN` log) on failure -
  matching `preserve_correlation_structure()`'s own graceful-failure contract exactly (same
  `n_depths < 2 || n_sims < 1` early-return check reused verbatim).
- [x] Exported in `NAMESPACE` (alphabetized, right after `preserve_correlation_structure`) and
  documented cleanly via `devtools::document()` (`preserve_correlation_structure_joint.Rd`).
- [x] Unit tests added to `tests/testthat/test-multivariate-adjustment.R` (new
  `make_independent_depth_property_matrices()` fixture helper: builds property matrices with a
  KNOWN property-correlation target but near-independent depths via per-depth Cholesky draws - the
  natural "before" state this function is meant to add vertical correlation to):
  - Drop-in signature test: calling with the exact same 4 positional arguments
    `preserve_correlation_structure()` takes works and returns correctly-shaped output; an
    unrecognized `primary_property` is silently ignored (not an error, unlike the original
    function's fallback-with-log-message behavior for the same case).
  - Insufficient-dimensions graceful degradation (`n_depths = 1` returns the original matrices
    unchanged via `expect_identical()`).
  - **The core acceptance test**: starting from data with near-zero depth correlation but a known
    property correlation, running the joint method (no `gp_models` - fallback length-scale path)
    and confirming via Phase 0's `validate_joint_correlation_structure()` that (a) property
    correlation stays within `0.15` of its original target (not distorted by adding depth
    correlation) and (b) depth-1-vs-depth-4 correlation increases from `< 0.15` (near-independent)
    to `> 0.15` for every property (real vertical structure actually added) - this is the exact
    "both structures satisfied simultaneously" property the whole redesign was built around,
    now verified end-to-end through the actual drop-in function.
  - `gp_models` supplied (a real `fit_local_gp_model_single()`-fitted model): confirms no error and
    correct output shape, exercising the `extract_depth_length_scale()` reuse path.
  - GP-mean recentering flows through `gp_predictions` correctly (per-depth mean shifts to the
    supplied target).
- [x] Full `test-multivariate-adjustment.R` run: **95/95 assertions pass, 0 failures, 0 errors**
  (all new tests passed on the first run). Full-package `testthat::test_dir()` run (scratch `.R`
  file, per this repo's `Rscript -e` segfault workaround): **1566/1566 assertions pass, 0
  failures, 0 errors** (16 more than Phase 2's 1550, matching the new test count; same
  pre-existing live-network/CRAN skips, none new).

**Notes for future reference**: this phase is still NOT wired into any pipeline call site -
`apply_gp_depth_trends()`/`process_single_cokey()` still call the original
`preserve_correlation_structure()` unconditionally. That wiring, behind a config flag defaulting
to the existing behavior, is Phase 4's job specifically so this can be A/B tested (Phase 5) before
ever affecting production output. The fallback length-scale (`diff(range(depths))`, used only when
`gp_models` isn't supplied) is deliberately conservative/generous (assumes correlation decays to
`exp(-1)` only at the full profile depth) - Phase 4's wiring should prefer passing the real
`gp_models` whenever they're available in that call chain (they already are, one level up in
`apply_local_gp_adjustments()`/`fit_local_gp_models()`) rather than relying on this fallback in
production.

## Phase 4 - Config-flag wiring

**Goal**: `config$monte_carlo$vertical_correlation_method` (`"gp_quantile_retrofit"` default vs.
`"joint_copula"`), dispatched in `apply_gp_depth_trends()`.

- [x] Added `vertical_correlation_method = "gp_quantile_retrofit"` to
  `get_monte_carlo_defaults()`'s `monte_carlo_config` list (`R/monte-carlo.R`, alongside
  `correlation_fallback`) - zero behavior change for any config built via this function, which is
  the normal path.
- [x] Added a corresponding entry to `validate_monte_carlo_config()`'s `param_specs`
  (`choices = c("gp_quantile_retrofit", "joint_copula")`) - **deliberately `required = FALSE`**
  (unlike `correlation_fallback`'s `required = TRUE`), so an older/ad-hoc config built without
  going through `get_monte_carlo_defaults()` and missing this newer key does not start failing
  validation; it's checked only when present.
- [x] `apply_gp_depth_trends()` (`R/multivariate-adjustment.R`) gained two new OPTIONAL trailing
  parameters: `config = NULL` and `gp_models = NULL`. At the point that previously called
  `preserve_correlation_structure()` unconditionally, now reads
  `config$monte_carlo$vertical_correlation_method %||% "gp_quantile_retrofit"` and dispatches:
  - `"gp_quantile_retrofit"` (or missing/`NULL` config - the default): calls
    `preserve_correlation_structure()` exactly as before, byte-for-byte unchanged code path.
  - `"joint_copula"`: calls `preserve_correlation_structure_joint()` (Phase 3), threading
    `gp_models` through directly, and extracting a per-unique-depth `boundary_distinctness` vector
    from `valid_data$bound_sd` when that column is present (Phase 1b's threaded-through OSD
    distinctness data) - the first non-NA `bound_sd` value at each unique depth, since
    `simulate_cokey_generalized()` broadcasts the same value across every realization at a given
    depth. `NULL` boundary_distinctness (no gating) when the column isn't present at all.
- [x] Unit tests added:
  - `tests/testthat/test-monte-carlo.R`: confirms `get_monte_carlo_defaults()` defaults to
    `"gp_quantile_retrofit"`, and `validate_monte_carlo_config()` accepts both valid choices,
    rejects an invalid one, and still validates cleanly when the key is absent entirely (locking
    in the `required = FALSE` contract).
  - `tests/testthat/test-multivariate-adjustment.R`: confirms `config = NULL`, `config = list()`,
    and an explicit `vertical_correlation_method = "gp_quantile_retrofit"` all produce
    `identical()` output at a fixed seed (the actual "zero behavior change by default" claim,
    verified bit-for-bit, not just "no error"); confirms `"joint_copula"` dispatches to a
    genuinely different code path (same seed, same inputs, `clay_pct` output is NOT
    `all.equal()` to the retrofit method's own output) while still returning valid,
    correctly-shaped output; confirms `bound_sd`-column extraction and `gp_models` pass-through
    both run end-to-end without error under `"joint_copula"`.
- [x] Full `test-multivariate-adjustment.R` (102/102) and `test-monte-carlo.R` (114/114) runs, all
  passing on the first try. Full-package `testthat::test_dir()` run (scratch `.R` file, per this
  repo's `Rscript -e` segfault workaround): **1578/1578 assertions pass, 0 failures, 0 errors** (12
  more than Phase 3's 1566, matching the new test count; same pre-existing live-network/CRAN
  skips, none new).
- [x] `devtools::document()` regenerated `apply_gp_depth_trends.Rd` cleanly for the two new
  parameters, no warnings.

**Notes for future reference / explicit scope boundary**: this phase wires the dispatch INTO
`apply_gp_depth_trends()` itself, exactly as the approved plan specified - it does NOT thread
`config`/`gp_models` further up the call chain, from the top-level `integrate_monte_carlo_with_gp()`
API down through `process_single_cokey()` -> `apply_local_gp_adjustments()`/
`apply_nrcs_trend_adjustments()` -> `apply_local_depth_trends()`/`apply_gp_depth_trends()`. Those
intermediate functions still call `apply_gp_depth_trends()` without passing `config`/`gp_models`,
so they always get the default `"gp_quantile_retrofit"` behavior regardless of what a top-level
caller's config says. This was a deliberate scope decision, not an oversight: the approved plan's
Phase 4 text says only "dispatched in `apply_gp_depth_trends()`," and threading two new parameters
through 4+ intermediate function signatures is meaningfully more surface area than "add a dispatch
point." **Phase 5's A/B validation can proceed without this** - it calls `apply_gp_depth_trends()`
directly (as this phase's own tests already do) rather than through the full top-level API. If
Phase 5 validates the joint method for production use, wiring the full top-level call chain (so an
end user can opt in via `simulate_ssurgo_mapunit_draws()`'s own config, not just by calling
`apply_gp_depth_trends()` directly) should be scoped as its own follow-up unit of work.

## Phase 5 - Validation on real data

**Goal**: run Phase 0's validator over representative real-AOI data, old vs. new method side by
side; flip the default only once the new method matches or beats the old method's fidelity at
acceptable runtime cost.

Implemented as a new "Step 9" section in `vignettes/gp-modeling-multivariate-adjustment.Rmd`,
rather than a live Salinas Valley SDA run (the original plan text's suggestion) - this vignette
already carries real, cached Amador County, CA SSURGO data (`inst/extdata/ssurgo_amador.rds`) and
real fitted GP models (`clay_gp_model`/`sand_gp_model` from the vignette's own Step 3), giving
offline-reproducible validation against real data without a live network dependency - a stronger
fit for this package's documentation/CI story than a live-network-only validation script would
have been.

- [x] Added Step 9 to the vignette: builds synthetic Monte Carlo-style realizations with a known
  clay-sand property correlation (`target_rho = -0.65`, matching Step 6's own established value)
  but genuinely INDEPENDENT draws at each of 5 depths (0-80cm) - unlike Step 6's own
  `clay_mat`/`sand_mat` fixture, which replicates one shared draw across every depth by
  construction (useful for the reference-quantile demo in Step 6, but already perfectly
  depth-correlated, so not a valid "before" state for testing whether a method adds vertical
  correlation).
- [x] Runs `apply_gp_depth_trends()` twice on identical input: once with the production default
  (`config = NULL`) and once with `vertical_correlation_method = "joint_copula"` plus the real
  `gp_models`, then uses Phase 0's `validate_joint_correlation_structure()` to report achieved
  property correlation (every depth) and achieved depth correlation (every property) for both.
- [x] **Real finding, verified numerically** (`clay_gp_model`'s own `extract_depth_length_scale()`
  = `12.8` cm on the real Amador data): the CURRENT PRODUCTION method
  (`gp_quantile_retrofit`/`preserve_correlation_structure()`) does not merely "fail to add" vertical
  correlation - it actively induces SPURIOUS near-perfect correlation (`> 0.99` at every depth
  tested, 0-80cm) between clay content at any two depths, regardless of physical distance. This
  is a direct, previously-undocumented consequence of its reference-quantile mechanism: every
  realization is pinned to ONE FIXED percentile (set once at the surface) and forced to occupy
  that same percentile at every subsequent depth, with no dependence on the GP model's own fitted
  smoothness. The `joint_copula` method, by contrast, reproduces a decay tied directly to the real
  fitted length-scale (`12.8` cm) - correlation near `1` at the surface, ~`0.15` at 20cm, and
  essentially zero by 40cm+ - matching an exponential kernel with that length-scale, and matching
  the physical behavior a vertical-correlation model should have.
  - Both methods preserved the clay-sand PROPERTY correlation equally well
    (`property_corr_max_diff` small for both - `~0.018` retrofit, `~0.044` joint, at `n_sims=400`)
    - confirming the joint method's core promise (both correlation axes satisfied simultaneously)
    holds on real data, not just the synthetic Kronecker-constructed fixtures used in
    Phases 0-3's own unit tests.
- [x] Vignette rendered end-to-end successfully via `rmarkdown::render()` (scratch `.R` file, per
  this repo's `Rscript -e` segfault workaround) against the real cached data, twice (once to
  discover the actual numeric result, once after correcting the vignette's own narrative text to
  match what was actually measured rather than an initial incorrect guess - see note below).
- [x] Full-package `testthat::test_dir()` still **1578/1578 assertions pass, 0 failures, 0 errors**
  (this phase touched only the vignette `.Rmd`, no `R/` source changes).

**Notes for future reference**: the first draft of this vignette section guessed the production
method would simply leave depth correlation near its (near-zero) starting point, reasoning that
"a mean-nudge shouldn't affect residual correlation." Running it against real data immediately
falsified that guess - the actual behavior is more interesting and more concerning: the
reference-quantile mechanism induces artificial full-profile rank-lock, not "no vertical
correlation modeling." The vignette text was corrected to report the measured result rather than
the original hypothesis - a good illustration of why this phase existed as a mandatory
verification step even though the code changes it exercises (Phases 0-4) were already fully unit
tested. **This finding meaningfully strengthens the case for the redesign**: the production
default isn't merely missing a nice-to-have feature, it's producing a correlation structure with no
grounding in the GP model's own fitted physical behavior. Two follow-ups worth flagging for anyone
picking this up next: (1) the discontinuity-gating defaults (`distinctness_range`/
`min_gate_weight`, Phase 1c) are still not calibrated against real KSSL/SSURGO profile lag
correlations - this phase validated the CORE joint-sampling mechanism, not the boundary-gating
feature specifically; (2) `vertical_correlation_method` is still only reachable by calling
`apply_gp_depth_trends()` directly (Phase 4's documented scope boundary) - not yet through the
top-level `integrate_monte_carlo_with_gp()`/`simulate_ssurgo_mapunit_draws()` APIs a real user
would actually call.

## Phase 6 - Regression tests

**Goal**: new-path fidelity tests (reusing Phase 0/2 synthetic-target pattern) + explicit
bit-identical-old-path-at-default-config test.

Most of this goal was already satisfied incrementally, as each of Phases 0-4 added its own
regression tests along the way (see those phases' own entries above) - specifically:
- **New-path fidelity** (reusing the Phase 0/2 synthetic-target pattern): Phase 2's
  `sample_joint_depth_property_copula()`/`apply_copula_to_marginals()` tests, and Phase 3's
  `preserve_correlation_structure_joint()` "induces real vertical correlation while preserving the
  existing property correlation" test, both already feed real function output into Phase 0's
  `validate_joint_correlation_structure()` against known targets.
- **Bit-identical-old-path-at-default-config**: Phase 4's
  `apply_gp_depth_trends()` test already asserts `identical()` output across `config = NULL`,
  `config = list()`, and an explicit `"gp_quantile_retrofit"`.

What Phase 6 added as genuinely new work, once that assessment was done:

- [x] **Closed a gap Phase 4 had explicitly deferred**: Phase 4 wired the
  `vertical_correlation_method` dispatch INTO `apply_gp_depth_trends()` itself, but noted that
  `config`/`gp_models` weren't threaded any further UP the call chain - so `"joint_copula"` was
  only reachable by calling `apply_gp_depth_trends()` directly, not through the real top-level API
  a user would actually call. Investigating this while writing Phase 6's tests found it was a
  smaller gap than Phase 4's note assumed: `apply_local_gp_adjustments()`
  (`R/multivariate-adjustment.R`) already receives `config` AND already fits real `local_gp_models`
  in scope (via `fit_local_gp_models()`) - both were simply going unused for this purpose. Only
  `apply_local_depth_trends()` (the one intermediate hop) needed new optional `config`/`gp_models`
  parameters threaded through to its own `apply_gp_depth_trends()` call, plus one line at
  `apply_local_gp_adjustments()`'s existing `apply_local_depth_trends()` call site passing
  `config = config, gp_models = local_gp_models` (previously it passed neither). Two small,
  low-risk, purely-additive edits - not the "own follow-up unit of work" Phase 4 anticipated.
  `apply_nrcs_trend_adjustments()` (the OTHER top-level path, "regional" GP models) was NOT
  similarly wired - `process_single_cokey()` doesn't thread `config` to it at all today, a
  pre-existing gap unrelated to this redesign; left as a documented follow-up, not addressed here.
- [x] Unit test added (`tests/testthat/test-multivariate-adjustment.R`): confirms
  `"joint_copula"` is reachable end-to-end from `apply_local_gp_adjustments()` (same seed/inputs,
  `retrofit_config` vs. `joint_config` produce genuinely different, non-`all.equal()` output - the
  dispatch actually threads all the way through, not silently dropped at either hop), AND that
  `apply_local_gp_adjustments()`'s own `config = NULL` default still reproduces its production
  default `identical()`ly - the same "zero behavior change unless opted in" guarantee Phase 4
  already established one level down, now confirmed one level higher too.
- [x] Full-package `testthat::test_dir()` run (scratch `.R` file, per this repo's `Rscript -e`
  segfault workaround): **1582/1582 assertions pass, 0 failures, 0 errors** (4 more than Phase 5's
  1578 - the vignette phase added no `testthat` assertions, only this phase's new test file
  changes did). `devtools::document()` regenerated `apply_local_depth_trends.Rd` cleanly for the
  two new parameters, no warnings.

**Notes for future reference**: `apply_nrcs_trend_adjustments()`'s config gap (noted above) means
a user relying on `integration_method = "nrcs_gp"` (rather than the default `"hybrid"`, which
applies both NRCS and local adjustments per cokey - so still benefits from this phase's fix via its
local-GP half) cannot yet reach `"joint_copula"` at all, even by setting the config flag. Worth a
small follow-up (`process_single_cokey()`'s NRCS branch would need the same two-parameter threading
this phase just did for the local branch) if that integration method needs the same capability.

## Phase 7 - Benchmark

**Goal**: cost/quality tradeoff table (old vs. new) at representative depths x properties x n_sims
scale, logged here in `PERFORMANCE_IMPROVEMENT_PLAN.md`-style table format.

- [x] Benchmarked `preserve_correlation_structure()` (production default) vs.
  `preserve_correlation_structure_joint()` (with REAL fitted `gp_models` supplied, so the "reuse
  the GP's own length-scale" path is exercised, not the fallback) across 5 scenarios anchored to
  this repo's own established representative scale
  (`data-raw/benchmark_performance.R::benchmark_merge_adjusted_data()`'s 10 depths x 1000
  realizations, and `simulate_ssurgo_mapunit_draws()`'s own `n_mc = 1000` default), 5 repeated
  calls per scenario, `system.time()`:

| Scenario | `preserve_correlation_structure()` | `preserve_correlation_structure_joint()` | Ratio (joint/retrofit) |
|---|---|---|---|
| Representative (10 depths, 6 props, 1000 sims) | 0.136s/call | 0.234s/call | 1.72x |
| All properties (10 depths, 9 props, 1000 sims) | 0.190s/call | 0.282s/call | 1.48x |
| Heavier (10 depths, 6 props, 2000 sims) | 0.244s/call | 0.428s/call | 1.75x |
| Many depths (20 depths, 6 props, 1000 sims) | 0.278s/call | 0.396s/call | 1.42x |
| Small (4 depths, 2 props, 500 sims) | 0.010s/call | 0.012s/call | 1.20x |

- [x] **Honest result: the joint method is genuinely, consistently slower** - roughly 1.2x-1.75x
  the retrofit method's wall-clock across every scale tested, not a performance win. This is
  expected given what it does per call that the retrofit method doesn't: estimate + PD-repair
  `R_prop` via `cor()`/`eigen()`, build + PD-repair `R_depth` via `build_depth_correlation_kernel()`
  (its own `chol()`+`ensure_positive_definite_matrix()`), draw the joint sample via two
  `chol()` factorizations plus array reshape/matrix-multiply passes
  (`sample_joint_depth_property_copula()`), then map to marginals via a per-depth `quantile()`
  loop (`apply_copula_to_marginals()`) - genuinely more numerical work than the retrofit's single
  sequential `gp_ratio`-nudge loop.
  - Still well within a workable per-call budget for production use (well under half a second per
    cokey even at the heaviest scale tested), so this is a real but non-prohibitive cost - not the
    kind of "should be a win, wasn't" result `PERFORMANCE_IMPROVEMENT_PLAN.md` documented for
    `fuse_texture_group()` (this phase never claimed a speedup was the goal - the design goal was
    correctness/fidelity, benchmarked here purely so that tradeoff is documented honestly rather
    than assumed free).
  - No optimization work was undertaken here (out of this phase's scope, which was to MEASURE the
    tradeoff, not chase it down) - a natural next step if `"joint_copula"` moves toward becoming
    the default would be an `Rprof()` pass (following this package's own established methodology)
    to find which of the four sub-steps above dominates at real-AOI scale, the same way
    `PERFORMANCE_IMPROVEMENT_PLAN.md`'s own Tier 1-4 items were diagnosed.
- [x] Benchmark script kept informally (scratch `.R` file, not committed to `data-raw/` - unlike
  `benchmark_performance.R`'s own committed, maintainer-run scenarios, this phase's script was a
  one-off measurement pass; promoting it to a permanent `data-raw/` benchmark would be a reasonable
  follow-up if `"joint_copula"` sees continued development).

---

## Open questions / risks (carried from the approved plan)

- **Single shared depth length-scale across properties** is a simplification vs. a full Linear
  Model of Coregionalization - simpler to implement, strictly more faithful than today's
  single-primary-property retrofit, but not the most general possible model.
- **OSD-derived distinctness is genhz-generic, not cokey-specific** - a cokey whose real horizons
  deviate from its OSD-typical boundary sharpness won't be gated correctly until the optional
  Phase 1e SDA-fetch lands.
- **Performance**: measured in Phase 7 - `preserve_correlation_structure_joint()` is consistently
  1.2x-1.75x slower than the production default across representative scales (never faster; this
  was a correctness/fidelity-motivated redesign, not a performance one), but stays well under
  half a second per cokey even at the heaviest scale tested - a real, now-documented cost rather
  than an assumed-free one.

## Overall status: all 7 core phases (0-7, including 1b/1c/1d) complete

The joint depth x property Gaussian-copula method exists, is fully unit-tested (1582/1582
assertions passing package-wide), is reachable both directly (`apply_gp_depth_trends()`) and
through the real top-level local-GP API (`apply_local_gp_adjustments()`/`process_single_cokey()`),
defaults to zero behavior change unless explicitly opted into via
`config$monte_carlo$vertical_correlation_method = "joint_copula"`, and has been validated against
real AOI data (Amador County SSURGO) via the package's own vignette - which surfaced a genuinely
important finding: the CURRENT production method induces spurious near-perfect (`>0.99`)
correlation across the entire depth profile regardless of physical distance, while the new method
reproduces a decay tied to the property's own fitted spatial smoothness.

That finding is a strong argument that the CURRENT default has a real, previously-undocumented
defect - but it is not, by itself, a green light to flip the default. Flipping a production
default silently changes the output distribution of every downstream AWC/PIW90/uncertainty number
this package produces, for every existing caller who hasn't opted into anything. That is a
consequential, hard-to-reverse-after-the-fact change (once users have built workflows, reports, or
downstream models on top of a given output distribution, un-flipping it later is its own
disruption). The decisions below are what should be resolved - and by whom - before that flip,
roughly in the order they'd block a maintainer from saying yes.

## Decisions required before flipping the production default

### 1. Is single-AOI, single-soil-type validation sufficient evidence?

**What's been done**: Phase 5 validated the core joint-sampling mechanism on ONE AOI (Amador
County, CA), using ONE soil group's fitted GP models (`clay_pct`/`sand_pct`, `xv_group`), with
SYNTHETIC Monte Carlo realizations (real correlation targets, but not real simulated property
values run through the actual `simulate_ssurgo_mapunit_draws()` pipeline end to end).

**Why it matters**: the Amador finding (spurious `>0.99` retrofit correlation) is compelling
precisely because it's a structural property of the retrofit algorithm, not a property of Amador
County's soil data specifically - so there's good reason to expect it generalizes. But "good
reason to expect" is not the same as "verified across the range of conditions this package
actually runs on": different soil orders (the vignette's own texture-correlation examples already
note residual/stratified/alluvial soils behave differently), different depth ranges, different
`sim_comppct` scales, and - important given Phase 1c/1d - profiles with real lithologic
discontinuities (buried horizons, stone lines) where the NEW method's discontinuity gating is the
untested piece, not the already-validated core copula mechanism.

**What would resolve this**: run the same before/after comparison
(`validate_joint_correlation_structure()` old vs. new) across a handful of AOIs spanning different
soil orders/regions already used elsewhere in this package's test/vignette data (the Salinas Valley
AOI `PERFORMANCE_IMPROVEMENT_PLAN.md` already established as a benchmark reference is a natural
second AOI - it wasn't used here specifically because the vignette's cached, offline-reproducible
Amador data was available and the Salinas Valley setup requires a live SDA call per that
benchmark's own documented usage), and ideally at least one AOI/cokey known to have genuine
lithologic discontinuities, to exercise Phase 1c's gating path for the first time on real
(not synthetic-injected) boundary distinctness.

**Decision owner**: whoever owns soilSIM's output-correctness sign-off (the pattern in
`PERFORMANCE_IMPROVEMENT_PLAN.md` suggests this has historically been the user directing this
work, benchmarking against a real Salinas Valley run before accepting changes).

**Risk if skipped**: the default could flip on the strength of one AOI's evidence and then produce
an unexpected regression (over- or under-correlation) for soil types/regions not represented in
that one validation pass - discovered downstream, in production, rather than here.

### 2. Is the discontinuity-gating feature (Phase 1c/1d) ready, or should the default flip exclude it?

**What's been done**: the gating mechanism itself is fully unit-tested with synthetic
`bound_sd` values (Phase 1d) and mathematically sound (compounds correctly across multiple
boundaries, degrades gracefully to the ungated kernel). What has NOT been done: calibrating
`distinctness_range`/`min_gate_weight`'s default numeric values against real KSSL/SSURGO profile
lag-correlations - they're grounded in `aqp::hzDistinctnessCodeToOffset()`'s literature-derived
offset scale (abrupt=1.0 ... diffuse=10.0), not in an empirical "how much should an abrupt boundary
actually suppress correlation" fit. This mirrors the exact same gap the sibling Python project's
own `VERTICAL_RHO` documentation flagged for ITS defaults ("literature-informed estimate, not
empirically calibrated from SSURGO data").

**Why it matters**: an uncalibrated gating strength could either (a) under-suppress, leaving some
of the retrofit's spurious-correlation problem intact across real discontinuities, or (b)
over-suppress, artificially decorrelating depths that are more continuous than the generic
OSD-genhz-level `bound_sd` estimate suggests for that specific cokey (see decision #3 below - this
compounds with the genhz-generic-vs-cokey-specific gap).

**What would resolve this**: a real calibration pass - pull multi-horizon KSSL profiles with known
boundary distinctness codes, compute empirical lag-1 correlation across labeled abrupt/clear/
gradual/diffuse boundaries (the same style of analysis the sibling Python project's
`scripts/estimate_rho.py` does for its own `VERTICAL_RHO`), and check whether the current
`distinctness_range = c(min=1, max=10)` / `min_gate_weight = 0.05` defaults produce gating strength
in the right ballpark, or need adjustment.

**Decision owner / risk**: same as #1 - an uncalibrated gating default shipped as part of a
flipped production default risks a second, independent source of miscalibrated output alongside
whatever the core copula mechanism itself introduces, making any downstream discrepancy harder to
diagnose (two unvalidated knobs moving at once instead of one).

**A narrower option worth considering**: flip the default for the CORE joint-copula mechanism (the
part Phase 5 actually validated) while leaving discontinuity gating off by default
(`boundary_distinctness = NULL`, already this function's own no-gating default) until #2's
calibration is separately complete. This decouples two decisions that don't need to be made
together.

### 3. Genhz-generic vs. cokey-specific boundary distinctness (Phase 1e)

**What's been done**: Phase 1b threads OSD-derived `bound_sd` through the pipeline, but it's
**one value per soil-series x generalized-horizon-group** (from the Official Series Description),
broadcast to every cokey of that series/genhz combination - not derived from that specific cokey's
own actual horizon boundaries.

**Why it matters**: a cokey whose real horizons deviate from its OSD-typical boundary sharpness
(e.g. a locally disturbed profile, or one with a discontinuity the series' "typical" description
doesn't capture) gets gated using the wrong distinctness value. This is a real, understood
limitation, not a bug - but it means the gating feature is only as good as how representative the
OSD series description is for a given cokey.

**What would resolve this**: Phase 1e, deferred and optional per the original plan - adding
`chorizon.hzdistinctness` to `execute_ssurgo_query_working()`'s SDA fetch would give real,
cokey-specific SSURGO boundary-distinctness codes instead of the OSD-genhz-generic proxy.

**Decision owner**: whoever is deciding the scope of a default-flip vs. accepting Phase 1e as a
"good enough for now, documented limitation" - this is a judgment call about how much the
genhz-generic approximation actually matters in practice (it may be a small effect for the
majority of gently-varying profiles and only matter for the minority with genuine local
discontinuities the series description doesn't capture) versus how much implementation effort a
cokey-specific fetch is worth before shipping.

**Risk if skipped**: lower than #1/#2 - this degrades gating QUALITY for a subset of cokeys, it
doesn't introduce a systematic bias, and the gating feature already degrades gracefully to
"no gating" when data is missing. Reasonable to explicitly defer/accept as a documented limitation
rather than a blocker, if the maintainer agrees.

### 4. `apply_nrcs_trend_adjustments()` config-threading gap

**What's been done**: Phase 6 wired `config`/`gp_models` through the LOCAL GP path
(`apply_local_gp_adjustments()` -> `apply_local_depth_trends()` -> `apply_gp_depth_trends()`).
The REGIONAL/NRCS GP path (`apply_nrcs_trend_adjustments()`) was NOT similarly wired -
`process_single_cokey()` doesn't thread `config` to it at all today (a pre-existing gap unrelated
to this redesign, not introduced by it).

**Why it matters**: under the default `integration_method = "hybrid"` (applies both NRCS and local
adjustments per cokey), a flipped default would still take effect via the local-GP half - but a
caller using `integration_method = "nrcs_gp"` specifically would see NO change at all regardless of
the config flag, silently. That's a confusing inconsistency if the default flips: two different
integration methods would disagree about which vertical-correlation method is "the default."

**What would resolve this**: the same two-parameter threading pattern Phase 6 already used for the
local path, applied to `process_single_cokey()`'s NRCS branch and
`apply_nrcs_trend_adjustments()`/`get_nrcs_gp_predictions()`'s GP-model access (the "regional" GP
models are already available in that code path, similar to how `local_gp_models` were already
available for the local path - likely comparable effort to Phase 6's fix).

**Decision owner / risk**: whoever decides whether `"nrcs_gp"`-only callers are a use case that
matters for the default-flip decision, or whether it's acceptable to close this gap in a follow-up
after the local-path default already flips (asymmetric coverage during a transition period is a
legitimate interim state, as long as it's documented, not a hard blocker).

### 5. Performance: is 1.2x-1.75x slower acceptable, or does it need optimization first?

**What's been done**: Phase 7 measured (not optimized) the cost - consistently slower than the
retrofit method across every scale tested, staying under ~0.5s/call even at the heaviest scale
benchmarked (10-20 depths, 6-9 properties, 1000-2000 realizations).

**Why it matters**: this cost is paid PER COKEY, and a full AOI run processes many cokeys (the
Salinas Valley small-AOI benchmark in `PERFORMANCE_IMPROVEMENT_PLAN.md` used 15 mukeys as a
reference scale) - the per-call overhead compounds across an AOI. Whether ~1.5x is acceptable
depends on what the current end-to-end AOI runtime budget is and how much headroom exists in it;
this plan did not re-benchmark a full AOI run end to end with the new method active, only the
isolated `preserve_correlation_structure()`/`preserve_correlation_structure_joint()` call.

**What would resolve this**:
- A full-AOI timing comparison (old vs. new method) using this package's own established
  `data-raw/benchmark_performance.R` methodology, to see the REAL end-to-end cost, not just the
  isolated per-call cost.
- If that reveals a genuine problem, an `Rprof()` profiling pass on the four sub-steps Phase 7
  identified (property-correlation estimation + PD repair, depth-kernel construction + PD repair,
  the joint draw itself, and the per-depth `quantile()`-based marginal mapping) to find which one
  dominates and is worth optimizing - following the same diagnose-before-optimizing discipline
  `PERFORMANCE_IMPROVEMENT_PLAN.md`'s own Tier 1-4 items used (several of which found the
  originally-assumed bottleneck was wrong once actually profiled).

**Decision owner / risk**: whoever owns AOI-run SLAs/turnaround expectations. Risk if skipped is
lower-severity than #1/#2 (a slower-but-correct result is recoverable by later optimization,
unlike a wrong result), but worth resolving before a flip if AOI runtime is a hard constraint for
any consumer of this package.

### 6. Two design choices flagged during implementation, not yet revisited

Two simplifications were made explicitly and flagged for revisiting once real evidence existed -
worth a maintainer decision now that Phase 5 has real numbers, even though neither blocks
correctness of what's shipped:

- **GP-mean recentering is a pure additive location shift** (Phase 2's design note): reasonable for
  properties like clay/sand percentage, but properties that are strictly positive and often
  right-skewed (e.g. bulk density) might warrant a multiplicative or log-scale shift instead, to
  avoid an additive shift pushing simulated values toward/past a physically implausible boundary
  (e.g. negative bulk density) at depths far from the GP training data. Not yet checked against
  real output for a non-texture property.
- **Single shared depth length-scale, blended across all properties in a cokey** (vs. a full Linear
  Model of Coregionalization where each property keeps its own vertical smoothness while still
  being jointly correlated): simpler and strictly more faithful than the old single-primary-property
  retrofit, but not the most general possible model. Worth a decision on whether this
  simplification is acceptable as a permanent design choice or should be revisited if a future
  validation pass (decision #1) finds it produces a poor fit for properties whose GP models have
  very different length-scales from each other within the same cokey.

### Summary: minimum bar for a default flip

At minimum, before flipping `vertical_correlation_method`'s default package-wide: (1) validate on
more than one AOI/soil type, including at least one profile with real discontinuities; (2) either
calibrate the discontinuity-gating defaults or ship with gating off by default as a decoupled
follow-up; (4) decide how to handle the `"nrcs_gp"`-only integration-method gap (fix it, or
document the interim asymmetry); (5) confirm full-AOI runtime overhead is acceptable for real
consumers. (3) and (6) are lower-severity and reasonable to accept as documented limitations for
an initial flip, revisited later if evidence suggests they matter in practice.

---

## Decision-resolution phases (8-13)

Source plan (approved): `C:\Users\jmaynard\.claude\plans\evaluate-these-files-from-quirky-pretzel.md`.
These phases resolve the six default-flip decision points above and execute the flip itself.
Placed at the end of this document, after the phases (0-7) and decision analysis they build on.

| Phase | Description | Status |
|---|---|---|
| 8 | Decouple core method from discontinuity gating (`vertical_correlation_gating` flag) | [x] Done |
| 9 | Multi-AOI validation (Salinas Valley + discontinuity profile) | [x] Done |
| 10 | Full-AOI performance benchmark | [x] Done |
| 11 | Thread `config`/`gp_models` through the NRCS GP path | [x] Done |
| 12 | Empirical check on GP-recentering for strictly-positive properties | [x] Done - no action needed |
| 13 | Execute the default flip | [x] Done |

## Phase 8 - Decouple core method from discontinuity gating

**Goal**: let the already-validated core copula mechanism flip independently of discontinuity
gating (not yet calibrated) via a new, separate config flag - resolves decision #2's blocking risk
without requiring calibration first.

- [x] Added `vertical_correlation_gating = FALSE` to `get_monte_carlo_defaults()`'s
  `monte_carlo_config` list (`R/monte-carlo.R`), alongside `vertical_correlation_method`. Zero
  behavior change for any existing config.
- [x] Added a corresponding `validate_monte_carlo_config()` entry (`type = "logical"`,
  `required = FALSE` - same non-breaking rationale already used for `vertical_correlation_method`).
- [x] `apply_gp_depth_trends()` (`R/multivariate-adjustment.R`) now reads
  `isTRUE(config$monte_carlo$vertical_correlation_gating)` and only extracts/passes
  `boundary_distinctness` to `preserve_correlation_structure_joint()` when that flag is `TRUE` -
  previously, `boundary_distinctness` was extracted automatically whenever a `bound_sd` column was
  present, which (since `attach_osd_boundary_distinctness()` runs unconditionally in
  `simulate_ssurgo_mapunit_draws()`) meant there was no way to use `"joint_copula"` WITHOUT gating
  whenever OSD lookup succeeded. Now gating requires its own explicit opt-in, independent of the
  core method choice.
- [x] `@param config` roxygen docs updated on `apply_gp_depth_trends()` to document the new
  sub-key relationship.
- [x] Unit tests added:
  - `tests/testthat/test-monte-carlo.R`: confirms the new key defaults to `FALSE`, validates as a
    logical (rejects a non-logical value), and validates cleanly when absent (non-required
    contract).
  - `tests/testthat/test-multivariate-adjustment.R`: confirms gating stays OFF under
    `"joint_copula"` when the flag is unset (identical output to explicitly `FALSE`, and identical
    - modulo the `bound_sd` passthrough column itself, which is unrelated to whether gating is
    active - to the case where `bound_sd` isn't present in `cokey_data` at all) even with a real
    injected discontinuity (`bound_sd = 1`, abrupt) in the data; confirms explicitly setting the
    flag `TRUE` actually changes the output when that same discontinuity is present (the flag does
    something when turned on).
  - **Test bug found and fixed during this phase**: an initial `expect_identical()` between the
    gating-off-with-bound_sd-present and no-bound_sd-at-all cases failed - not because gating was
    incorrectly active, but because the `bound_sd` column itself (untouched passthrough metadata,
    unrelated to the adjustment) was present in one output and absent in the other. Fixed by
    comparing only the columns both outputs share.
- [x] Full `test-multivariate-adjustment.R` (109/109) and `test-monte-carlo.R` (119/119) runs, all
  passing. Full-package `testthat::test_dir()` run (scratch `.R` file, per this repo's
  `Rscript -e` segfault workaround): **1590/1590 assertions pass, 0 failures, 0 errors** (8 more
  than Phase 7's 1582, matching the new test count).
- [x] `devtools::document()` regenerated `apply_gp_depth_trends.Rd` cleanly, no warnings.

**Notes for future reference**: this phase makes NO claim about what the right gating default
SHOULD be once calibrated (decision #2 remains genuinely open) - it only ensures the core-method
decision and the gating-strength decision can be made on separate timelines. Phase 9's
discontinuity-profile check (below) explicitly set `vertical_correlation_gating = TRUE` to
actually exercise this path on real (not synthetic-injected) data for the first time.

## Phase 9 - Multi-AOI validation (Salinas Valley)

**Goal**: extend Phase 5's single-AOI (Amador County) validation to a second, genuinely different
AOI/soil type, and exercise Phase 1c/8's discontinuity gating on real (not synthetic-injected)
`bound_sd` data for the first time - resolves decision #1.

- [x] Used the package's own established Salinas Valley reference area
  (`data-raw/benchmark_performance.R::small_aoi()`'s ~200m x 200m benchmark box was tried first
  and found too sparse - 58 rows, 11 cokeys, 0 adequate GP-fitting groups, all-NA `bound_sd` -
  widened to a ~3km x 3km box around the same center point:
  `POLYGON((-121.665 36.598, -121.635 36.598, -121.635 36.622, -121.665 36.622, -121.665 36.598))`).
  Live `download_and_prepare_ssurgo()`/SDA call, same pattern the vignette's own Step 0 uses for
  Amador County, run as a scratch script (not a vignette chunk, per the plan) since this AOI
  requires live network access rather than cached data.
- [x] Fetched 271 real horizon rows across 83 cokeys; `prepare_nrcs_training_data()` +
  `build_stratified_gp_models()` (properties `clay_pct`/`sand_pct`, `min_profiles_per_group = 3`,
  `min_observations_per_group = 15` - same thresholds the vignette's Amador validation uses) fit 4
  real GP models across 2 adequate soil-series groups (`Placentia`, `Mocho`).
- [x] Ran the same old-vs-new `apply_gp_depth_trends()` comparison Phase 5 used (independent-
  per-depth synthetic realizations with a known `target_rho = -0.65` clay-sand correlation, real
  fitted `clay_gp_model`/`sand_gp_model` for the `Placentia` group) and checked via
  `validate_joint_correlation_structure()`:

| Method | Property corr max diff | Clay surface-vs-deepest corr |
|---|---|---|
| `gp_quantile_retrofit` | 0.0034 | **0.9938** |
| `joint_copula` (gating off) | 0.0559 | **0.0140** |

  - **Phase 5's core finding replicates on a second, independent AOI/soil type**: the retrofit
    method again produces near-perfect (`>0.99`) correlation across the entire depth profile
    (`0.997, 0.997, 0.998, 0.994` at depths 20/40/60/80cm vs. the surface), regardless of physical
    distance - confirming this is a structural property of the retrofit algorithm, not an
    Amador-County-specific artifact.
  - **New finding this AOI surfaced that Amador's smoother-fitting group didn't**: `Placentia`'s
    real fitted GP models came out with an extremely short length-scale
    (`extract_depth_length_scale()` = `0.0082` cm, vs. Amador's `12.8` cm) - meaning this
    particular group's own training data is fit as essentially having NO smooth depth trend (very
    noisy clay/sand-vs-depth relationship in the underlying field data for this series). The joint
    method, which faithfully reuses the GP's own fitted smoothness, correctly reproduces
    near-zero correlation at every depth in this case (`0.011, -0.010, 0.055, 0.014`) - arguably
    still more honest than the retrofit's spurious full-profile lock-in, but a genuinely different
    regime than Amador's moderate, visually-legible decay. Property correlation was preserved by
    both methods either way (`property_corr_max_diff` small for both).
- [x] **Real production bug found and fixed via this validation**:
  `attach_osd_boundary_distinctness()` (Phase 1b) initially returned all-`NA` `bound_sd` for every
  real compname/genhz combination in this AOI, despite `query_osd_distinctness()` itself
  successfully returning real, varied distinctness data (confirmed via a direct diagnostic call -
  e.g. real `Placentia` OSD data shows `abrupt`/`clear`/`gradual` boundaries at different depths,
  bound_sd 1.0-7.5). Root cause: `soilDB::fetchOSD()` returns its `id` column UPPER-CASED
  regardless of requested case (a quirk already documented in this file's own
  `fetch_osd_horizons_cached()` docs, which notes the ORIGINAL caller,
  `simulate_and_perturb_soil_profiles()`, was unaffected because it joins by `genhz` only, never
  by `id`) - but `attach_osd_boundary_distinctness()` DOES join by `(compname, genhz)`, so it
  silently inherited a case-sensitivity bug the original code was specifically documented as
  immune to. Not caught by Phase 1b's own unit tests because their mocked fixtures happened to use
  matching case throughout (`"amador"`/`"amador"`), never exercising the mismatch.
  - **Fix**: `attach_osd_boundary_distinctness()` (`R/depth-simulation.R`) now joins on
    `toupper(compname)` against `toupper(id)`, matching the case-insensitive pattern
    `fetch_osd_horizons_cached()` already uses elsewhere in this same file for the identical
    reason.
  - **Regression test strengthened** (`tests/testthat/test-depth-simulation.R`): the existing
    mocked-`query_osd_distinctness()` test now mocks upper-cased `id` values (`"AMADOR"`/
    `"PENTZ"`, matching what a real `soilDB::fetchOSD()` call actually returns) against
    lower-cased `compname` in `hz_data` - the exact mismatch that silently broke before the fix,
    now locked in as the test's own realistic fixture rather than an artificially matching one.
  - Full `test-depth-simulation.R` re-run after the fix: **37/37 assertions pass, 0 failures, 0
    errors**. Full-package `testthat::test_dir()` re-run: still 0 failures/errors (see below).
- [x] Re-ran the Salinas Valley comparison after the fix - real `bound_sd` values now resolve
  correctly (e.g. `Placentia` genhz A = 1.5, B = 5.5; genhz-averaged from real per-horizon OSD
  distinctness codes). Ran `apply_gp_depth_trends()` with
  `vertical_correlation_gating = TRUE` and a real per-depth `bound_sd` sequence
  (`1.5, 5.5, 5.83, 7.5, 7.5`, derived from this AOI's own observed genhz/distinctness data) - the
  call completed successfully with no error, producing clay depth-correlation values of
  `-0.020, -0.040, -0.057, -0.001`. Because the ungated baseline for this particular group was
  ALREADY near-zero (per the short-length-scale finding above), gating's own suppression effect
  wasn't visually distinguishable from that already-low baseline for this specific group -
  gating's synthetic-data behavior (Phase 1d) is unaffected by this and remains the primary
  evidence that the mechanism itself works correctly; what's still missing is a real-data
  demonstration on a soil group with BOTH a real discontinuity AND enough underlying smoothness
  for gating's suppression to be visually distinguishable from a near-zero baseline - not achieved
  in this pass, worth a targeted follow-up (e.g. explicitly search for a real profile with a
  larger fitted length-scale AND a documented abrupt/clear boundary) if a cleaner illustration is
  needed before a final decision on gating defaults.
- [x] Full-package `testthat::test_dir()` run after the bugfix (scratch `.R` file, per this repo's
  `Rscript -e` segfault workaround): **1591/1591 assertions pass, 0 failures, 0 errors** (1 more
  than Phase 8's 1590, matching the one new assertion added to the strengthened regression test).

**Notes for future reference**: this phase's most valuable outcome was arguably the bug it found,
not the numbers it produced - a stark demonstration of why Phase 9 (real multi-AOI validation) was
correctly identified as a hard blocker rather than something to skip on the strength of Phase 5's
single-AOI result alone: `attach_osd_boundary_distinctness()` had been "working" in every unit test
and in the Amador vignette (which never happened to need real cross-case matching) but was
silently non-functional against real production data from a second AOI. Decision #3
(genhz-generic vs. cokey-specific distinctness) is unaffected by this fix - it's a bug in USING the
genhz-generic data correctly, not a change to what data is available.

## Phase 10 - Full-AOI performance benchmark

**Goal**: time a full `simulate_ssurgo_mapunit_draws()` run over a real AOI, old vs. new method,
not just the isolated `preserve_correlation_structure()`/`preserve_correlation_structure_joint()`
call Phase 7 benchmarked - resolves decision #5.

- [x] **Prerequisite gap closed**: `simulate_ssurgo_mapunit_draws()` had no `config` parameter at
  all, and its internal call chain (`maybe_adjust_soil_data_depth_trend()` ->
  `adjust_one_cokey_depth_trend()` -> `apply_local_gp_adjustments()`) dropped it even where it
  existed - meaning `"joint_copula"` was unreachable through the actual top-level production API
  a real user calls, only through `apply_local_gp_adjustments()` directly (Phase 6's fix) or
  `apply_gp_depth_trends()` directly (Phase 4's original wiring). Added `config = NULL` to all
  three functions (`R/ssurgo-simulation.R`), threaded straight down to the existing
  `apply_local_gp_adjustments(..., config = config)` call `adjust_one_cokey_depth_trend()` already
  had in scope - the same low-risk, purely-additive pattern Phase 6 used one level down.
  Regression test added (`tests/testthat/test-ssurgo-simulation.R`): confirms
  `maybe_adjust_soil_data_depth_trend()` reaches `"joint_copula"` end-to-end and that
  `config = NULL` still reproduces the exact production default.
- [x] Benchmarked `simulate_ssurgo_mapunit_draws()` end-to-end (live SDA call) over the same
  widened Salinas Valley AOI Phase 9 used (83 cokeys, `n_mc = 1000`), once per method:

| Method | Elapsed | `nrow(draws)` |
|---|---|---|
| `gp_quantile_retrofit` | 32.35s | 23659 |
| `joint_copula` (gating off) | 31.03s | 23663 |

  **Ratio: 0.96x - joint_copula was essentially on par with (marginally faster than) the
  production default at real full-AOI scale**, a materially different picture than Phase 7's
  isolated-call benchmark (1.2x-1.75x slower). Likely explanation: `apply_local_gp_adjustments()`'s
  own docs already note GP FITTING itself (`fit_local_gp_models()`) is the dominant cost of the
  whole per-cokey depth-trend-adjustment step at real-AOI scale - the
  `preserve_correlation_structure()`/`preserve_correlation_structure_joint()` call Phase 7 isolated
  is a much smaller fraction of that total once real GP fitting is included, diluting the
  measured-in-isolation overhead into noise. **Conclusion: performance is not a blocker for a
  default flip** - no `Rprof()` optimization pass is needed per the plan's own "only if it's a
  meaningful fraction" branch condition.
- [x] **Real bug found and fixed via this benchmark** (the actual point of running it against real
  data rather than trusting Phase 7's synthetic numbers): the FIRST full-AOI run crashed
  `preserve_correlation_structure_joint()`'s copula-to-marginal mapping step with "missing value
  where TRUE/FALSE needed" for multiple cokeys, silently caught by that function's own `tryCatch()`
  and logged as `WARN`-level "Copula-to-marginal mapping failed" (correctly avoiding a full
  pipeline crash, but silently skipping the adjustment for every affected cokey/property - a real
  robustness gap, and one that also made the FIRST run's timing numbers untrustworthy, since a
  meaningful number of cokeys were short-circuiting via the error path rather than doing the real
  computation).
  - **Root cause**: `apply_copula_to_marginals()` (`R/multivariate-adjustment.R`) guarded its
    GP-recentering and quantile-mapping steps with `if (stats::var(x, na.rm = TRUE) > 0)` -
    `stats::var()` returns `NA` (not `FALSE`) when `x` is entirely `NA`, and `if (NA > 0)` errors
    rather than falling through to the intended "no variation - keep original values" else branch.
    An entirely-`NA` property column at a given depth is a real, non-rare case on messy field data
    (the same rows `simulate_cokey_generalized()`'s own per-row texture `tryCatch()` already skips
    when a texture triplet is incomplete) - no synthetic unit-test fixture across Phases 0-9 had
    ever happened to include one.
  - **Fix**: both `if` conditions wrapped in `isTRUE()` (`isTRUE(stats::var(x, na.rm = TRUE) > 0)`),
    so an `NA` result from an all-`NA` input is treated as `FALSE` (no variation), matching the
    function's already-documented intent.
  - **Regression test added** (`tests/testthat/test-multivariate-adjustment.R`): an all-`NA` depth
    row for one property, confirmed `apply_copula_to_marginals()` no longer errors, the all-`NA`
    depth stays all-`NA` (not fabricated), unaffected depths are untouched, and the same case is
    also confirmed at the `preserve_correlation_structure_joint()` call-site level (where the real
    bug actually surfaced) rather than only at the lower-level function.
  - Re-ran the full-AOI benchmark after the fix: **zero errors** (confirmed via log grep), and the
    corrected timing numbers above are what's reported (the pre-fix run's numbers, which showed an
    even larger apparent "speedup" for `joint_copula`, are NOT reported here since they were an
    artifact of the crash-and-skip behavior, not real computation).
- [x] Full-package `testthat::test_dir()` run after the fix (scratch `.R` file, per this repo's
  `Rscript -e` segfault workaround): **1600/1600 assertions pass, 0 failures, 0 errors** (9 more
  than Phase 9's 1591, matching the new test count from this phase's prerequisite-gap fix and the
  `apply_copula_to_marginals()` bugfix regression test).

**Notes for future reference**: like Phase 9, this phase's most valuable output was arguably the
bug it found rather than the timing number itself - a second consecutive real-data validation
phase that caught a production-breaking gap no synthetic unit test had ever exercised. This
reinforces the same lesson Phase 9 already surfaced: real-AOI validation (Phases 9-10) earns its
place as a hard blocker before any default flip, not a nice-to-have on top of Phase 5's
single-AOI/synthetic-fixture result.

## Phase 11 - Thread `config`/`gp_models` through the NRCS GP path

**Goal**: mirror Phase 6's local-GP-path fix for the REGIONAL/NRCS path, so
`integration_method = "nrcs_gp"` callers aren't silently stuck on the old method regardless of the
config flag once it's the default - resolves decision #4.

- [x] `apply_nrcs_trend_adjustments()` (`R/multivariate-adjustment.R`) gained a `config = NULL`
  parameter, threaded through to its own `apply_gp_depth_trends()` call (previously called
  unconditionally with no `config`/`gp_models`, unlike the local-GP path's already-fixed
  `apply_local_depth_trends()`).
- [x] Also extracts the actual fitted NRCS GP model OBJECTS (not just their predictions) for the
  joint method's depth-kernel length-scale reuse - the same
  `gp_models[[nrcs_prop]]$models[[model_group]]` lookup `get_nrcs_gp_predictions()` already
  performs internally, re-keyed here by `prop` (the `cokey_data`/`property_matrices` column name)
  to match `apply_gp_depth_trends()`'s `gp_models` contract, instead of `nrcs_prop` (the NRCS
  property-mapping name) `get_nrcs_gp_predictions()` itself uses internally.
- [x] `process_single_cokey()` now passes its own `config` parameter to the
  `apply_nrcs_trend_adjustments()` call site - previously omitted entirely even though `config` was
  already in scope and already passed to the sibling local-GP branch two lines below (the same
  "gap smaller than expected once actually investigated" pattern Phase 6 found for the local path).
- [x] `@param config`/`@param gp_models` roxygen docs added.
- [x] Unit test added (`tests/testthat/test-multivariate-adjustment.R`): fits real NRCS-shaped GP
  models for two mappable properties (`claytotal`/`sandtotal` -> `clay_pct`/`sand_pct` via
  `get_nrcs_property_mapping()`) - two properties specifically, since
  `apply_gp_depth_trends()`'s method dispatch only branches with `>= 2` properties; a first-draft
  single-property version of this test passed trivially without actually exercising the dispatch
  (both methods fall through the same "individual adjustment" path with only one property) and had
  to be corrected once the test failure exposed that. Confirms `"joint_copula"` produces genuinely
  different output than the retrofit at the `apply_nrcs_trend_adjustments()` level, confirms
  `config = NULL` still reproduces the exact production default, and confirms the same reachable
  end-to-end through `process_single_cokey()` itself (the actual call site the gap existed in).
- [x] Full `test-multivariate-adjustment.R` run: **123/123 assertions pass, 0 failures, 0 errors**
  (one bug in the test itself, not the source fix, was caught and corrected - see above). Full-package
  `testthat::test_dir()` run (scratch `.R` file, per this repo's `Rscript -e` segfault
  workaround): **1607/1607 assertions pass, 0 failures, 0 errors** (7 more than Phase 10's 1600,
  matching the new test count).
- [x] `devtools::document()` regenerated `apply_nrcs_trend_adjustments.Rd` cleanly, no warnings.

**Notes for future reference**: with Phases 6 and 11 both done, `"joint_copula"` is now reachable
end-to-end through BOTH top-level integration paths (`use_local_gp` and `use_nrcs_gp`) and through
`simulate_ssurgo_mapunit_draws()`'s own top-level API (Phase 10's prerequisite fix) - decision #4
is fully resolved, not just documented as an accepted gap. The default `integration_method =
"hybrid"` applies both branches per cokey regardless, so this phase's practical impact is mainly
for callers who explicitly restrict themselves to `"nrcs_gp"` only.

## Phase 12 - Empirical check on GP-recentering for strictly-positive properties

**Goal**: check whether `apply_copula_to_marginals()`'s additive GP-mean recentering (Phase 2's
design choice, flagged as unrevisited) ever pushes a strictly-positive property (e.g. bulk
density) to a non-physical value, before deciding whether a multiplicative/log-scale alternative
is actually needed - resolves the first half of decision #6, empirically rather than
speculatively.

- [x] Reused the same Salinas Valley AOI data Phases 9/10 already fetched (live SDA call). Fit a
  REAL `bulk_density` GP model (`build_stratified_gp_models()`, `Placentia` soil-series group,
  training data: depths 0-106cm, bulk density 1.62-1.745 g/cm3).
- [x] Predicted the GP trend at a depth 100cm BEYOND the training data's own range (206cm total) -
  deliberately the case where the additive shift is largest, per the plan's own framing (`"far
  from the training data"`). `predict_gp_depth_trends()` returned `1.63` g/cm3 there - GPfit's own
  extrapolation behavior kept the predicted mean close to the training range rather than
  diverging, which is itself informative (GP extrapolation tends toward stability, not runaway
  drift, for a reasonably-fit model).
- [x] Ran `apply_copula_to_marginals()`'s real GP-recentering step directly (2000 synthetic
  realizations, realistic ~0.15 g/cm3 marginal spread, recentered onto the far-depth GP mean):
  achieved mean `1.627` g/cm3 (matching the `1.63` target closely), full output range
  `[1.18, 2.09]` g/cm3.
  - **No non-physical values**: zero realizations `<= 0`; zero realizations outside a generous
    `[0.3, 2.5]` g/cm3 sanity bound; `assess_realistic_values()` - this package's own established
    plausible-range check, already used by `validate_gp_models()`, with its own `[0.3, 3.0]`
    g/cm3 bound for `bulk_density` - returned `TRUE` (plausible) for every recentered value.
- [x] **Conclusion: no real problem found - closing this decision as "no action needed"**, per the
  plan's own "if no real problem is found: document and close, avoids unnecessary complexity"
  branch. No `recentering_mode` parameter added.

**Notes for future reference**: this check used one real AOI/soil-group/property combination and
one "far" depth (100cm beyond training). It does not prove the additive shift is safe in every
conceivable case (e.g. a property with a much steeper/more divergent fitted trend, or a
depth-window request far more extreme than 100cm beyond training) - but combined with GPfit's own
observed tendency to keep extrapolated means near the training range rather than diverging
unboundedly, there's no concrete evidence this design choice needs revisiting. If a future
real-AOI run (or production use) surfaces an actual non-physical output for a strictly-positive
property, revisit via the `recentering_mode = c("additive", "multiplicative")` design sketched in
Phase 2's original notes - not built preemptively here.

## Phase 13 - Execute the default flip

**Goal**: flip `vertical_correlation_method`'s default from `"gp_quantile_retrofit"` to
`"joint_copula"` now that Phases 0-12 resolved every blocking decision point, update every test
that asserted the old default, and update the vignette's narrative - `"gp_quantile_retrofit"`
remains fully supported as an explicit opt-out throughout.

- [x] **The flip itself**: `get_monte_carlo_defaults()`'s `monte_carlo_config` list
  (`R/monte-carlo.R`) now sets `vertical_correlation_method = "joint_copula"` (was
  `"gp_quantile_retrofit"`). `vertical_correlation_gating` stays `FALSE` - a deliberately separate,
  still-open decision (Phase 8's whole point), not flipped alongside the core method.
- [x] **A second, more consequential fallback also flipped for consistency**: `apply_gp_depth_trends()`
  (`R/multivariate-adjustment.R`) has its OWN internal `%||%` fallback for when `config` is `NULL`
  entirely (not routed through `get_monte_carlo_defaults()`) - this was separately hardcoded to
  `"gp_quantile_retrofit"` and would otherwise have left "no config passed" meaning something
  DIFFERENT depending on whether a caller went through `get_monte_carlo_defaults()` first. Flipped
  in lockstep so "no config" means the same thing everywhere in the package, matching every
  regression test across Phases 4/6/8/10/11 (all of which call with `config = NULL` to test
  "default behavior").
- [x] Roxygen docs updated on all four functions whose `@param config` describes the default
  (`apply_gp_depth_trends()`, `apply_nrcs_trend_adjustments()`, `apply_local_depth_trends()`,
  `maybe_adjust_soil_data_depth_trend()`, `simulate_ssurgo_mapunit_draws()`) to describe
  `"joint_copula"` as current and `"gp_quantile_retrofit"` as the explicit opt-out.
  `devtools::document()` regenerated all affected `.Rd` files cleanly.
- [x] **Every test asserting the old default as a baseline was found and updated** (not just
  patched to pass - each one re-verified to still test something meaningful post-flip):
  - `tests/testthat/test-monte-carlo.R`: `get_monte_carlo_defaults()`'s default assertion flipped
    to `"joint_copula"`; explicit `"gp_quantile_retrofit"` opt-out still validates cleanly.
  - `tests/testthat/test-multivariate-adjustment.R`: `apply_gp_depth_trends()`'s "defaults to
    gp_quantile_retrofit" test rewritten as "defaults to joint_copula - config NULL/empty/explicit
    all agree"; a NEW separate test locks in that the `"gp_quantile_retrofit"` opt-out still
    reaches the original, unmodified `preserve_correlation_structure()` code path and still
    produces genuinely different output than the new default. The Phase 6 (local-GP,
    `apply_local_gp_adjustments()`) and Phase 11 (NRCS-GP, `apply_nrcs_trend_adjustments()`)
    end-to-end tests both had a `config = NULL` assertion pinned to the old default - both flipped
    to assert against the `joint_result`/new default instead, with the explicit-opt-out comparison
    preserved unchanged.
  - `tests/testthat/test-ssurgo-simulation.R`: Phase 10's `maybe_adjust_soil_data_depth_trend()`
    test's `config = NULL` assertion flipped the same way.
- [x] **Two real bugs found while fixing the tests above** (both pre-existing gaps, surfaced only
  because flipping the default meant `config = NULL` test paths finally exercised code that used
  to only run under an explicit opt-in):
  - Phase 11's NRCS-path test and Phase 10's `maybe_adjust_soil_data_depth_trend()` test had both
    used only ONE mappable/adjustable property - `apply_gp_depth_trends()`'s method dispatch only
    branches with `>= 2` properties, so both tests were passing "for the wrong reason" (both
    methods trivially fell through the same single-property fallback path and were never actually
    exercising the divergence being tested). Fixed both fixtures to use two properties, matching
    the pattern already used elsewhere.
  - Phase 10's `maybe_adjust_soil_data_depth_trend()` test used `simulation_number = 1` (a single
    realization per depth) - with `n_sims = 1`, quantile-based remapping is a no-op for EITHER
    algorithm (nothing to remap against), so both methods degenerately produced the same
    unchanged output regardless of which was used, and `preserve_correlation_structure_joint()`'s
    own `R_prop` estimation step additionally errored on the single-realization input
    (`sapply()` over 2 properties each contributing a length-1 vector returns a plain vector, not
    a matrix, and `stats::cor()` requires a matrix - already caught gracefully by the function's
    own `tryCatch()`, logging a `WARN` and falling back to an identity `R_prop` exactly as
    designed, so this was NOT a crash, just a genuinely degenerate test case). Fixed by using a
    realistic multi-realization fixture (`n_sims = 30`) instead of patching the source - the
    graceful degradation for single-realization data is correct, intentional behavior, not a bug.
- [x] **Vignette updated** (`vignettes/gp-modeling-multivariate-adjustment.Rmd`): Step 9's
  narrative reframed from "a newer, opt-in method" to "why the production default changed,"
  keeping the exact same real measured numbers (unaffected by the flip, since that section's own
  code always explicitly set `config` for both methods being compared - never relied on the
  ambient default). Step 7's description updated to stop calling `preserve_correlation_structure()`
  "the production path" (no longer accurate). **One real vignette bug caught and fixed**: the
  `retrofit_result <- apply_gp_depth_trends(...)` call in Step 9 had relied on the (now-flipped)
  default with no explicit `config` - after the flip this would have silently computed the
  `joint_copula` result while a variable named `retrofit_result` and all the surrounding
  commentary described it as the retrofit method, invalidating the whole comparison. Fixed by
  making both calls' `config` explicit, matching the already-explicit `joint_result` call.
  Re-rendered end-to-end (`rmarkdown::render()`, scratch `.R` file per this repo's `Rscript -e`
  segfault workaround) - success, all 41 chunks processed.
- [x] Full-package `testthat::test_dir()` run (scratch `.R` file, same workaround): **1608/1608
  assertions pass, 0 failures, 0 errors** (1 more than Phase 12's 1607 - Phase 12 added no test
  assertions of its own, only this phase's test-fixture corrections did).
- [x] Full `devtools::check()` run (scratch `.R` file, ~16 minutes): **0 errors, 0 warnings, 3
  NOTEs**. Two pre-existing/environmental (`unable to verify current time`; the untracked
  `pkgdown/` directory, already noted as pre-existing in earlier phases' own check runs). The
  THIRD was real and new: `attach_osd_boundary_distinctness(): no visible binding for global
  variable 'compname_upper'` - a `dplyr::mutate()`-created NSE column from Phase 9's
  case-insensitive-join fix, not a genuine undefined global. Fixed by adding `"compname_upper"` to
  `R/soilSIM-package.R`'s existing `utils::globalVariables()` declaration - the same established
  pattern already used ~60 times in that same file for identical dplyr-NSE false positives, not a
  new convention. Verified via `devtools::load_all()` + a full `test-depth-simulation.R` re-run
  (37/37 passing) rather than a second full 16-minute `devtools::check()` pass, given the fix's
  well-established low-risk nature in this codebase.

**Notes for future reference**: this phase is a strong illustration of why "update the tests to
match the new behavior" is a genuinely dangerous instruction if taken literally - several of the
tests that needed updating had a SECOND, independent problem (single-property or single-realization
fixtures that never really exercised what they claimed to) that a blind find-and-replace of
`"gp_quantile_retrofit"` -> `"joint_copula"` would have silently preserved as broken tests that
merely pass. Each one was re-verified to test something real post-flip, not just edited to stop
failing.

---

## Final status: all 13 phases complete - `joint_copula` is now the package default

`config$monte_carlo$vertical_correlation_method` defaults to `"joint_copula"` as of this work.
`"gp_quantile_retrofit"` (the original algorithm) remains fully supported as an explicit opt-out,
not deprecated or removed - every function in the dispatch chain still accepts and correctly
routes to it. Every decision point raised in this document's own "Decisions required before
flipping the production default" section was resolved, not merely documented as accepted risk:

| # | Decision | Resolution |
|---|---|---|
| 1 | Single-AOI validation sufficiency | Resolved - Phase 9 validated a second, independent AOI (Salinas Valley), replicating the core finding |
| 2 | Discontinuity-gating readiness | Resolved by decoupling (Phase 8) - core method flipped; gating stays a separate, still-`FALSE`-by-default decision, not blocking |
| 3 | Genhz-generic vs. cokey-specific distinctness | Accepted as a documented limitation, per the original risk assessment (Phase 1e still deferred/optional) |
| 4 | `apply_nrcs_trend_adjustments()` config gap | Resolved - Phase 11 closed it, mirroring Phase 6's local-path fix |
| 5 | Performance | Resolved - Phase 10's full-AOI benchmark found no meaningful overhead at real scale (0.96x - within noise of the retrofit method) |
| 6 | GP-recentering / length-scale-blending design choices | GP-recentering resolved empirically (Phase 12 - no real problem found); length-scale blending remains an accepted simplification per the original risk assessment |

**Real bugs found and fixed along the way** (the concrete payoff of insisting on real-AOI
validation rather than trusting synthetic-fixture test coverage alone): a case-sensitivity bug in
`attach_osd_boundary_distinctness()` that silently broke real boundary-distinctness lookups
(Phase 9); a crash in `apply_copula_to_marginals()` on all-`NA` property columns, a real and
non-rare condition on messy field data (Phase 10); a dispatch-never-branches gap in two different
tests that were passing without actually exercising what they claimed to (Phase 13); a stale
reliance on the (now-flipped) default inside the vignette's own comparison code, which would have
silently invalidated its own headline result (Phase 13); and a genuine `R CMD check` NOTE from an
unlabeled dplyr NSE column (Phase 13).

Total package-wide test suite: **1608/1608 assertions pass, 0 failures, 0 errors**.
`devtools::check()`: **0 errors, 0 warnings, 3 NOTEs** (all accounted for - 2 pre-existing/
environmental, 1 fixed).
