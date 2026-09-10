# Multi-Source Raster Fusion Pipeline

## Overview

This pipeline is soilSIM's core-vs-adapter architecture made concrete. A generic **fusion core**
(`R/core-distributions-raster.R`, `R/core-fusion.R`, `R/cache.R`) implements
distribution fitting, Bayesian prior/likelihood fusion, and disk caching purely in terms of
`terra::SpatRaster` percentile-value inputs — it has no knowledge of SSURGO, SOLUS, or any other
specific data source, and will fuse *any* prior/likelihood `SpatRaster` pair sharing a common
probability grid. Two **adapters** currently plug into that core over a shared area of interest
(AOI): `R/adapter-ssurgo-simulate.R` supplies the **prior** side (SSURGO map units, Monte Carlo-simulated
to per-cell percentiles), and `R/adapter-solus.R` supplies the **likelihood** side (SOLUS100's
own published low/prediction/high percentiles, fetched directly). Neither adapter is special-cased
by the core: both simply produce a `list(values = <named list of percentile-value SpatRasters>,
probs = <matching probabilities>)` structure, which is all `fuse_property_adaptive()` requires from
either side. This is deliberate: a future data source (e.g. HWSD, SoilGrids) could plug in as a new
adapter supplying either the prior or the likelihood side — or both, fusing two entirely new
sources — without any change to the fitting/fusion/caching core itself.

## Core Functions

### Generic core (core-distributions-raster.R, core-fusion.R, cache.R)

#### 1. `fit_normal_raster()` / `quantile_normal_raster()`
**Purpose**: Fit a Normal distribution's `mu`/`sigma` rasters from three percentile-value rasters
(low/median/high), then evaluate the fitted distribution at any target quantile — pure raster
arithmetic, EXACT relative to the closed-form Normal quantile formula.

**Parameters**:
```r
fit_normal_raster(p_lo_r, p50_r, p_hi_r, p_lo, p_hi)
# p_lo_r, p50_r, p_hi_r - SpatRasters for the low/median/high percentile values
# p_lo, p_hi            - the probabilities p_lo_r/p_hi_r represent (e.g. 0.025/0.975)

quantile_normal_raster(fit, q)
# fit - output of fit_normal_raster()
# q   - target quantile probability
```
**Returns**: `fit_normal_raster()` returns `list(mu = SpatRaster, sigma = SpatRaster)`;
`quantile_normal_raster()` returns a single `SpatRaster` of quantile values.

**Behavior**: `mu` is simply the median raster; `sigma` is `(p_hi_r - p_lo_r) / (qnorm(p_hi) -
qnorm(p_lo))`, i.e. the value-per-z-unit spacing implied by the two tail percentiles. Evaluation is
`mu + sigma * qnorm(q)`. Both are single arithmetic expressions across the whole raster, so cost is
independent of `q` or of how many times a fitted `fit` is re-queried.

#### 2. `quantile_linear_cdf_raster()`
**Purpose**: Evaluate a piecewise-linear inverse-CDF at a fixed quantile `q`, interpolating between
whichever two percentile-value rasters bracket `q`. Serves as the nonparametric fallback target for
every other fit routine in this file.

**Parameters**:
```r
quantile_linear_cdf_raster(value_rasters, probs, q)
# value_rasters - list of SpatRasters, sorted the same way as probs
# probs         - numeric probabilities matching value_rasters' order
# q             - target quantile probability (must satisfy min(probs) <= q <= max(probs))
```
**Returns**: A single `SpatRaster` of interpolated values.

**Behavior**: Because the percentile breakpoints (`probs`) are identical at every cell, the
bracketing pair of indices and the scalar interpolation weight are resolved once in plain R (no
per-cell branching), then only `value_rasters[[i]] + weight * (value_rasters[[i+1]] -
value_rasters[[i]])` runs as raster arithmetic.

#### 3. `fit_beta_mom_raster()` / `quantile_beta_mom_raster()`
**Purpose**: Method-of-moments Beta fit — APPROXIMATE relative to MLE, used only as the
Newton-Raphson seed in `fit_beta_mle_newton_raster()` below (not intended as a standalone fit).

**Parameters**:
```r
fit_beta_mom_raster(mean_r, var_r)
# mean_r, var_r - SpatRasters of the rescaled-to-unit-interval mean/variance

quantile_beta_mom_raster(fit, q)
# fit - output of fit_beta_mom_raster()
# q   - target quantile probability
```
**Returns**: `fit_beta_mom_raster()` returns `list(alpha = SpatRaster, beta = SpatRaster)`;
`quantile_beta_mom_raster()` returns a `SpatRaster` (via `terra::lapp()` calling `stats::qbeta()`
per cell, since `qbeta()` itself has no raster-native form).

**Behavior**: `common <- mean*(1-mean)/var - 1`; `alpha = mean*common`, `beta = (1-mean)*common` —
the standard closed-form Beta method-of-moments identity, applied elementwise.

#### 4. `fit_beta_mle_newton_raster()` / `quantile_beta_mle_newton_raster()`
**Purpose**: Fit a Beta distribution's `(alpha, beta)` via vectorized Newton-Raphson MLE across an
entire raster at once — `fitdistrplus::fitdist()`'s generic `optim()`-based MLE has no vectorized
form, so this exploits the Beta log-likelihood's known closed-form score/Hessian in terms of
`digamma()`/`trigamma()` (both vectorized base-R functions) to update all cells' estimates
simultaneously per iteration. Validated upstream against `fitdistrplus::fitdist()` to ~1e-3–1e-5.

**Parameters**:
```r
fit_beta_mle_newton_raster(value_rasters, bounds, n_iter = 15, eps = 1e-6)
# value_rasters - list of percentile-value SpatRasters (any count >= 3)
# bounds        - c(lower, upper) physical bounds
# n_iter        - fixed Newton-Raphson iteration count (validated sufficient at 15)
# eps           - clamp epsilon away from the exact 0/1 rescaled boundary

quantile_beta_mle_newton_raster(fit, q)
# fit - output of fit_beta_mle_newton_raster()
# q   - target quantile probability
```
**Returns**: `list(alpha = SpatRaster, beta = SpatRaster)`; the quantile function returns a
`SpatRaster` (again via `terra::lapp()` + `qbeta()`).

**Behavior**: Values are rescaled to `[eps, 1-eps]` within `bounds` (clamping cells that land
exactly on 0/1, where `log(x)`/`log(1-x)` would be `-Inf`). The Beta method-of-moments fit seeds
`(alpha, beta)`, clamped to a `0.1` floor. Each of `n_iter` iterations computes digamma/trigamma
rasters for `alpha`, `beta`, and `alpha+beta`, forms the 2x2 score vector (`S1`, `S2`) and Hessian
(`H11`, `H12`, `H21`, `H22`), solves the 2x2 linear system in closed form (`det`, `d_alpha`,
`d_beta`), and updates `alpha_r`/`beta_r` with a `0.01` floor. Two internal (non-exported) helpers,
`raster_digamma()`/`raster_trigamma()`, wrap `terra::app()` around base R's `digamma`/`trigamma`.

#### 5. `fit_metalog_linear_raster()` / `quantile_metalog_linear_raster()`
**Purpose**: Fit a metalog distribution via an exactly-determined linear solve, vectorized as
raster arithmetic — the raster counterpart of `R/core-distributions.R`'s `fit_metalog_linear()`. EXACT
when feasible (implied density non-negative everywhere); matches `rmetalog::metalog()` to ~1e-12, but reproduces none of `rmetalog`'s LP-based feasibility-correction
fallback for when it isn't feasible — see `check_metalog_feasibility_raster()` below.

**Parameters**:
```r
fit_metalog_linear_raster(value_rasters, probs, bounds, boundedness)
# value_rasters - list of INTERIOR percentile-value rasters (excluding p=0/p=1 if present)
# probs         - matching interior probabilities
# bounds        - c(lower, upper); required unless boundedness = "u"
# boundedness   - one of "u"/"sl"/"su"/"b" (see core-distributions.R's metalog_to_z())

quantile_metalog_linear_raster(fit, q, bounds, boundedness)
# fit                  - output of fit_metalog_linear_raster()
# q                    - target quantile probability
# bounds, boundedness  - same as fit_metalog_linear_raster()'s arguments
```
**Returns**: `list(a = <list of coefficient SpatRasters>, term = <number of metalog terms>)`; the
quantile function returns a single `SpatRaster`.

**Behavior**: The number of interior percentiles must equal the number of metalog terms, making the
fit an exactly-determined linear system rather than an optimization. The fixed, non-spatial
probability-grid basis matrix `Y` (from `R/core-distributions.R`'s existing, reused
`metalog_basis_matrix()`) and its inverse `Y_inv` are computed once in plain R; each percentile
value raster is transformed to z-space via the existing, reused `metalog_to_z()`; each coefficient
raster `a[[j]]` is then `Reduce(+, Y_inv[j,i] * z_rasters[[i]])` across `i`. Quantile evaluation
builds the basis row at `q`, recombines it with the coefficient rasters, and inverts back via the
existing, reused `metalog_from_z()`.

#### 6. `check_metalog_feasibility_raster()`
**Purpose**: Vectorized per-cell feasibility probe for a metalog fit — a valid metalog quantile
function must be monotonically increasing (equivalent to non-negative density everywhere); this
flags cells where it isn't.

**Parameters**:
```r
check_metalog_feasibility_raster(fit, bounds, boundedness, y_grid = seq(0.02, 0.98, by = 0.02))
# fit                 - output of fit_metalog_linear_raster()
# bounds, boundedness - same as fit_metalog_linear_raster()'s arguments for this fit
# y_grid              - probability grid to probe
```
**Returns**: A boolean `SpatRaster`, `TRUE` where infeasible.

**Behavior**: Probes `quantile_metalog_linear_raster()` at each successive `y_grid` value and flags
cells where the probe value *decreases* from the previous probe — a probe, not a proof, but fully
vectorized. Streams one probe raster at a time (rather than materializing the whole grid at once),
set to avoid an allocation failure at large cell counts.

#### 7. `quantile_metalog_linear_with_fallback()`
**Purpose**: Metalog quantile evaluation with automatic per-cell fallback to
`quantile_linear_cdf_raster()` on infeasible cells — verified to have zero effect on
feasible cells and an exact `linear_cdf` match on infeasible ones.

**Parameters**:
```r
quantile_metalog_linear_with_fallback(fit, infeasible_r, full_value_rasters, full_probs,
                                       q, bounds, boundedness)
# fit                              - output of fit_metalog_linear_raster() (fit on INTERIOR percentiles)
# infeasible_r                     - output of check_metalog_feasibility_raster() for this same fit
# full_value_rasters, full_probs   - the FULL percentile set (including p=0/p=1 if available)
# q                                - target quantile probability
# bounds, boundedness              - same as fit_metalog_linear_raster()'s arguments
```
**Returns**: `list(value = <blended quantile SpatRaster>, used_fallback = infeasible_r verbatim)`.

**Behavior**: Feasibility (`infeasible_r`) is computed once per fit and passed in, since it depends
only on the fitted coefficients, not on which `q` is evaluated. `terra::ifel()` picks the metalog
quantile where feasible and the linear-CDF fallback (computed against the *full* percentile set,
not just the interior ones used to fit the metalog) where infeasible.

#### 8. `fit_gamma_mom_raster()`
**Purpose**: Method-of-moments Gamma fit from a list of percentile-value rasters — a thin raster
wrapper around `R/core-fusion.R`'s existing scalar `moments_to_gamma()`.

**Parameters**:
```r
fit_gamma_mom_raster(value_rasters)
# value_rasters - list of percentile-value SpatRasters
```
**Returns**: `list(shape = SpatRaster, rate = SpatRaster)`.

**Behavior**: Computes `mean_r`/`var_r` across the list via `Reduce(+, ...)` (the only genuinely
raster-specific part), then delegates to `moments_to_gamma()`, which is pure elementwise arithmetic
and therefore already works unchanged on `SpatRaster` inputs.

#### 9. `align_percentile_probs()`
**Purpose**: Re-express one side's percentile-value rasters onto the other side's probability grid
when the prior and likelihood sources report different percentiles (e.g. one source's 5th/95th vs.
another's 2.5th/97.5th) — every downstream fusion route fits both sides using one shared `probs`
vector, so mismatched probabilities passed through unchanged would silently mis-fit whichever
side's actual percentiles don't match the assumed labels.

**Parameters**:
```r
align_percentile_probs(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs)
# prior_value_rasters, lik_value_rasters - lists of percentile-value SpatRasters
# prior_probs, lik_probs                 - matching probability vectors
```
**Returns**: `list(prior_value_rasters =, lik_value_rasters =, probs =)` — both sides re-expressed
onto a single shared `probs` vector.

**Behavior**: If `prior_probs`/`lik_probs` are already identical, returns the inputs unchanged. If
not, the *narrower*-ranging side's probabilities become the shared target (`diff(range(...))`
comparison), and the wider-ranging side is reinterpolated onto that target via
`quantile_linear_cdf_raster()` — the narrower side is always safely interpolatable from the wider
one, never the reverse, since `quantile_linear_cdf_raster()` requires the target probability to
fall within the source's own range. Reinterpolating outside that range `stop()`s.

#### 10. `fuse_adaptive()`
**Purpose**: Fuse a prior and likelihood belief distribution (normal/beta/gamma families only),
choosing the fusion *route* by AOI cell count rather than requiring the caller to pick — the
size-adaptive engine underneath `fuse_property_adaptive()`.

**Parameters**:
```r
fuse_adaptive(prior_value_rasters, lik_value_rasters, percentile_probs,
              family = c("normal", "beta", "gamma"), bounds = NULL,
              threshold_cells = 80000, n_samples = 500,
              grid_resolution = NULL, verbose = TRUE)
# prior_value_rasters, lik_value_rasters - named lists of percentile-value SpatRasters
#                                           (e.g. list(P0=, P5=, P50=, P95=, P100=)), same
#                                           layer count/order on both sides
# percentile_probs  - numeric probabilities matching the value-raster list order
# family            - "normal" | "beta" | "gamma" - fit family for the closed-form route, and
#                      the family the general route's output is moment-matched back into
# bounds            - required for family = "beta" (c(lower, upper))
# threshold_cells   - AOI cell count at or below which the general bayesian_update() route
#                      is used; above it, the closed-form route
# n_samples         - samples drawn per side for the general route's KDE fusion
# grid_resolution   - passed to bayesian_update() for the general route
# verbose           - if TRUE (default), print the chosen route and why
```
**Returns**: `list(posterior =, route =, route_detail =, n_fallback_cells =, diagnostics =)`.
`posterior` is always a family-native parameter list (`mu`/`sigma` for "normal", `alpha`/`beta` for
"beta", `shape`/`rate` for "gamma") regardless of which route ran. `route` is one of
`"bayesian_update_general"`, `"closed_form_normal"`, `"closed_form_beta"`, `"closed_form_gamma"`.
`route_detail` is a `SpatRaster` of 1 (native same-family fusion)/0 (per-cell fallback) for
`family %in% c("beta","gamma")` on the closed-form route, `NULL` otherwise. `diagnostics` is
`list(ncell=, threshold_cells=, family=, elapsed_sec=)`.

**Behavior**: `use_general <- ncell <= threshold_cells`. Below/at threshold, delegates to the
internal `fuse_general_kde()`: per cell (via `terra::app()` over a combined prior+likelihood
percentile stack), draws `n_samples` samples from each side via
`simulate_from_percentiles(method = "linear_cdf")`, fuses via `bayesian_update()`, and
moment-matches the posterior samples back into the requested family — wrapped in a per-cell
`tryCatch()` so a cell with too many `NA` percentile values (e.g. at raster edges after
`terra::resample()`/`align_percentile_probs()` alignment) degrades to `NA` output instead of
aborting the whole raster operation. Above threshold, delegates to the internal
`fuse_closed_form()`: fits both sides in the requested family (`fit_normal_raster()` for normal,
`fit_beta_mle_newton_raster()` for beta, `fit_gamma_mom_raster()` for gamma) and fuses via
`core-fusion.R`'s `bayes_update_normal_normal()`/`fuse_beta()`/`fuse_gamma()`; for beta/gamma,
cells where the native same-family fusion is infeasible fall back per-cell to Normal-moment fusion
(`bayes_update_normal_normal()` on each side's moment-derived mean/sd), re-expressed back into the
requested family via `moments_to_beta()`/`moments_to_gamma()`, with `n_fallback_cells` counting how
many cells needed it. `verbose` prints the chosen route (and, if any cells fell back, how many).

#### 11. `resolve_property_dist()`
**Purpose**: Resolve a property's distribution family, refining a `dist = "auto"` config value from
the AOI's own percentile skew rather than trusting a fixed per-property label — never overrides an
explicit (non-`"auto"`) `dist`.

**Parameters**:
```r
resolve_property_dist(property_config, prior_value_rasters, prior_probs,
                       aggregate_fun = stats::median)
# property_config       - list with dist and optionally bounds/auto_skew_threshold
# prior_value_rasters   - prior percentile-value rasters
# prior_probs           - prior percentile probabilities
# aggregate_fun         - function to aggregate the AOI-wide skew proxy (default stats::median)
```
**Returns**: `list(dist =, dist_source = one of "config"/"auto", skew_proxy = NA_real_ unless
dist_source == "auto")`.

**Behavior**: If `property_config$dist` isn't literally `"auto"`, returns it unchanged with
`dist_source = "config"`. Otherwise computes a skewness proxy from the prior's low/median/high
percentile rasters, `((high-median)-(median-low))/(high-low)`, aggregated once across the whole AOI
via `aggregate_fun` (default `median`) — deliberately *not* per-cell, which would create map
discontinuity artifacts at family-boundary cells. Resolves to `"beta"` if `property_config$bounds`
is set, else `"normal"` if `abs(skew_proxy) < threshold` (default `0.15`, overridable via
`property_config$auto_skew_threshold`), else `"lognormal"`. Deliberately narrow: never resolves to
`"metalog"` or `"gamma"`, since distinguishing those from a 3-point proxy isn't considered reliable.

#### 12. `fuse_property_adaptive()` — main pipeline entry point
**Purpose**: Fuse one property's prior and likelihood percentile rasters, dispatching automatically
to the right distribution family and fusion route from a `property_config` list. This is the
raster-native counterpart of `core-fusion.R`'s `fuse_property()`, the toolkit's top-level
entry point for non-compositional properties (see `fuse_texture_group()` for the compositional
case), and the function every adapter (SSURGO prior + SOLUS likelihood, or any future pair) is
ultimately built to feed.

**Parameters**:
```r
fuse_property_adaptive(prior_value_rasters, prior_probs,
                        lik_value_rasters, lik_probs,
                        property_config, threshold_cells = 80000, ...)
# prior_value_rasters, lik_value_rasters - named lists of percentile-value SpatRasters
# prior_probs, lik_probs   - matching probability vectors (reconciled via
#                            align_percentile_probs() if they differ)
# property_config          - list with dist (one of "normal"/"beta"/"gamma"/"lognormal"/
#                            "metalog"/"auto"), and optionally bounds/boundedness/
#                            auto_skew_threshold
# threshold_cells          - AOI cell count at or below which fuse_adaptive() uses its
#                            general route
# ...                      - additional arguments passed to fuse_adaptive()/
#                            fuse_lognormal_adaptive() (e.g. n_samples, grid_resolution, verbose)
```
**Returns**: `list(posterior =, route =, route_detail =, n_fallback_cells =, dist =, dist_source =,
skew_proxy =)`.

**Behavior**: First calls `align_percentile_probs()` so both sides share one probability grid, then
`resolve_property_dist()` to pick the working `dist`. Routing then branches on `dist`: `"normal"`/
`"beta"`/`"gamma"` go straight to `fuse_adaptive(family = dist)` (Function 10 above, itself
size-adaptive); `"lognormal"` goes to the internal `fuse_lognormal_adaptive()`, which — below
`threshold_cells` — forces `fuse_adaptive()` into its general route directly on raw (non-log)
values (since `bayesian_update()`'s general route makes no distributional assumption and needs no
log-space detour), or — above threshold — fits both sides Normal, transforms to lognormal parameter
space via `normal_to_lognormal_params()`, fuses via `bayes_update_normal_normal()`, and transforms
back via `lognormal_to_normal_params()`; `"metalog"` goes to the internal `fuse_metalog_adapter()`,
which fits both sides via `fit_metalog_linear_raster()` on interior percentiles (always, regardless
of AOI size), checks feasibility via `check_metalog_feasibility_raster()`, computes each side's
mean/sd via `metalog_moments_raster()` (quadrature over the raw metalog quantile function, with a
Normal-moment fallback for infeasible cells), and fuses those moments via
`bayes_update_normal_normal()`. Any other `dist` value `stop()`s. The returned list always carries
`dist`/`dist_source`/`skew_proxy` from `resolve_property_dist()` alongside whichever route's
posterior/route/route_detail/n_fallback_cells.

#### 13. `group_members()`
**Purpose**: Look up a compositional group's member property ids, in configured order, reusing
soilSIM's own `config$monte_carlo$composition_groups` convention (see `R/core-distributions.R`'s
`resolve_composition_groups()`) rather than a separate config schema.

**Parameters**:
```r
group_members(group, composition_groups)
# group               - composition group name (e.g. "texture")
# composition_groups  - a config$monte_carlo$composition_groups-shaped list, e.g.
#                        list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))
```
**Returns**: A character vector of member property ids, in configured order.

**Behavior**: Looks up `composition_groups[[group]]`, `stop()`s if not found or if it has fewer
than 2 members; otherwise returns `group_def$members` verbatim (order trusted from config, not
hardcoded to literal clay/sand/silt names).

#### 14. `fuse_texture_group()`
**Purpose**: Fuse a compositional group's members (e.g. clay/sand/silt) JOINTLY via ILR
(isometric log-ratio) fusion, rather than independently fusing each member with `fuse_beta()` —
independent fusion measurably breaks sum-to-100 (up to about 10 percentage points on realistic
synthetic data). The raster counterpart of the scalar `fuse_texture_group_from_triplets()`
(`core-fusion.R`).

**Parameters**:
```r
fuse_texture_group(fetched)
# fetched - a list (one entry per group member, in group_members() order) of
#           list(id=, prior=<named list of ALIGNED percentile-value rasters>, prior_probs=,
#                lik=<value rasters>, lik_probs=)
```
**Returns**: A named list keyed by each member's `id`: `list(posterior = list(value = <inverse-ILR
point-estimate raster for that fraction>, ilr_mu = <2-layer raster, shared across the group>,
ilr_Sigma = <3-layer raster (S11,S12,S22), shared>), dist = "texture_ilr", route =
"closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)`.

**Behavior**: Requires exactly 3 members (`stopifnot(length(fetched) == 3)` — clay/sand/silt).
Builds one combined `terra::rast()` stack of every member's prior/likelihood
low/median/high-percentile layers, then per cell (via `terra::app()`) calls `R/core-distributions.R`'s
`estimate_ilr_moments_mc()` on each side's clay/sand/silt triplet to get `mu`/`Sigma` in
2-dimensional ILR space, fuses the two sides via `core-fusion.R`'s `fuse_bivariate_normal()`,
and returns the fused `mu1`/`mu2`/`S11`/`S12`/`S22` per cell. Afterward, `ilr_inverse()` converts the
fused `mu1`/`mu2` raster values back to a 3-column point estimate (one column per fraction) for a
point-estimate output layer per member; the full `ilr_mu`/`ilr_Sigma` rasters are also returned
(shared across all three members) so callers can draw posterior samples for a specific fraction/
cell via `sample_ilr_posterior()`.

#### 15. `run_stage1_fusion()` — top-level AOI orchestrator (non-compositional)
**Purpose**: Run the complete prior-fetch → likelihood-fetch → align → fuse pipeline for one
property/depth window over an AOI, with disk caching at every fetch step.

**Parameters**:
```r
run_stage1_fusion(aoi_vect, property_config, top_depth, bottom_depth,
                   composition_groups = NULL, property_configs = NULL)
# aoi_vect            - a terra::SpatVector AOI (projected, e.g. EPSG:5070)
# property_config     - list with id (used as fetch_ssurgo_percentiles()'s property id and as
#                        the cache key), solus_variable (a soilDB::fetchSOLUS()-recognized
#                        variable name), and dist/bounds/boundedness/auto_skew_threshold/
#                        composition_group as needed by fuse_property_adaptive()
# top_depth, bottom_depth - numeric depth bounds in cm
# composition_groups, property_configs - only needed when property_config$composition_group
#                        is set; passed through to run_stage1_fusion_group()
```
**Returns**: `list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=,
route_detail=, n_fallback_cells=)`, or `NULL` if the SSURGO or SOLUS side failed. `posterior`'s
shape depends on `dist` — see `fuse_property_adaptive()`.

**Behavior**: If `property_config$composition_group` is set, delegates entirely to
`run_stage1_fusion_group()` and slices out this property's own result. Otherwise: builds a cache key
via `build_cache_key(..., "ssurgo")`, checks `cache_get()`/`unwrap_percentile_list()`, and on a miss
calls `fetch_ssurgo_percentiles()` and caches the result (`wrap_percentile_list()` first, per
`cache.R`'s known limitation); repeats the same cache-or-fetch pattern for the SOLUS side
via `fetch_solus_percentiles()` under kind `"solus"`. Because the SSURGO prior (~30 m, from
`mukey.wcs()`) and the SOLUS likelihood (100 m) sit on different grids, the prior is
`terra::resample(..., method = "bilinear")`'d onto the SOLUS grid before fusion, since every fusion
route requires cell-aligned rasters. Finally calls `fuse_property_adaptive()` with the aligned
prior and the SOLUS likelihood, and returns both raw sides plus the fused posterior.

The SSURGO simulation is restricted to *just this call's own property* (plus the full
sand/silt/clay draw for a texture member) via `simulate_ssurgo_mapunit_draws()`'s
`requested_properties` argument — the per-cokey depth-trend GP fit is the pipeline's dominant cost
and scales with the property count. Two consequences: the prior is *statistically equivalent* to,
not bit-identical to, a full-property simulation's matching column (the pipeline is unseeded — a
fresh draw either way), and a map unit whose components carry no data for *this one* property no
longer benefits from those components' other properties being present, so a single-property prior
can be `NA` for a few cells a full-property run would have covered. The restriction is honoured
only under the default `vertical_correlation_method = "joint_copula"`; `run_stage1_fusion()` never
forwards a `config`, so it always gets that method. An optional `seed` argument makes the whole
call reproducible given identical upstream SSURGO/SOLUS data; `NULL` (default) is stochastic.

#### 16. `run_stage1_fusion_group()` — top-level AOI orchestrator (compositional)
**Purpose**: Run the same fetch → align → fuse pipeline as `run_stage1_fusion()`, but for a whole
compositional group (currently only `"texture"`: clay/sand/silt) jointly, since independent fusion
per member measurably breaks sum-to-100.

**Parameters**:
```r
run_stage1_fusion_group(aoi_vect, group, composition_groups, property_configs,
                         top_depth, bottom_depth)
# aoi_vect                 - a terra::SpatVector AOI
# group                    - a composition group name (e.g. "texture")
# composition_groups       - a config$monte_carlo$composition_groups-shaped list (see
#                             group_members()) giving this group's member ids in order
# property_configs         - a named list, keyed by each member id composition_groups lists,
#                             of per-member config lists (each needs at least id/solus_variable)
# top_depth, bottom_depth  - numeric depth bounds in cm
```
**Returns**: A named list keyed by member id, each element `list(posterior = list(value=, ilr_mu=,
ilr_Sigma=), dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL,
n_fallback_cells = 0)` — NOT the uniform `(mu,sigma)`/`(alpha,beta)` contract non-compositional
properties get — or `NULL` if any member's SSURGO/SOLUS fetch failed.

**Behavior**: Resolves member ids via `group_members()`. Checks a group-level cache entry (kind
`"texture_group"`) first; on a miss, fetches each member's SSURGO prior and SOLUS likelihood
(cache-or-fetch, exactly as in `run_stage1_fusion()`, one `"ssurgo"`/`"solus"` cache entry per
member so a later single-property request still hits cache), resamples each member's prior onto its
own SOLUS grid, and — only if every member's fetch succeeded — calls `fuse_texture_group()` once for
the whole group. The joint result is cached under the group key, and each member's own posterior is
*also* seeded into a per-property `"posterior"` cache kind, so three sequential
`run_stage1_fusion()` calls for clay, then sand, then silt trigger the joint fetch+fusion exactly
once, not three times. The shared simulation is itself restricted to the group's members (which for
`"texture"` normalizes to "just the sand/silt/clay draw"), skipping the other properties' GP fits.

#### 16b. `run_stage1_fusion_multi()` — many properties, one simulation

**Purpose**: Fuse many properties (and many depth windows) over an AOI from a **single** SSURGO
Monte Carlo simulation, instead of one `run_stage1_fusion()` call per property/window — each of
which re-runs the full (dominant-cost) simulation only to keep one of its ~9 jointly-simulated
columns.

**Parameters**:
```r
run_stage1_fusion_multi(aoi_vect, property_configs, depth_windows,
                         composition_groups = NULL, n_mc = 1000,
                         parallel = FALSE, n_cores = NULL,
                         verbose = FALSE, simplify = FALSE)
# property_configs  - a NAMED list, keyed by each config's own id, of the same property_config
#                     lists run_stage1_fusion() takes; a config carrying composition_group is fused
#                     jointly per run_stage1_fusion_group(), and every member of a referenced group
#                     must have its own entry
# depth_windows     - a non-empty list of c(top, bottom) numeric pairs
# simplify          - TRUE reduces each leaf to list(percentiles = <named SpatRasters>)
```
**Returns**: a nested named list `result[[config_id]][["<top>-<bottom>"]]`, each leaf a
`run_stage1_fusion()`-shaped list (or the `texture_ilr` shape for group members), `NULL` on that
leaf's own failure; the whole call returns `NULL` only if the shared simulation itself fails.

**Behavior**: Validates the configs, partitions them into standalone vs compositional-group, then
runs `simulate_ssurgo_mapunit_draws(span, depth_windows =, requested_properties = <union of every
config's property>)` **once** (skipping it entirely when every needed `"ssurgo"` cache is already
warm and no config uses raw-draws fusion). SOLUS is fetched via one batched
`fetch_solus_percentiles_multi()` call **per window** (covering every variable any standalone
config or group member needs at once, S1), falling back to the per-variable scalar
`fetch_solus_percentiles()` only if that batched request itself errors; results are memoized. Per
window then per property it calls the same fusion tail `run_stage1_fusion()` uses
(`stage1_fuse_from_prior_solus()`), and for texture groups the same one `run_stage1_fusion_group()`
uses (`stage1_fuse_texture_group_from_fetched()`), **seeding the same per-property
`"ssurgo"`/`"solus"`/`"posterior"` disk caches** — so a later single-property `run_stage1_fusion()`
for the same AOI/property/window is a cache hit. Every leaf therefore shares one draw set, so
(unlike independent `run_stage1_fusion()` calls) the properties' `NA` masks are mutually
consistent. The per-window draws are released as each window finishes.

#### 17. `build_cache_key()` / `CACHE_TTL_SECONDS`
**Purpose**: Build a deterministic, filename-safe cache key for one AOI/property/depth-window/kind
combination, and the default cache max-age (30 days) used throughout the pipeline.

**Parameters**:
```r
build_cache_key(aoi_vect, id, top_depth, bottom_depth, kind)
# aoi_vect               - a terra::SpatVector (single-feature AOI)
# id                     - a property id or composition-group name (e.g. "ph", "texture")
# top_depth, bottom_depth - numeric depth bounds in cm
# kind                   - a short string distinguishing what's cached for this key
#                          (e.g. "ssurgo", "solus", "texture_group", "posterior")

CACHE_TTL_SECONDS  # = 30 * 24 * 3600 (30 days)
```
**Returns**: A character string, safe for use as a filename.

**Behavior**: Extracts the AOI's WKT geometry (`terra::geom(aoi_vect, wkt = TRUE)[1]`), hashes it
via `digest::digest()`, and concatenates a `"raster_"` prefix, the first 12 hash characters, `id`,
`top_depth`, `bottom_depth`, and `kind`, then sanitizes any character outside `[A-Za-z0-9_-]` to an
underscore.

#### 18. `cache_get()` / `cache_set()`
**Purpose**: Retrieve/store a cached value by key, under `tools::R_user_dir("soilSIM", "cache")`
(the standard R >= 4.0 per-package cache location) — the general disk-cache mechanism every
adapter and orchestrator above uses.

**Parameters**:
```r
cache_get(key, ttl_seconds = CACHE_TTL_SECONDS)
# key         - a cache key from build_cache_key()
# ttl_seconds - maximum cache age in seconds before a hit is treated as stale

cache_set(key, kind, value)
# key   - a cache key from build_cache_key()
# kind  - unused; documents intent at call sites. The value is stored keyed only by key, since
#         build_cache_key() already encodes kind
# value - the R object to cache (anything saveRDS() can serialize)
```
**Returns**: `cache_get()` returns the cached value, or `NULL` on a cache miss/stale entry/read
failure. `cache_set()` returns `TRUE` on success, `FALSE` on a write failure (never errors).

**Behavior**: `cache_get()` checks file existence, computes age from `file.mtime()`, returns `NULL`
if stale or missing, otherwise `readRDS()`s (wrapped in `tryCatch`). `cache_set()` creates the cache
directory if needed and `saveRDS()`s the value (wrapped in `tryCatch`) — see Known Limitations for
the `SpatRaster` wrap/unwrap caveat that both this function and its callers must handle manually.

### SSURGO adapter (adapter-ssurgo-simulate.R)

#### 19. `fetch_ssurgo_mukey_raster()`
**Purpose**: Fetch a categorical (factor) raster of SSURGO map unit keys covering an AOI — the
first step of the SSURGO prior pipeline.

**Parameters**:
```r
fetch_ssurgo_mukey_raster(aoi_vect)
# aoi_vect - a terra::SpatVector (projected, e.g. EPSG:5070) - the AOI footprint itself, not
#            widened to its bounding box
```
**Returns**: A single-layer, categorical `terra::SpatRaster` named `"mukey"`, or `NULL` if either
the grid or the polygon fetch returns nothing.

**Behavior**: Calls `soilDB::mukey.wcs(aoi = aoi_vect, db = "gssurgo")` for the grid, then
`soilDB::SDA_spatialQuery(aoi_vect, what = "mupolygon", db = "SSURGO", geomIntersection = TRUE)`
for the map-unit polygons (wrapped in `tryCatch` — a real, encountered failure mode is
`SDA_spatialQuery()` erroring outright when its `sf` dependency isn't loadable in-session, not just
a hypothetical). The polygons are reprojected to the grid's CRS, rasterized onto it by `mukey`
value via `terra::rasterize()`, and cast to a factor raster via `terra::as.factor()`.

#### 20. `infill_soil_data()`
**Purpose**: Infill missing values across the standard SSURGO property set for one horizon data
frame.

**Parameters**:
```r
infill_soil_data(df, water_retention_method = c("saxton_rawls", "generic"))
# df                    - a horizon data frame, as returned by download_ssurgo_tabular()
# water_retention_method - how to fill wthirdbar/wfifteenbar gaps (see Behavior)
```
**Returns**: `df` with missing values infilled where possible.

**Behavior**: A thin wrapper around `process_soil_properties_comprehensive()` (the single
infilling orchestrator - see `docs/01_data_acquisition_processing.md` entry 15), restricting it
to whichever of the standard SSURGO property set (`sandtotal`/`claytotal`/`silttotal`/`dbovendry`/
`om`/`rfv`/`wthirdbar`/`wfifteenbar`/`ph1to1h2o`/`cec7` + the 5 P2 chemistry properties) is
actually present in `df`, and passing `water_retention_method` through. That orchestrator runs
three phases - foundation (texture/BD/OM/RFV), water retention (Saxton-Rawls pedotransfer where
texture + BD allow, `overwrite = FALSE`; generic hierarchy for the rest or under
`water_retention_method = "generic"`), then the remaining chemical properties - each property
going through `infill_soil_property()`'s six-strategy hierarchy.

#### 21. `maybe_adjust_soil_data_depth_trend()`
**Purpose**: Apply depth-trend GP (Gaussian process) adjustment per cokey, guarded by whether the
`GPfit` package is installed.

**Parameters**:
```r
maybe_adjust_soil_data_depth_trend(sim_long, properties, min_depths = 2)
# sim_long   - long-format simulated data with cokey, hzdept_r, and property columns
# properties - character vector of property column names to adjust
# min_depths - minimum distinct depths required to attempt GP fitting (default 2)
```
**Returns**: `sim_long`, depth-trend-adjusted where possible.

**Behavior**: Warns and returns `sim_long` unchanged if `GPfit` isn't installed. Otherwise, groups
by `cokey` and, for each cokey with at least `min_depths` distinct non-`NA` depths, applies
`R/core-gp.R`'s `apply_local_gp_adjustments()` (which fits its own local GP per
property from that cokey's own within-simulation depth trend — no pre-supplied GP models needed);
cokeys below the threshold pass through unadjusted.

#### 22. `aggregate_depth_window_by_replicate()`
**Purpose**: Collapse long-format simulated property data down to one thickness-weighted mean row
per replicate, over a requested depth window.

**Parameters**:
```r
aggregate_depth_window_by_replicate(sim_long, top_depth, bottom_depth, property_cols,
                                     replicate_cols = c("mukey", "cokey", "simulation_number"))
# sim_long                - long-format simulated data with hzdept_r/hzdepb_r and replicate_cols
# top_depth, bottom_depth - numeric depth window bounds in cm
# property_cols           - character vector of property columns to aggregate
# replicate_cols          - columns identifying one replicate
```
**Returns**: One row per replicate, with `top`/`bottom` set to the requested window and each
property column replaced by its overlap-weighted mean (`NA` if the replicate has no overlap with
the window).

**Behavior**: Computes each horizon's overlap with `[top_depth, bottom_depth]` as
`pmax(0, pmin(hzdepb_r, bottom_depth) - pmax(hzdept_r, top_depth))`, drops zero-overlap horizons,
and for each `replicate_cols` group takes `stats::weighted.mean()` of every property column weighted
by that overlap.

#### 23. `property_to_sim_column()`
**Purpose**: Map a caller-facing property id to `simulate_cokey_generalized()`'s output column name.

**Parameters**:
```r
property_to_sim_column(property_id)
# property_id - one of "ph", "bulk_density", "soc", "cec", "clay", "sand", "silt",
#               "rock_fragments", "wthirdbar", "wfifteenbar", "caco3", "ec", "ecec",
#               "gypsum", "sar" (+ synonyms)
```
**Returns**: The corresponding column name (e.g. `"clay"` -> `"clay_total"`); `stop()`s on an
unrecognized id.

**Behavior**: A fixed lookup table (`ph`, `bulk_density -> db`, `soc`, `cec`, `clay -> clay_total`,
`sand -> sand_total`, `silt -> silt_total`, `rock_fragments -> rfv`,
`wthirdbar/water_retention_third_bar -> wr_3b`, `wfifteenbar/water_retention_15_bar -> wr_15b`,
`caco3`/`ec`/`ecec`/`gypsum`/`sar` -> themselves). The water-retention rows were added 2026-09-02 -
`simulate_cokey_generalized()` already emitted those columns (`SSURGO_SIM_PROPERTY_COLUMNS` lists
them) but the id map was incomplete, blocking water-retention `run_stage1_fusion()` (needed by
`remarginalized_awc()`). The 5 chemistry-property rows self-map, since the SOLUS variable name, the SSURGO
chorizon column stem, and the id are all identical spellings for these 5.

#### 24. `simulate_ssurgo_mapunit_draws()`
**Purpose**: Orchestrate the full SSURGO Monte Carlo pipeline for one AOI/depth window — fetch
tabular SSURGO data, infill, simulate, adjust, and aggregate — producing the replicate-level draws
that `fetch_ssurgo_percentiles()` summarizes into percentile rasters.

**Parameters**:
```r
simulate_ssurgo_mapunit_draws(aoi_vect, top_depth, bottom_depth, n_mc = 1000,
                              ..., depth_windows = NULL, requested_properties = NULL)
# aoi_vect                 - a terra::SpatVector AOI
# top_depth, bottom_depth  - numeric depth window bounds in cm
# n_mc                     - number of triangular draws sim_component_comp() uses per
#                            component (default 1000)
# depth_windows            - optional list of c(top, bottom) pairs; simulate once, aggregate per
#                            window, return a named list of data frames (top_depth/bottom_depth
#                            then just give the overall span)
# requested_properties     - NULL simulates every property (unchanged); otherwise a vector (any
#                            vocabulary normalize_requested_properties() accepts) restricting the
#                            simulation - the per-cokey depth-trend GP fit then runs once per kept
#                            property. Ignored (with a warning) under
#                            vertical_correlation_method = "gp_quantile_retrofit", which is not
#                            subset-invariant
# seed                     - optional integer for opt-in determinism (NULL = current stochastic
#                            behavior); set.seed(seed) once up front + future.seed on the parallel
#                            depth-trend path. Conditional on unchanged upstream SSURGO data, and
#                            does not make a requested_properties subset bit-identical to a full run
```
**Returns**: A data frame, one row per `mukey`/`cokey`/`simulation_number` replicate, with simulated
property columns aggregated over the depth window — or a named list of such data frames when
`depth_windows` is supplied — or `NULL` if the tabular fetch fails.

**Behavior**: Reprojects `aoi_vect` to EPSG:4326 first (`download_ssurgo_tabular()`'s
`aoi_wkt` argument is hardcoded to assume lon/lat, regardless of `aoi_vect`'s actual CRS — extracting
WKT directly from an already-projected AOI, e.g. EPSG:5070, would silently mislabel meter
coordinates as degrees). Checks a `"ssurgo_tabular"`-kind disk cache first, keyed by **AOI alone**
via `ssurgo_tabular_cache_key()` (depth-independent as of task L1 — `download_ssurgo_tabular()`
takes no depth-window argument and fetches every horizon for every AOI mukey regardless, so one
AOI's tabular download is shared across every depth window/call requested for it); on a miss, calls
`R/adapter-ssurgo-acquire.R`'s `download_ssurgo_tabular()` and unwraps its `$ssurgo_data`, caching that
data frame directly (no `wrap()` needed — it's tabular, not a raster). Then: `infill_soil_data()`;
removes organic horizons (`remove_organic_layer()`) if any `hzname` contains `"O"`; derives `genhz`
via `classify_genhz()`; simulates component composition via `sim_component_comp()` and left-joins
`sim_comppct` back onto the horizon data by `cokey`; loads the KSSL reference correlation matrices
via `.kssl_property_matrices()`/`.kssl_texture_matrices()` (already built into `R/sysdata.rda`, not
re-read from raw files at runtime); per unique `cokey`, calls `simulate_cokey_generalized()`
(wrapped in `tryCatch`, skipping and messaging on a per-cokey failure — passing
`requested_properties` through so each `simulate_cokey_generalized()` call only draws the kept
subset) and row-binds all results; applies `maybe_adjust_soil_data_depth_trend()` to whichever of
`SSURGO_SIM_PROPERTY_COLUMNS` actually came back (so a restricted simulation fits fewer GP models);
and finally aggregates to the requested depth window(s) via `aggregate_depth_window_by_replicate()`.
`property_matrices` (the genhz-keyed correlation matrices passed to `simulate_cokey_generalized()`)
is now built via `build_kssl_fallback_matrix()` per genhz key, sized to the full 14-name
`param_order` vocabulary, rather than indexing `.kssl_property_matrices()` directly - see #25's
P2 note for why.

#### 25. `SSURGO_SIM_PROPERTY_COLUMNS`
**Purpose**: The fixed set of simulated property column names this pipeline aggregates/adjusts.

**Value**: `c("db", "wr_3b", "wr_15b", "rfv", "ph", "cec", "soc", "sand_total", "silt_total",
"clay_total", "caco3", "ec", "ecec", "gypsum", "sar")`.

**Note on `soc`**: SSURGO's `chorizon` table reports organic matter (`om_l/_r/_h`), not organic
carbon. `simulate_cokey_generalized()` converts it to an estimated SOC value
(`om * OM_TO_SOC_FACTOR`, the Van Bemmelen factor, `1/1.724`) before it enters the correlated
draw, so this `soc` column - and everything fused against SOLUS100's real `soc` variable via
`solus_variable = "soc"` - is SOC-scale, not raw OM. The KSSL
reference correlation matrix's `soc` correlations were still fit on OM data and needed no change
(correlation is invariant under a linear rescale of one variable) - treat those *correlations* as
an OM-based approximation, but the simulated *values* as genuine SOC estimates.

**Note on the 5 chemistry properties (`caco3`/`ec`/`ecec`/`gypsum`/`sar`)**: SSURGO chorizon
column names confirmed against a real gSSURGO query. They have
**no** entry in the static KSSL reference matrices (`.kssl_property_matrices()`) - `simulate_ssurgo_mapunit_draws()` builds its per-genhz correlation
matrices via `build_kssl_fallback_matrix()` (sized to the full `param_order` vocabulary) so
these 5 simulate as *uncorrelated* with everything else (and each other) rather
than erroring: identity fallback now, real correlations after a future KSSL-data re-fit that
needs no code changes here (see `build_kssl_fallback_matrix()`'s own docs). SSURGO coverage for these 5 is real but noticeably sparser in spot-checks (e.g.
`ecec` was populated in only ~25% of sampled horizons for one AOI) - expect more `NA` in their fused output - a SSURGO data-completeness gap, not a pipeline defect.

#### 26. `rasterize_mukey_percentiles()`
**Purpose**: Merge a per-mukey percentile-value table onto a mukey raster, producing one
`SpatRaster` per percentile column.

**Parameters**:
```r
rasterize_mukey_percentiles(mukey_raster, percentile_by_mukey)
# mukey_raster        - a categorical (factor) mukey terra::SpatRaster from
#                       fetch_ssurgo_mukey_raster()
# percentile_by_mukey - a data frame with a mukey column and one column per percentile
#                       (e.g. P05, P25, P50, P75, P95)
```
**Returns**: A named list of single-layer `terra::SpatRaster`s, one per percentile column.

**Behavior**: Reads the raster's factor category table (`terra::cats()`), whose first column is the
internal numeric factor ID (must stay first, since `levels<-` requires it) and whose second column
holds the mukey codes as character. Left-joins (`dplyr::left_join()`, deliberately not base
`merge()`, which reorders the join column to the front and would break the "ID must stay first"
requirement) the percentile columns onto that category table by the mukey column, reassigns
`levels(mukey_raster) <- merged`, then calls `terra::catalyze()` to expand the categorical raster
into one continuous layer per joined percentile column, returned as a named list.

#### 27. `fetch_ssurgo_percentiles()` — SSURGO prior entry point
**Purpose**: The top-level SSURGO "prior" entry point for `fuse_property_adaptive()`: rasterizes
map units, Monte Carlo-simulates the requested property, and returns per-cell percentile-value
rasters in the `list(values=, probs=)` shape the fusion core expects.

**Parameters**:
```r
fetch_ssurgo_percentiles(aoi_vect, property_id, top_depth, bottom_depth,
                          probs = c(0.05, 0.25, 0.5, 0.75, 0.95), n_mc = 1000,
                          ..., requested_properties = NULL)
# aoi_vect                 - a terra::SpatVector AOI
# property_id              - one of property_to_sim_column()'s recognized ids
# top_depth, bottom_depth  - numeric depth window bounds in cm
# probs                    - percentile probabilities to compute
# n_mc                     - number of triangular draws per component (default 1000)
# requested_properties     - NULL means "just property_id" (this function only reads one column
#                            back out), so a bare call simulates only what it needs; pass a vector
#                            to keep extra properties in the shared draws
```
**Returns**: `list(values = <named list of percentile-value SpatRasters>, probs = probs)`, or
`NULL` if the mukey raster or the Monte Carlo draws are unavailable for this AOI.

**Behavior**: Calls `fetch_ssurgo_mukey_raster()` (returning `NULL` early if it fails), then
`simulate_ssurgo_mapunit_draws()` for the replicate-level Monte Carlo draws — restricted to
`requested_properties` (defaulting to just `property_id`). Maps `property_id` to its simulated
column via `property_to_sim_column()`, computes `stats::quantile(..., probs = probs)` per `mukey`
group, reshapes into a `mukey` + one-column-per-percentile data frame, and rasterizes via
`rasterize_mukey_percentiles()`.

### SOLUS100 adapter (adapter-solus.R)

#### 28. `closest_solus_depth_slice()`
**Purpose**: Snap a `[top_depth, bottom_depth]` depth window to the nearest native SOLUS depth
slice, since `soilDB::fetchSOLUS()`'s `depth_slices` are single depth points, not ranges - the
simplest representative-point choice for an otherwise point-based data source.

**Parameters**:
```r
closest_solus_depth_slice(top_depth, bottom_depth, available_slices = c(0, 5, 15, 30, 60, 100, 150))
# top_depth, bottom_depth - numeric depth window bounds in cm
# available_slices        - native SOLUS depth points to snap to
```
**Returns**: A single numeric value from `available_slices`.

**Behavior**: Computes the window midpoint `(top_depth + bottom_depth) / 2` and returns whichever
`available_slices` value minimizes `abs(available_slices - midpoint)`.

**Relationship to the fusion pipeline**: `fetch_solus_low_pred_high()` uses
`solus_depth_window_weights()` for a trapezoidal depth-average across native slices, so the
SOLUS likelihood and SSURGO prior both represent a window mean. This function is public API for
callers that want a single nearest slice.

#### 28b. `solus_depth_window_weights()`
Trapezoidal depth-average weights over the native SOLUS slices spanning a window: `f(d)`
piecewise-linear between slices (`rule = 2` clamp outside `[0, 150]`), trapezoidal integration over
window ends + interior slices, normalized by window thickness. Returns a named numeric vector
summing to 1. Exact for a linear profile (== midpoint value); reduces SOLUS point predictions to the
same estimand the SSURGO prior side computes via `aggregate_depth_window_by_replicate()`.

#### 29. `fetch_solus_low_pred_high()`
**Purpose**: Fetch SOLUS100's low/prediction/high rasters for one variable and depth window — the
raw three-way fetch underlying `fetch_solus_percentiles()`.

**Parameters**:
```r
fetch_solus_low_pred_high(aoi_vect, solus_variable, top_depth, bottom_depth)
# aoi_vect                 - a terra::SpatVector AOI
# solus_variable           - a soilDB::fetchSOLUS()-recognized variable name (e.g. "claytotal",
#                            "dbovendry", "ph1to1h2o")
# top_depth, bottom_depth  - numeric depth window bounds in cm
```
**Returns**: `list(pred=, low=, high=)`, each a single-layer `terra::SpatRaster` or `NULL` if that
output type wasn't returned.

**Behavior**: Computes `solus_depth_window_weights(top_depth, bottom_depth)`, then calls
`soilDB::fetchSOLUS()` once per `output_type` (`"prediction"`, `"95% low prediction interval"`,
`"95% high prediction interval"`, each wrapped in `tryCatch` with a warning on failure) requesting
every contributing native slice, then reduces the returned layers
(`paste0(solus_variable, "_", <slice>, "_cm_", suffix)`, suffixes `p`/`l`/`h`) to one raster per
output type via `terra::app(stack * weights, "sum")` — the trapezoidal depth-average. `NULL` for an
output type if any required slice layer is missing. The output-layer-naming assumption was verified
against a live, network-connected `fetchSOLUS()` call before porting.

#### 30. `fetch_solus_percentiles()` — SOLUS likelihood entry point
**Purpose**: The top-level SOLUS "likelihood" entry point for `fuse_property_adaptive()`.

**Parameters**:
```r
fetch_solus_percentiles(aoi_vect, solus_variable, top_depth, bottom_depth)
# aoi_vect                 - a terra::SpatVector AOI
# solus_variable           - a soilDB::fetchSOLUS()-recognized variable name
# top_depth, bottom_depth  - numeric depth window bounds in cm
```
**Returns**: `list(values = list(P025 = <low raster>, P50 = <pred raster>, P975 = <high raster>),
probs = c(0.025, 0.5, 0.975))`, or `NULL` if any of low/pred/high is unavailable.

**Behavior**: Calls `fetch_solus_low_pred_high()` and returns `NULL` if any of `pred`/`low`/`high`
is missing; otherwise packages the three rasters directly as the `P025`/`P50`/`P975` percentile
triplet (SOLUS100's published 95% prediction interval bounds plus its point prediction), matching
the `list(values=, probs=)` shape every fusion route expects.

#### 29b/30b. `fetch_solus_low_pred_high_multi()` / `fetch_solus_percentiles_multi()` — batched multi-variable SOLUS fetch (S1)

**Purpose**: Batched siblings of #29/#30 used by `run_stage1_fusion_multi()`. `soilDB::fetchSOLUS()`
accepts `variables`, `depth_slices`, and `output_type` all as vectors simultaneously
(live-verified 2026-09-04), so every needed variable's low/prediction/high rasters for one window
can be fetched in **one** `fetchSOLUS()` call instead of `3 * length(solus_variables)` separate
ones.

**Parameters**:
```r
fetch_solus_low_pred_high_multi(aoi_vect, solus_variables, top_depth, bottom_depth)
fetch_solus_percentiles_multi(aoi_vect, solus_variables, top_depth, bottom_depth)
# solus_variables          - character vector of soilDB::fetchSOLUS()-recognized variable names
# (aoi_vect, top_depth, bottom_depth same as #29/#30)
```
**Returns**: A named list, one entry per `solus_variables` element, in the same per-variable shape
`fetch_solus_low_pred_high()` / `fetch_solus_percentiles()` return (`list(pred=,low=,high=)`, or
`list(values=,probs=)` / `NULL`).

**Behavior**: Computes the depth-window weights once (shared across every variable and output
type, since the needed native slices don't depend on either), issues **one**
`soilDB::fetchSOLUS(aoi_vect, depth_slices=, variables = solus_variables, output_type = <all
three>, grid = TRUE)` call, then reduces per `(variable, output_type)` exactly as #29 does. A
single `fetchSOLUS()` error fails the whole batch (every variable's entry becomes
`list(pred=NULL,low=NULL,high=NULL)`, with a `warning()` — callers fall back to the scalar path,
as `run_stage1_fusion_multi()`'s A4 SOLUS memo does); a missing layer for just one
`(variable, output_type)` combo within an otherwise-successful response only nulls that one piece,
with its own `warning()` naming it. The scalar `fetch_solus_low_pred_high()` /
`fetch_solus_percentiles()` (#29/#30) are unchanged and remain what single-property
`run_stage1_fusion()` uses.

### Per-pixel ensemble bridge (core-fusion.R)

Wires the fused posterior back into the modelling framework as a **per-pixel** product, without
touching either existing pipeline.

#### 31. `extract_mukey_joint_ensemble()` *(adapter-ssurgo-simulate.R)*
Runs `simulate_ssurgo_mapunit_draws(depth_windows=)` once and returns, per mukey, the retained bag
of **joint** multivariate realizations as `[n_replicate x n_property]` matrices, one per depth
window, row-aligned on `(cokey, simulation_number)` so cross-property and cross-window rank
structure is preserved. The multi-property / multi-window analogue of `mukey_draws_lookup()`. Takes
`properties` (which columns to keep) and, separately, `requested_properties` (which properties the
simulation itself draws — defaulting to `properties`, so a call restricted to the Saxton-Rawls
inputs does not pay for the properties it will discard).

#### 32. `remarginalize_ensemble_to_posterior()`
Rank-preserving marginal transform: for each property x window, `v_new = F_post^-1(F_prior(v_r))`
where `F_prior` is the mukey ensemble's empirical CDF and `F_post` the pixel's fused posterior
(piecewise-linear on its 9 percentile knots, `rule = 2`). One `terra::subst()` builds the
`n_kept`-layer per-cell `u` stack; one `invert_posterior_cdf_raster()` call inverts it in a fixed
`~4*(k-1)` terra ops regardless of `n_out`. Each realization keeps its rank in every property and
window, so the ensemble's empirical copula carries over; only the marginals become the posterior.
`summarize = TRUE` -> per-pixel percentile rasters; `FALSE` -> the raw realization stacks.

**Approximations**: copula stays mukey-level; within-window vertical shape from the source
realization; comppct weighting stays mukey-level.

#### 33. `zonal_distribution_from_posterior()` *(internal)*
Percentile-only `terra::zonal()` reduction of a posterior to per-mukey `low/rep/high` (the shape
`fuse_observed_data_into_priors(observed_data_by_mukey=)` accepts). No point-estimate path. The
mukey-collapsed comparison arm for the benchmark; gated `internal` pending
`compare_perpixel_vs_zonal_awc()` numbers.

#### 34. `remarginalized_awc()`
Per-pixel available water capacity from the re-marginalized stacks, per realization, then per-pixel
percentiles, clamped at 0. Does **not** use `calculate_aws_df()` (ROSETTA network POST + own MC -
can't run per pixel).

- `method = "saxton_rawls"` (**default**): `fetchSOLUS()` publishes no water-retention variable, so
  `wr_3b`/`wr_15b` are not fusable. Instead the fusable `sand_total`/`silt_total`/`clay_total`/`db`
  (+ optional `soc`, `rfv`) are re-marginalized, then `saxton_rawls_raster()` (a `terra`-native
  vectorized port of `calculate_saxton_rawls_single()`) derives field capacity / wilting point per
  (pixel, realization). `AWC = sum_windows (fc - wp)/100 * thickness`.
- `method = "direct"`: the original `(wr_3b - wr_15b)/100 * thickness * (1 - rfv/100)` path, for a
  caller with a non-SOLUS water-retention likelihood.

`tile_rows =` processes the AOI in row-strips for bounded memory, numerically identical to the
whole-grid run.

**`restriction_depth =`** - optional
bedrock/restriction-depth truncation.** `NULL` (default) preserves every prior behavior exactly -
each window uses its full nominal `thickness = bottom - top` regardless of bedrock, which is
physically wrong for shallow soils (SOLUS predicts every property at every depth even where
bedrock is shallower) but was the only behavior available before this task. When supplied (a
single-layer `terra::SpatRaster` of restriction depth in cm - typically
`fetch_solus_restriction_depth(aoi_vect)`), resampled once (bilinear) onto the working grid, `NA`
cells recoded to `Inf` (no evidence -> don't truncate), and each window's thickness becomes the
per-pixel **effective** thickness `pmax(0, pmin(bottom, restriction_depth) - top)` instead - via
`terra::clamp(restriction_depth, upper = bottom)` then `terra::clamp(... - top, lower = 0)`, since
`clamp(x, upper=)`/`clamp(x, lower=)` are exactly `pmin`/`pmax` against a scalar. A window entirely
above the restriction keeps full thickness; one straddling it gets partial credit; one entirely
below contributes exactly **0**, not `NA` (`NA` would incorrectly blank the whole pixel's AWC,
including valid shallower windows, rather than just that one window's contribution). The
single-layer effective-thickness raster is recycled across the multi-layer `fc`/`wp` realization
stack by `terra`'s own arithmetic recycling, in place of a scalar `thickness` factor.

#### 34b. `fetch_solus_site_level()` / `fetch_solus_restriction_depth()` *(adapter-solus.R)*
**Purpose**: Fetch a SOLUS100 **site-level** (depth-independent) variable - `anylithicdpt`
(depth to bedrock) / `resdept` (depth to restriction) are the two relevant ones - requested via
`fetchSOLUS()`'s special `depth_slices = "all"`, and turn it into the censoring-guarded
`restriction_depth` raster `remarginalized_awc()` consumes.

**`fetch_solus_site_level(aoi_vect, variable, output_types = c("prediction", "95% low prediction
interval", "95% high prediction interval"))`**: one `fetchSOLUS(depth_slices = "all", ...)` call,
parses the confirmed `<variable>_all_cm_<p|l|h>` layer-naming convention (live-verified during the
S1 spike). Returns `list(pred=, low=, high=)` (deliberately not built on
`solus_depth_window_weights()`, which only applies to numeric depth slices).

**`SOLUS_RESTRICTION_CENSOR_CM <- 201`**: both `anylithicdpt` and `resdept` are right-censored at
201 cm in SOLUS100's own training data (per gSSURGO convention, confirmed against the published
SOLUS100 methodology, section 2.1.3) - a predicted value at/above this point means "no restriction
detected within the training range," not a literal depth.

**`fetch_solus_restriction_depth(aoi_vect, variable = "anylithicdpt")`**: calls
`fetch_solus_site_level()` for the prediction only, then recodes every cell
`>= SOLUS_RESTRICTION_CENSOR_CM` to `Inf` (never truncates). Default `"anylithicdpt"` - trained
(per SOLUS100's own methodology) only on `reskind %in% c("Lithic bedrock", "Paralithic bedrock")`
records, so it's hard-bedrock-specific **by construction** - no separate SSURGO-side hard/soft
classifier is needed, since SOLUS already encodes that split across its two site-level variables.
`"resdept"` is trained on *every* restriction record regardless of kind (fragipans, duripans,
petrocalcic, etc. included) - offered as a broader, more conservative (more truncation-prone)
non-default alternative.

#### 35. `saxton_rawls_raster()` *(internal)*
The vectorized `terra`-native counterpart of `calculate_saxton_rawls_single()` (`R/adapter-ssurgo-infill.R`)
- identical equations and clamps (including the `|sand+silt+clay-100| > 5` texture renormalization),
but every operation is a `terra` `Arith`/`Math`/`clamp`/`ifel`, so it runs on a whole multi-layer
realization stack in one pass. Returns `list(fc =, wp =)` volumetric % rasters.

### Statistical structure of the per-pixel bridge — and how to read the uncertainty

The bridge is a **rank-preserving re-marginalization** of a single joint Monte Carlo ensemble, not
an independent per-property resample. Understanding that is what makes the per-pixel AWC (or
per-pixel property percentiles) interpretable as a genuine probability distribution.

**1. One joint source ensemble.** `extract_mukey_joint_ensemble()` runs the SSURGO Monte Carlo
**once**. `simulate_cokey_generalized()` draws every property jointly under the KSSL correlation
matrices; the joint-copula vertical-correlation step runs down each simulated profile; the result is
aggregated to the requested depth windows **row-aligned** - realization `r` in window `0-5` and
realization `r` in window `5-15` are the *same* simulated profile. So the source ensemble carries
(a) cross-property correlation within each window and (b) vertical correlation of each property
across windows.

**2. A single realization index threads through the transform.**
`remarginalize_ensemble_to_posterior()` computes, for realization `r`, property `p`, window `w`:
`u = (rank of r in the source ensemble for (p, w) - 0.5) / n`, then `v_new = F_posterior^-1(u)` at
each pixel. The **same index `r`** is used for every property and every window (`n_kept` is one
count across all mukeys and windows). Realization `r` therefore keeps its rank position everywhere,
so the transformed ensemble preserves the source's **Spearman** cross-property correlation *within*
a window and vertical correlation *across* windows. Only the marginals change - each becomes that
pixel's SSURGO x SOLUS fused posterior. (Pinned by a unit test: transformed-ensemble Spearman `cor`
between two properties equals the source's to `1e-6`.)

**3. Derived quantities are computed per realization, then combined.**
`remarginalized_awc(method = "saxton_rawls")` takes, for realization `r` / pixel `c` / window `w`,
the jointly-consistent `(sand_r, clay_r, silt_r, db_r, om_r, rfv_r)` -> `saxton_rawls_raster()` ->
`(fc_r, wp_r)`, then `AWC_r = sum_windows (fc_r - wp_r)/100 * thickness_w` (summed over windows for
the **same** realization, so shallow and deep contributions are vertically consistent). The result
is an `n_kept`-layer AWC raster; `terra::quantile()` per pixel gives the AWC probability distribution
(P5..P95) -> the per-pixel uncertainty band.

**What the uncertainty band does and does not include:**

| Aspect | Treatment |
|---|---|
| Cross-property & cross-depth **dependence** | The SSURGO tabular copula (KSSL + joint-copula vertical correlation), **held fixed**. SOLUS carries no joint information - only independent per-property, per-depth point predictions - so there is no fused joint structure to use. |
| Correlation type preserved | **Rank (Spearman)**, not linear (Pearson). "Wet stays wet", "high-clay stays high-clay" is exact; Pearson structure is approximate if a marginal shifts substantially under fusion. |
| **Marginal** uncertainty of each input property | The full SSURGO x SOLUS fused posterior spread, per pixel, propagated through Saxton-Rawls. |
| Pedotransfer-function error | **Not propagated.** Saxton-Rawls is applied deterministically per realization; the +/-15% `_l`/`_h` band `calculate_saxton_rawls_single()` returns is not used. |
| Sub-mukey variation in the **dependence** structure | None - every pixel in a mukey shares the source ensemble (same copula, same rank ordering); only the fused marginals vary pixel-to-pixel. |
| Texture compositional closure | Independent re-marginalization can leave `sand + silt + clay != 100`; `saxton_rawls_raster()` renormalizes to 100 when off by more than 5 (the equations use only the sand and clay fractions). |

In short: the per-pixel AWC distribution reflects **SOLUS-fused marginal uncertainty threaded
through the SSURGO tabular copula**, not a fully re-estimated joint posterior and not
pedotransfer-model error.

## Internal Connections

```
                                   AOI (terra::SpatVector, e.g. EPSG:5070)
                                                  |
                    +-----------------------------+-----------------------------+
                    |                                                           |
                    v                                                           v
   SSURGO adapter (adapter-ssurgo-simulate.R)                        SOLUS100 adapter (adapter-solus.R)
   fetch_ssurgo_mukey_raster()                                 closest_solus_depth_slice()
     -> soilDB::mukey.wcs() + soilDB::SDA_spatialQuery()        fetch_solus_low_pred_high()
   simulate_ssurgo_mapunit_draws()                                -> soilDB::fetchSOLUS() x3
     -> download_ssurgo_tabular() [R/adapter-ssurgo-acquire.R, cached]    (low / pred / high)
     -> infill_soil_data() -> classify_genhz()                 fetch_solus_percentiles()
     -> sim_component_comp() [R/core-simulation.R]            -> list(values=list(P025,P50,P975),
     -> simulate_cokey_generalized() [property-simulation.R]           probs=c(.025,.5,.975))
     -> maybe_adjust_soil_data_depth_trend()
     -> aggregate_depth_window_by_replicate()
   fetch_ssurgo_percentiles()
     -> rasterize_mukey_percentiles()
     -> list(values=<percentile SpatRasters>, probs=)
                    |                                                           |
                    v                                                           v
        prior percentile-value rasters                         likelihood percentile-value rasters
             (SSURGO grid, ~30 m)                                       (SOLUS grid, 100 m)
                    |                                                           |
                    +--------------------- terra::resample() ------------------+
                              (prior bilinear-resampled onto the SOLUS grid,
                               inside run_stage1_fusion()/run_stage1_fusion_group())
                                                  |
                                                  v
                             fuse_property_adaptive()  [raster-fusion.R - GENERIC CORE]
                             ├── align_percentile_probs()
                             ├── resolve_property_dist()
                             └── dispatch on dist:
                                   ├── normal/beta/gamma -> fuse_adaptive()
                                   │      ├── fuse_closed_form()  (large AOI)
                                   │      │     uses core-distributions-raster.R's
                                   │      │     fit_normal_raster()/fit_beta_mle_newton_raster()/
                                   │      │     fit_gamma_mom_raster(), and core-fusion.R's
                                   │      │     bayes_update_normal_normal()/fuse_beta()/fuse_gamma()
                                   │      └── fuse_general_kde()   (small AOI)
                                   │            uses simulate_from_percentiles()/bayesian_update()
                                   ├── lognormal -> fuse_lognormal_adaptive()
                                   └── metalog   -> fuse_metalog_adapter()
                                          uses fit_metalog_linear_raster()/
                                          check_metalog_feasibility_raster()/
                                          metalog_moments_raster()
                                                  |
                                                  v
                                  fused posterior SpatRaster(s)
                            (mu/sigma, alpha/beta, or shape/rate, per resolved `dist`)

  Compositional (texture) path: run_stage1_fusion_group() fetches clay/sand/silt priors+
  likelihoods per member (same adapters as above), then fuse_texture_group() fuses all three
  JOINTLY via ILR (estimate_ilr_moments_mc()/fuse_bivariate_normal()/ilr_inverse()) instead of
  three independent fuse_property_adaptive() calls, to preserve sum-to-100.

  cache.R (build_cache_key() / cache_get() / cache_set() / wrap_percentile_list() /
  unwrap_percentile_list()) sits ALONGSIDE both adapters and both top-level orchestrators,
  disk-caching SSURGO tabular data, SSURGO/SOLUS percentile rasters, and texture-group/per-
  property posterior results under tools::R_user_dir("soilSIM", "cache"), keyed by
  AOI + property/group id + depth window + kind.


  PER-PIXEL ENSEMBLE BRIDGE (core-fusion.R) - optional, purely additive layer on top:

  extract_mukey_joint_ensemble(aoi, depth_windows)          run_stage1_fusion_multi(aoi, configs,
    -> simulate_ssurgo_mapunit_draws(depth_windows=)          depth_windows)  -- ONE simulation for
         (ONE joint Monte Carlo; KSSL + joint-copula           every property x window, seeding the
          vertical correlation; row-aligned across windows)     same per-property caches
    -> per mukey: list(windows = [n_replicate x n_property]                     |
                       matrices, replicate_key)                                v
                                                              posterior$percentiles rasters
                                                               (marginal, per property x window)
                    \                                           /
                     \                                         /
                      v                                       v
              remarginalize_ensemble_to_posterior(ensemble, posterior_by_property_window)
                - u[r] = rank(r in source ensemble)/n         (ONE realization index r for
                - v_new[cell, r] = F_post^-1(u[r])             every property AND window)
                - terra::subst() -> n_kept-layer u stack; invert_posterior_cdf_raster()
                => transformed ensemble: fused marginals, SOURCE Spearman copula preserved
                                                  |
                    +-----------------------------+------------------------------+
                    v                                                            v
        summarize = TRUE:                                        remarginalized_awc(method=)
        per-pixel percentile rasters                               "saxton_rawls" (default):
        [[property]][[window]]                                       re-marg sand/silt/clay/db
                                                                     -> saxton_rawls_raster()
        zonal_distribution_from_posterior()  (internal)              -> fc/wp per (pixel, realization)
        cheap mukey-collapsed comparison arm ->                     "direct": wr_3b/wr_15b path
        fuse_observed_data_into_priors(observed_data_by_mukey=)      => AWC_r = sum_w (fc-wp)/100*thk
                                                                     => terra::quantile() per pixel
                                                                        = per-pixel AWC distribution
```

## Dependencies

**External packages**:
- `terra` — used throughout every file in this pipeline (`SpatRaster`/`SpatVector` arithmetic,
  `terra::app()`/`terra::lapp()`/`terra::ifel()`/`terra::rasterize()`/`terra::resample()`/
  `terra::catalyze()`/`terra::wrap()`/`terra::unwrap()`).
- `soilDB::fetchSOLUS()` — the live network API `adapter-solus.R`'s
  `fetch_solus_low_pred_high()` calls for SOLUS100 low/prediction/high rasters.
- `soilDB::mukey.wcs()` / `soilDB::SDA_spatialQuery()` — the live network APIs
  `adapter-ssurgo-simulate.R`'s `fetch_ssurgo_mukey_raster()` calls for SSURGO map unit grids/polygons.
- `digest` — `cache.R`'s `build_cache_key()` hashes the AOI's WKT via `digest::digest()`.
- `dplyr` — used throughout the SSURGO adapter for grouping/joining/summarizing
  (`dplyr::group_by()`/`left_join()`/`summarise()`/`group_modify()`).
- `GPfit` (optional) — guarded by `requireNamespace()` in
  `maybe_adjust_soil_data_depth_trend()`; skipped with a warning if unavailable.

**soilSIM dependencies**:
- `R/core-distributions.R` — `core-distributions-raster.R` directly reuses
  `metalog_basis_matrix()`/`metalog_to_z()`/`metalog_from_z()` unchanged, since they're pure
  elementwise arithmetic that already works on `SpatRaster` inputs (confirmed by a dedicated smoke
  test before this file was written, not merely assumed); `core-fusion.R`'s `fuse_texture_group()`
  similarly reuses `estimate_ilr_moments_mc()`/`ilr_inverse()` directly, and
  `resolve_composition_groups()`'s composition-group config convention is reused rather than
  duplicated.
- `R/core-fusion.R` — `core-fusion.R` reuses `bayes_update_normal_normal()`, `fuse_beta()`,
  `fuse_gamma()`, `moments_to_gamma()`/`moments_to_beta()`/`beta_to_moments()`/`gamma_to_moments()`,
  `normal_to_lognormal_params()`/`lognormal_to_normal_params()`, `bayesian_update()`, and
  `fuse_bivariate_normal()` directly, unmodified — all are pure elementwise arithmetic, so they
  already work unchanged on `SpatRaster` inputs exactly as they do on plain numerics.
- `R/core-simulation.R` — feeds the SSURGO adapter's Monte Carlo step: `sim_component_comp()`
  and `simulate_cokey_generalized()` are called directly from
  `simulate_ssurgo_mapunit_draws()`.
- `R/adapter-ssurgo-acquire.R` — `download_ssurgo_tabular()` is called by
  `simulate_ssurgo_mapunit_draws()` for the underlying tabular SSURGO fetch (cached separately from
  the raster-fusion cache's own `"ssurgo"`/`"solus"` kinds, under a `"ssurgo_tabular"` kind, keyed
  by AOI alone - see `ssurgo_tabular_cache_key()`, `R/cache.R`).
- `R/adapter-ssurgo-infill.R` — `infill_soil_property()` is called by `infill_soil_data()`.
- `R/core-correlations.R` — `.kssl_property_matrices()`/`.kssl_texture_matrices()`
  (already built into `R/sysdata.rda`) and `classify_genhz()` are called by
  `simulate_ssurgo_mapunit_draws()`.
- `R/core-gp.R` — `apply_local_gp_adjustments()` is called by
  `maybe_adjust_soil_data_depth_trend()`.
- `R/core-montecarlo.R` — `simulate_from_percentiles()` is called (per cell, inside `terra::app()`) by
  `fuse_general_kde()`. The per-pixel bridge additionally reuses `fuse_observed_data_into_priors()`'s
  `observed_data_by_mukey=` parameter (the A.1 opt-in) as the sink for `zonal_distribution_from_posterior()`.
- `R/adapter-ssurgo-infill.R` — the per-pixel bridge's `saxton_rawls_raster()` is a `terra`-native port of
  this file's `calculate_saxton_rawls_single()` (the equations are duplicated, not called, because
  the scalar function's `max()`/`min()`/`if` do not vectorize over a raster stack).

`run_stage1_fusion()`/`run_stage1_fusion_group()` (`R/core-fusion.R`) are the top-level AOI
orchestrators tying all of the above together — the single entry points a downstream caller
needs, hiding the fetch/cache/align/fuse sequence behind one call per property (or
per compositional group) per AOI/depth window.

## Data Flow In/Out

**In**:
- An AOI, as a `terra::SpatVector` (projected, e.g. EPSG:5070).
- A property config (`property_config`): `id`, `solus_variable`, `dist` (`"normal"`/`"beta"`/
  `"gamma"`/`"lognormal"`/`"metalog"`/`"auto"`), and optionally `bounds`/`boundedness`/
  `auto_skew_threshold`/`composition_group`.
- Depth window bounds (`top_depth`/`bottom_depth`, cm).
- For compositional properties: a composition-group config (`composition_groups`, in
  `config$monte_carlo$composition_groups` shape) and a per-member `property_configs` list.

**Out**:
- For a single, non-compositional property: a fused posterior `SpatRaster` set (family-native
  parameter rasters — `mu`/`sigma`, `alpha`/`beta`, or `shape`/`rate` — per the resolved `dist`),
  plus the raw prior/likelihood rasters, routing diagnostics (`route`, `route_detail`,
  `n_fallback_cells`), and distribution-resolution diagnostics (`dist`, `dist_source`,
  `skew_proxy`) — from `run_stage1_fusion()`.
- For a compositional group: each member's fused texture posterior (`value` point-estimate raster
  plus shared `ilr_mu`/`ilr_Sigma` rasters for posterior sampling) — from `run_stage1_fusion_group()`.
- Intermediate cached artifacts on disk (under `tools::R_user_dir("soilSIM", "cache")`): raw SSURGO
  tabular data, SSURGO/SOLUS percentile rasters (`terra::wrap()`ped), texture-group fusion results,
  and per-property posterior results, each keyed by AOI + id + depth window + kind via
  `build_cache_key()`.

## Known Limitations

1. **`metalog_moments_raster()` is new, less-validated glue code** (`R/core-fusion.R`). Unlike the
   rest of this pipeline's math — validated against `fitdistrplus`/closed-form references —
   the mean/variance-via-quadrature computation in `metalog_moments_raster()` is less validated
   and should be spot-checked against real data before being trusted in production.

2. **`SpatRaster`/`SpatVector` objects do not survive a plain `saveRDS()`/`readRDS()` round-trip**
   (`R/cache.R`). These objects hold an external pointer to in-memory/on-disk GDAL state.
   `cache_get()`/`cache_set()` themselves do not wrap/unwrap automatically, since not every cached
   value is a raster (e.g. Monte Carlo draw data frames aren't). Callers caching raster values must
   call `terra::wrap()` before `cache_set()` and `terra::unwrap()` after `cache_get()` — which is
   exactly what `wrap_percentile_list()`/`unwrap_percentile_list()` do for the
   `list(values=<SpatRasters>, probs=)` percentile shape, and what `run_stage1_fusion()`/
   `run_stage1_fusion_group()` call around every SSURGO/SOLUS cache read and write.

3. **Per-pixel bridge — dependence structure is not re-estimated from SOLUS.** The bridge fuses only
   the *marginals*; the cross-property and cross-depth dependence carried into the per-pixel product
   is the SSURGO tabular copula (KSSL + joint-copula vertical correlation), held fixed. SOLUS
   provides independent per-property, per-depth point predictions with no joint information, so
   there is nothing fused to use. See "Statistical structure of the per-pixel bridge" above.

4. **Per-pixel bridge — Saxton-Rawls pedotransfer error is not propagated.** In
   `remarginalized_awc(method = "saxton_rawls")` the pedotransfer is applied deterministically per
   realization; the +/-15% band `calculate_saxton_rawls_single()` returns is unused. The per-pixel
   AWC distribution reflects fused-posterior *input* uncertainty only.

5. **`remarginalized_awc()` cannot fuse water retention directly.** `fetchSOLUS()` publishes no
   water-retention variable, so `wr_3b`/`wr_15b` have no SOLUS likelihood. `method = "saxton_rawls"`
   (default) derives them from the fusable sand/silt/clay/`dbovendry`/`soc`/`fragvol`;
   `method = "direct"` is only usable with a non-SOLUS water-retention likelihood.

## Usage Example

```r
library(terra)

# AOI, already projected (e.g. EPSG:5070)
aoi_vect <- vect("aoi_polygon.gpkg")

# --- Using the top-level orchestrator (fetches + caches both sides internally) ---
result <- run_stage1_fusion(
  aoi_vect = aoi_vect,
  property_config = list(id = "clay", solus_variable = "claytotal", dist = "auto"),
  top_depth = 0, bottom_depth = 30
)
if (!is.null(result)) {
  posterior_mu    <- result$posterior$mu     # or $alpha/$beta, or $shape/$rate, per result$dist
  posterior_sigma <- result$posterior$sigma
  writeRaster(posterior_mu, "clay_posterior_mean_0_30cm.tif", overwrite = TRUE)
}

# --- Using fuse_property_adaptive() directly with pre-fetched percentile rasters ---
ssurgo_prior <- fetch_ssurgo_percentiles(aoi_vect, "clay", top_depth = 0, bottom_depth = 30)
solus_lik    <- fetch_solus_percentiles(aoi_vect, "claytotal", top_depth = 0, bottom_depth = 30)

# Align the SSURGO (~30m) prior onto the SOLUS (100m) grid before fusing.
prior_aligned <- lapply(ssurgo_prior$values, resample, y = solus_lik$values[[1]], method = "bilinear")

fused <- fuse_property_adaptive(
  prior_value_rasters = prior_aligned, prior_probs = ssurgo_prior$probs,
  lik_value_rasters   = solus_lik$values, lik_probs = solus_lik$probs,
  property_config     = list(dist = "beta", bounds = c(0, 100))
)

cat("route:", fused$route, "\n")
plot(fused$posterior$alpha)

# --- Per-pixel joint ensemble + per-pixel AWC (core-fusion.R) ---
windows <- list(c(0, 5), c(5, 15), c(15, 30))

# 1. one joint SSURGO Monte Carlo, retained per mukey, aligned across windows
ens <- extract_mukey_joint_ensemble(aoi_vect, windows)

# 2. a SOLUS-fused posterior for each Saxton-Rawls input, per property x window - one simulation
#    for all 5 x 3 leaves, not 15
sr_solus <- c(sand_total = "sandtotal", silt_total = "silttotal",
              clay_total = "claytotal", db = "dbovendry", soc = "soc")
configs <- setNames(
  lapply(names(sr_solus), function(nm) list(id = nm, solus_variable = sr_solus[[nm]], dist = "auto")),
  names(sr_solus)
)
pbpw <- run_stage1_fusion_multi(aoi_vect, configs, windows, simplify = TRUE)

# 3. per-pixel AWC probability distribution (P5..P95 rasters), clamped at 0
awc <- remarginalized_awc(ens, pbpw, n_out = 250)         # method = "saxton_rawls" by default
plot(awc$awc_cm$P50)                                       # median AWC (cm) over 0-30 cm
plot(awc$awc_cm$P95 - awc$awc_cm$P5)                       # 90% credible-interval width = uncertainty
```
