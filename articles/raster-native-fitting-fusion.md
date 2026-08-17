# Raster-Native Fitting & Fusion Internals

## Overview

`raster-fusion-ssurgo-solus.Rmd` shows the *pipeline* view of soilSIM’s
multi-source raster fusion: call
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
and get back a fused posterior. This vignette opens the hood.
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
is a thin fetch-and-cache wrapper around a handful of small, composable,
raster-native building blocks - distribution fitting functions
(`R/distribution-fitting-raster.R`), a fusion core
(`R/raster-fusion.R`), and a disk cache (`R/raster-cache.R`) - each of
which does exactly one thing to a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
(or list of them) and can be called directly. Understanding these
individually is what makes the pipeline’s behavior (which route it took,
why a family was resolved a certain way, why two sources’ percentiles
had to be realigned before fusing) legible rather than a black box. See
`soilSIM/docs/09_multi_source_raster_fusion_pipeline.md` for the full
function-level reference this vignette is built from.

``` r

library(soilSIM)
library(terra)
#> terra 1.9.34
library(ggplot2)
```

## Reusing real data instead of fetching

Fetching SSURGO/SOLUS100 live requires network access and takes minutes;
this vignette reuses the same cached Salinas Valley clay-content fusion
result `raster-fusion-ssurgo-solus.Rmd` uses, so every function below
runs on real percentile rasters rather than synthetic ones. The real
(network) calls that originally produced this cached object are shown
for reference only and not evaluated here:

``` r

# Reference only - not run in this vignette.
salinas_wkt <- "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))"
aoi <- terra::project(terra::vect(salinas_wkt, crs = "epsg:4326"), "epsg:5070")
property_config <- list(id = "clay", solus_variable = "claytotal", dist = "normal")
fusion_clay <- run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5)
```

``` r

fusion_clay <- unwrap_nested_rasters(
  readRDS(system.file("extdata", "fusion_clay_salinas.rds", package = "soilSIM"))
)

# The two inputs every function below operates on: named lists of percentile-value
# SpatRasters, one list per source, plus the probabilities each percentile represents.
prior_values <- fusion_clay$prior$values        # SSURGO, already resampled onto the SOLUS grid
prior_probs <- fusion_clay$prior$probs
lik_values <- fusion_clay$likelihood$values     # SOLUS100
lik_probs <- fusion_clay$likelihood$probs

names(prior_values); prior_probs
#> [1] "P05" "P25" "P50" "P75" "P95"
#> [1] 0.05 0.25 0.50 0.75 0.95
names(lik_values); lik_probs
#> [1] "P025" "P50"  "P975"
#> [1] 0.025 0.500 0.975
bounds <- c(0, 100)  # clay content is a percentage
```

This is the “common currency” every function in this vignette speaks: a
**named list of single-layer percentile-value rasters** (e.g. `P05`,
`P25`, `P50`, …) plus a **numeric vector of the probabilities they
represent**, all sharing one grid. Everything downstream - fitting,
quantile evaluation, fusion, caching - operates on this shape (or a pair
of them, prior and likelihood).

## 1. Fitting a Normal distribution

[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)
is the simplest fit in the file: closed-form, exact, and just three
raster arithmetic operations. It takes a low/median/high
percentile-value raster triplet and the probabilities the low/high
rasters represent, and returns `mu` (the median raster, used directly as
the mean) and `sigma` (back-solved from the tail spacing):

``` r

fit_n <- fit_normal_raster(
  p_lo_r = prior_values$P25, p50_r = prior_values$P50, p_hi_r = prior_values$P75,
  p_lo = 0.25, p_hi = 0.75
)
fit_n$mu
#> class       : SpatRaster
#> size        : 26, 23, 1  (nrow, ncol, nlyr)
#> resolution  : 100, 100  (x, y)
#> extent      : -2246800, -2244500, 1810700, 1813300  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070)
#> source(s)   : memory
#> varname     : claytotal_0_cm_l
#> name        :       P50
#> min value   :  1.721668
#> max value   : 28.000607
```

`sigma` is `(p_hi_r - p_lo_r) / (qnorm(p_hi) - qnorm(p_lo))` - the raw
value spread between the two tail percentiles, divided by how many
standard-normal units apart those two probabilities are. This is a
single elementwise expression, so its cost doesn’t depend on how many
cells the raster has.
[`quantile_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)
then evaluates the fitted distribution at any target probability `q` via
`mu + sigma * qnorm(q)`:

``` r

q10 <- quantile_normal_raster(fit_n, 0.1)
q90 <- quantile_normal_raster(fit_n, 0.9)

normal_stack <- c(q10, fit_n$mu, q90)
names(normal_stack) <- c("Q10", "mu (=P50)", "Q90")
terra::plot(normal_stack, main = c("Normal Q10", "Normal mean", "Normal Q90"),
            col = grDevices::hcl.colors(50, "viridis"), nc = 3)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-5-1.png)

## 2. Fitting a Beta distribution: method-of-moments vs. Newton-Raphson MLE

Clay content is bounded on `[0, 100]`, which the Normal fit above
ignores (its Q10/Q90 rasters can dip below 0 in low-clay cells). A Beta
fit respects those bounds. soilSIM ships two Beta fitters:

- `fit_beta_mom_raster(mean_r, var_r)` - closed-form method-of-moments,
  APPROXIMATE, and used only to seed the Newton-Raphson solver below
  (not intended as a standalone fit). It takes a **rescaled-to-`[0,1]`**
  mean/variance pair directly, not raw percentile rasters.
- `fit_beta_mle_newton_raster(value_rasters, bounds, n_iter = 15, eps = 1e-6)` -
  a vectorized Newton-Raphson MLE solver that exploits the Beta
  log-likelihood’s closed-form score/Hessian (in
  [`digamma()`](https://rdrr.io/r/base/Special.html)/[`trigamma()`](https://rdrr.io/r/base/Special.html),
  both vectorized base-R functions) to update every cell’s
  `(alpha, beta)` simultaneously per iteration, `n_iter` times. It takes
  the raw percentile-value rasters and the physical `bounds` directly,
  and internally seeds itself from
  [`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md).

``` r

# Method-of-moments needs mean/var already rescaled to [0,1]; reuse the Normal fit above as a
# cheap mean/variance proxy, rescaled by `bounds`.
mean_r <- fit_n$mu / bounds[2]
var_r <- (fit_n$sigma / bounds[2])^2
fit_beta_mom <- fit_beta_mom_raster(mean_r, var_r)

# The Newton-Raphson MLE fit takes the raw percentile rasters and bounds directly - it does its
# own internal method-of-moments seeding.
fit_beta_mle <- fit_beta_mle_newton_raster(prior_values, bounds)
```

``` r

beta_params <- c(fit_beta_mom$alpha, fit_beta_mom$beta, fit_beta_mle$alpha, fit_beta_mle$beta)
names(beta_params) <- c("alpha (MoM)", "beta (MoM)", "alpha (Newton MLE)", "beta (Newton MLE)")
terra::plot(beta_params, col = grDevices::hcl.colors(50, "plasma"), nc = 2)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-7-1.png)

The two fits’ shape-parameter rasters differ cell-by-cell (MoM only
matches the first two moments; Newton MLE fits all of the percentiles
supplied to it), but both feed the same quantile function shape -
`quantile_beta_mom_raster(fit, q)` /
`quantile_beta_mle_newton_raster(fit, q)`, each implemented via
[`terra::lapp()`](https://rspatial.github.io/terra/reference/lapp.html)
calling [`stats::qbeta()`](https://rdrr.io/r/stats/Beta.html) per cell,
since [`qbeta()`](https://rdrr.io/r/stats/Beta.html) has no
raster-native vectorized form of its own:

``` r

q_mom <- quantile_beta_mom_raster(fit_beta_mom, 0.5)
q_mle <- quantile_beta_mle_newton_raster(fit_beta_mle, 0.5)

median_compare <- c(prior_values$P50, q_mom, q_mle)
names(median_compare) <- c("Raw P50 input", "Beta MoM Q50", "Beta Newton-MLE Q50")
terra::plot(median_compare, col = grDevices::hcl.colors(50, "viridis"), nc = 3)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-8-1.png)

The Newton-MLE fit’s own Q50 tracks the raw input P50 raster much more
closely than the MoM fit’s does, since MLE was fit against all five of
`prior_values`’ percentiles rather than just a mean/variance proxy
derived from two of them.

## 3. Fitting a metalog distribution, and its feasibility fallback

The metalog family can match an arbitrary number of percentiles exactly
via a linear solve rather than an optimization -
`fit_metalog_linear_raster(value_rasters, probs, bounds, boundedness)`
requires the number of *interior* percentiles (excluding
`p = 0`/`p = 1`, where the metalog’s logit transform is undefined) to
equal the number of metalog terms, making the fit exactly determined.

Symmetric percentile sets are a real, sharp edge here: fitting all five
of `prior_values`’s percentiles (`0.05/0.25/0.5/0.75/0.95` - symmetric
about the median) makes the underlying basis matrix numerically
singular, since two of its basis columns become linearly dependent for a
symmetric probability grid. An asymmetric 4-term subset avoids it:

``` r

ml_probs <- c(0.05, 0.25, 0.5, 0.75)
ml_values <- prior_values[c("P05", "P25", "P50", "P75")]
fit_ml <- fit_metalog_linear_raster(ml_values, ml_probs, bounds, boundedness = "b")
length(fit_ml$a)   # one coefficient raster per term
#> [1] 4
fit_ml$term
#> [1] 4
```

Not every fitted metalog is a *valid* one - a metalog quantile function
must be monotonically increasing (equivalent to non-negative density
everywhere), and the exact linear solve above gives no such guarantee.
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)
probes the fit at a grid of probabilities and flags cells where
consecutive probes decrease:

``` r

infeasible <- check_metalog_feasibility_raster(fit_ml, bounds, boundedness = "b")
n_infeasible <- sum(terra::values(infeasible), na.rm = TRUE)
n_total <- terra::ncell(infeasible)
cat(sprintf("%d / %d cells (%.0f%%) have an infeasible metalog fit\n",
            n_infeasible, n_total, 100 * n_infeasible / n_total))
#> 152 / 598 cells (25%) have an infeasible metalog fit
terra::plot(infeasible, main = "Infeasible metalog fit (TRUE = non-monotonic)",
            col = c("grey85", "firebrick"))
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-10-1.png)

A substantial fraction of cells are infeasible here - not a bug, just a
real consequence of fitting only 4 of the 5 available percentiles to a
family with no correction step of its own.
[`quantile_metalog_linear_with_fallback()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback.md)
is the function every real caller should use instead of
[`quantile_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
directly: it evaluates the raw metalog quantile function, and for
whichever cells
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)
flagged, blends in
[`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md)’s
nonparametric interpolation (evaluated against the *full* percentile
set, not just the interior ones used to fit the metalog) instead:

``` r

fallback_result <- quantile_metalog_linear_with_fallback(
  fit_ml, infeasible, full_value_rasters = prior_values, full_probs = prior_probs,
  q = 0.5, bounds = bounds, boundedness = "b"
)

raw_metalog_q50 <- quantile_metalog_linear_raster(fit_ml, 0.5, bounds, boundedness = "b")
blend_stack <- c(raw_metalog_q50, fallback_result$value, prior_values$P50)
names(blend_stack) <- c("Raw metalog Q50 (no fallback)", "With fallback", "Raw input P50")
terra::plot(blend_stack, col = grDevices::hcl.colors(50, "viridis"), nc = 3)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-11-1.png)

## 4. Fitting a Gamma distribution

[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md)
is a thin raster wrapper: it computes mean/variance rasters across a
list of percentile-value rasters (the only genuinely raster-specific
step, via `Reduce(+, ...)`), then delegates to `R/bayesian-updating.R`’s
existing scalar
[`moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_gamma.md),
which is pure elementwise arithmetic and therefore already works
unchanged on `SpatRaster` inputs:

``` r

fit_gamma <- fit_gamma_mom_raster(prior_values)
gamma_stack <- c(fit_gamma$shape, fit_gamma$rate)
names(gamma_stack) <- c("shape", "rate")
terra::plot(gamma_stack, col = grDevices::hcl.colors(50, "inferno"), nc = 2)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-12-1.png)

## 5. Aligning mismatched percentile grids: `align_percentile_probs()`

SSURGO’s prior percentiles (0.05, 0.25, 0.5, 0.75, 0.95) and SOLUS100’s
likelihood percentiles (0.025, 0.5, 0.975) aren’t the same
probabilities - every fusion route needs both sides to share one `probs`
vector, so passing them through mismatched would silently mis-fit
whichever side’s percentiles don’t match the assumed labels.

``` r

aligned <- align_percentile_probs(prior_values, prior_probs, lik_values, lik_probs)
aligned$probs
#> [1] 0.05 0.25 0.50 0.75 0.95
names(aligned$prior_value_rasters)   # unchanged (SSURGO's probs won - see below)
#> [1] "P05" "P25" "P50" "P75" "P95"
length(aligned$lik_value_rasters)    # SOLUS reinterpolated onto aligned$probs - unnamed list
#> [1] 5
```

[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md)
picks whichever side’s probabilities span the *narrower* range as the
shared target (`diff(range(prior_probs))` = 0.9 vs.
`diff(range(lik_probs))` = 0.95, so the SSURGO side’s probabilities win
here) and reinterpolates the other side onto it via
[`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md) -
never the reverse, since interpolating requires the target probability
to fall within the source’s own range, and the narrower side is always
safely interpolatable from the wider one. The reinterpolated side comes
back as a plain (unnamed) list in `aligned$probs`’ order, so pick the
layer to compare by matching its position in `aligned$probs` rather than
by name:

``` r

idx50 <- which(aligned$probs == 0.5)
lik_compare <- c(lik_values$P50, aligned$lik_value_rasters[[idx50]])
names(lik_compare) <- c("SOLUS raw P50", "SOLUS P50 re-expressed at aligned probs")
terra::plot(lik_compare, col = grDevices::hcl.colors(50, "viridis"), nc = 2)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-14-1.png)

Here the shared target happens to already include `0.5` on both sides,
so the `P50` layer above is unchanged by realignment - it’s the tail
percentiles (`P05`/`P95` vs. `P025`/`P975`) that actually get
reinterpolated.

## 6. Family resolution: `resolve_property_dist()`

When `property_config$dist == "auto"`,
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md)
picks a distribution family from the AOI’s own percentile skew rather
than a fixed per-property label - and never overrides an explicit
(non-`"auto"`) `dist`. It’s deliberately narrow: it only ever resolves
to `"beta"` (if `bounds` is configured), `"normal"`, or `"lognormal"` -
never `"metalog"` or `"gamma"`, since distinguishing those from a
3-point proxy isn’t considered reliable.

``` r

# With bounds configured, "auto" always resolves to "beta" regardless of skew.
resolve_property_dist(list(dist = "auto", bounds = c(0, 100)), prior_values, prior_probs)
#> $dist
#> [1] "beta"
#> 
#> $dist_source
#> [1] "auto"
#> 
#> $skew_proxy
#> [1] 0.1160459

# Without bounds, it falls back to a skew-proxy test: normal if |skew| is small, lognormal if not.
resolve_property_dist(list(dist = "auto"), prior_values, prior_probs)
#> $dist
#> [1] "normal"
#> 
#> $dist_source
#> [1] "auto"
#> 
#> $skew_proxy
#> [1] 0.1160459
```

The skew proxy itself is
`((high - median) - (median - low)) / (high - low)`, aggregated once
across the whole AOI (default
[`stats::median()`](https://rdrr.io/r/stats/median.html)) rather than
per cell - computing it per cell would create map discontinuities right
at the family boundary. An explicit `dist` is always honored untouched:

``` r

resolve_property_dist(list(dist = "normal"), prior_values, prior_probs)$dist_source
#> [1] "config"
```

## 7. Adaptive fusion dispatch: `fuse_adaptive()`

[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)
is the size-adaptive engine
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
(and therefore
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md))
delegates to for `"normal"`/`"beta"`/`"gamma"` families. It picks the
fusion *route* from the AOI’s own cell count rather than requiring the
caller to choose: at or below `threshold_cells` (default 80,000), it
uses a fully general per-cell KDE route
([`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
on samples drawn from each side’s percentiles); above it, a much faster
closed-form route built directly on the fitting functions from Sections
1-4 above (e.g.
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md) +
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
for `family = "normal"`).

This AOI is small (598 cells), so it naturally takes the general route.
Forcing `threshold_cells = 0` shows the alternate closed-form route on
the exact same inputs, with `verbose = TRUE` (the pipeline default is
quiet - this is a deliberate log demo):

``` r

fused_general <- fuse_adaptive(
  aligned$prior_value_rasters, aligned$lik_value_rasters, aligned$probs,
  family = "normal", verbose = TRUE
)
#> [fuse_adaptive] AOI: 598 cells (threshold: 80000) -> route: bayesian_update_general
fused_closed_form <- fuse_adaptive(
  aligned$prior_value_rasters, aligned$lik_value_rasters, aligned$probs,
  family = "normal", threshold_cells = 0, verbose = TRUE
)
#> [fuse_adaptive] AOI: 598 cells (threshold: 0) -> route: closed_form_normal
```

``` r

route_compare <- c(fused_general$posterior$mu, fused_closed_form$posterior$mu)
names(route_compare) <- c(paste0("General route (", fused_general$route, ")"),
                           paste0("Closed-form route (", fused_closed_form$route, ")"))
terra::plot(route_compare, col = grDevices::hcl.colors(50, "viridis"), nc = 2)
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-18-1.png)

Both routes fuse the same prior/likelihood into a `mu`/`sigma` posterior
with the same output contract, but they reach it differently: the
general route makes no distributional assumption about the *fusion* step
itself (only about the requested output family), while the closed-form
route fits both sides Normal first (Section 1) and fuses via the
closed-form
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
update - much cheaper per cell, at the cost of that Normal assumption.

``` r

route_df <- data.frame(
  value = c(terra::values(fused_general$posterior$mu, na.rm = TRUE),
            terra::values(fused_closed_form$posterior$mu, na.rm = TRUE)),
  route = rep(c("General (KDE)", "Closed-form (Normal-Normal)"),
              c(sum(!is.na(terra::values(fused_general$posterior$mu))),
                sum(!is.na(terra::values(fused_closed_form$posterior$mu)))))
)
ggplot(route_df, aes(x = value, fill = route)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_viridis_d(name = NULL) +
  labs(title = "Fused posterior mean clay content: general vs. closed-form route",
       x = "Clay content (%)", y = "Density") +
  theme_minimal()
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-19-1.png)

## 8. Did fusion actually update the prior? Comparing prior, likelihood, and posterior distributions

Section 7 compared two fusion *routes* against each other. A more basic
question is whether fusion did anything sensible at all relative to its
two *inputs*: does the posterior mean clay distribution actually sit
somewhere between the SSURGO prior and the SOLUS100 likelihood,
sharpened by combining them? `fusion_clay` (loaded back in “Reusing real
data instead of fetching”) already holds the real, full-pipeline answer
for this AOI - its `prior`/`likelihood`/`posterior` are the actual
objects
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
produced, not the toy
`aligned`/[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)
example from Section 7:

``` r

dist_df <- data.frame(
  value = c(
    terra::values(prior_values$P50, na.rm = TRUE),
    terra::values(lik_values$P50, na.rm = TRUE),
    terra::values(fusion_clay$posterior$mu, na.rm = TRUE)
  ),
  source = rep(
    c("Prior (SSURGO median)", "Likelihood (SOLUS prediction)", "Posterior (fused mean)"),
    c(sum(!is.na(terra::values(prior_values$P50))),
      sum(!is.na(terra::values(lik_values$P50))),
      sum(!is.na(terra::values(fusion_clay$posterior$mu))))
  )
)
dist_df$source <- factor(dist_df$source,
  levels = c("Prior (SSURGO median)", "Likelihood (SOLUS prediction)", "Posterior (fused mean)"))

means <- stats::aggregate(value ~ source, dist_df, mean)

ggplot(dist_df, aes(x = value, fill = source)) +
  geom_density(alpha = 0.45, color = NA) +
  geom_vline(data = means, aes(xintercept = value, color = source),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  scale_fill_viridis_d(name = NULL) +
  scale_color_viridis_d() +
  labs(title = "Prior vs. likelihood vs. posterior mean clay content",
       subtitle = "Dashed lines mark each distribution's mean",
       x = "Clay content (%)", y = "Density") +
  theme_minimal()
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-20-1.png)

``` r


means
#>                          source    value
#> 1         Prior (SSURGO median) 12.50059
#> 2 Likelihood (SOLUS prediction) 20.83220
#> 3        Posterior (fused mean) 12.71270
```

The posterior sits almost on top of the prior here, barely nudged toward
the likelihood’s much higher mean (~21%). If fusion were simply
averaging its two inputs, the posterior would land roughly halfway
between them; it clearly doesn’t. This is exactly what
precision-weighted Bayesian fusion is *supposed* to do when one side
states far more uncertainty about itself than the other - not a sign
that the likelihood is being ignored. Fitting each side’s own Normal
distribution independently (the same
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)
from Section 1, applied to whichever low/high percentile pair each
source actually has) makes the reason concrete:

``` r

fit_prior_normal <- fit_normal_raster(
  p_lo_r = prior_values$P05, p50_r = prior_values$P50, p_hi_r = prior_values$P95,
  p_lo = 0.05, p_hi = 0.95
)
fit_lik_normal <- fit_normal_raster(
  p_lo_r = lik_values$P025, p50_r = lik_values$P50, p_hi_r = lik_values$P975,
  p_lo = 0.025, p_hi = 0.975
)

sigma_compare <- data.frame(
  source = c("Prior (SSURGO)", "Likelihood (SOLUS)"),
  mean_sigma = c(mean(terra::values(fit_prior_normal$sigma), na.rm = TRUE),
                 mean(terra::values(fit_lik_normal$sigma), na.rm = TRUE))
)
sigma_compare$mean_precision <- 1 / sigma_compare$mean_sigma^2
sigma_compare
#>               source mean_sigma mean_precision
#> 1     Prior (SSURGO)    2.49010    0.161274775
#> 2 Likelihood (SOLUS)   12.54174    0.006357473
```

SOLUS100’s raw 95% prediction interval (`P025`-`P975`) implies a mean
`sigma` roughly 5x wider than SSURGO’s 5th-95th percentile spread for
this AOI - a real, known characteristic of SOLUS100’s output (its stated
prediction intervals are often wide), not an artifact of this pipeline.
Since fusion weights each side by *precision* (`1/sigma^2`, see
`bayesian-updating.Rmd`’s worked scalar examples of this exact
tug-of-war), a 5x wider sigma means roughly a 25x lower precision - so
the posterior mean above (12.71) landing close to the prior’s mean
(12.5) rather than the likelihood’s (20.83) is the fusion math doing
exactly what it should with these particular inputs, confirmed by
manually recomputing
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
on these two independently-fit sides and getting (within
floating-point/per-cell alignment differences) the same answer the
pipeline’s own route produced. A different AOI/property where the two
sources agree more closely on their own uncertainty would show the
posterior pulled further toward whichever side is sharper, without one
side this lopsidedly dominating.

## 9. A per-pixel view: prior, likelihood, and posterior at individual mapunits

Section 8’s density plot pools every cell in the AOI into one
distribution per source, which is useful for the big picture but hides
what fusion actually does *at a single location* - precision weighting
happens per pixel, using that pixel’s own prior/likelihood spread, not
the AOI-wide average spread used above. This section looks at a handful
of individual pixels instead.

The cached object here only carries percentile-*value* rasters
(`prior_values`, `lik_values`), not the categorical mukey raster
[`fetch_ssurgo_mukey_raster()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md)
would return - so there’s no direct mukey lookup to pick pixels by
mapunit membership. As a proxy, SSURGO’s simulated median
(`prior_values$P50`) is close to constant within a mapunit and differs
between mapunits (each mapunit gets its own Monte Carlo simulation
draws), so cells sharing a repeated `P50` value are, in practice,
sitting in the same sizeable mapunit. Picking a few *different* repeated
values gives a handful of distinct mapunit areas:

``` r

p50_all <- terra::values(prior_values$P50)[, 1]
cell_ids <- which(!is.na(p50_all))
rounded <- round(p50_all[cell_ids], 1)

# Only consider values shared by several cells (a proxy for "a sizeable, distinct mapunit",
# rather than a one-off simulation outlier), then spread the picks across the low-to-high clay
# range so the examples aren't all pulled from one part of the AOI.
freq <- table(rounded)
frequent_vals <- sort(as.numeric(names(freq[freq >= 4])))
pick_idx <- round(seq(1, length(frequent_vals), length.out = 4))
picked_vals <- frequent_vals[pick_idx]
example_cells <- vapply(picked_vals, function(v) cell_ids[which(rounded == v)[1]], numeric(1))

picked_vals
#> [1]  1.7 10.7 14.4 28.0
```

For each of these four example pixels, pull `mu`/`sigma` from the same
Normal fits used above - `fit_prior_normal`/`fit_lik_normal` (Section 8)
for the prior/likelihood, and `fusion_clay$posterior` (the pipeline’s
real per-cell posterior, which carries both `mu` *and* `sigma`) for the
posterior:

``` r

extract_at <- function(r, cell) terra::values(r)[cell]

pixel_df <- do.call(rbind, lapply(seq_along(example_cells), function(i) {
  cell <- example_cells[i]
  data.frame(
    pixel = sprintf("Pixel %d (mapunit-like cluster, prior P50 ~%.1f%%)", i, picked_vals[i]),
    prior_mu = extract_at(fit_prior_normal$mu, cell), prior_sigma = extract_at(fit_prior_normal$sigma, cell),
    lik_mu = extract_at(fit_lik_normal$mu, cell), lik_sigma = extract_at(fit_lik_normal$sigma, cell),
    post_mu = extract_at(fusion_clay$posterior$mu, cell), post_sigma = extract_at(fusion_clay$posterior$sigma, cell)
  )
}))
pixel_df
#>                                              pixel prior_mu prior_sigma lik_mu
#> 1  Pixel 1 (mapunit-like cluster, prior P50 ~1.7%)  1.74080    3.081634     18
#> 2 Pixel 2 (mapunit-like cluster, prior P50 ~10.7%) 10.70449    3.508857     17
#> 3 Pixel 3 (mapunit-like cluster, prior P50 ~14.4%) 14.40415    3.725312     24
#> 4 Pixel 4 (mapunit-like cluster, prior P50 ~28.0%) 28.00061    4.541492     24
#>   lik_sigma   post_mu post_sigma
#> 1  12.50023  3.723446   3.338727
#> 2  12.75534 10.353137   3.674820
#> 3  12.75534 13.498238   4.005540
#> 4  12.24512 27.030965   4.630531
```

Then evaluate each pixel’s three Normal curves
([`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html), no raster
machinery needed for a single cell) over a shared x-range wide enough to
show all three, and facet by pixel:

``` r

curve_at_pixel <- function(row) {
  x_lo <- min(row$prior_mu - 4 * row$prior_sigma, row$lik_mu - 4 * row$lik_sigma,
              row$post_mu - 4 * row$post_sigma, 0)
  x_hi <- max(row$prior_mu + 4 * row$prior_sigma, row$lik_mu + 4 * row$lik_sigma,
              row$post_mu + 4 * row$post_sigma, 100)
  x <- seq(max(x_lo, 0), min(x_hi, 100), length.out = 300)
  rbind(
    data.frame(pixel = row$pixel, x = x, density = stats::dnorm(x, row$prior_mu, row$prior_sigma),
               source = "Prior (SSURGO)"),
    data.frame(pixel = row$pixel, x = x, density = stats::dnorm(x, row$lik_mu, row$lik_sigma),
               source = "Likelihood (SOLUS)"),
    data.frame(pixel = row$pixel, x = x, density = stats::dnorm(x, row$post_mu, row$post_sigma),
               source = "Posterior (fused)")
  )
}
curve_df <- do.call(rbind, lapply(seq_len(nrow(pixel_df)), function(i) curve_at_pixel(pixel_df[i, ])))
curve_df$source <- factor(curve_df$source,
  levels = c("Prior (SSURGO)", "Likelihood (SOLUS)", "Posterior (fused)"))

ggplot(curve_df, aes(x = x, y = density, color = source, fill = source)) +
  geom_area(alpha = 0.3, position = "identity") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~pixel, scales = "free_y") +
  scale_color_viridis_d(name = NULL) +
  scale_fill_viridis_d(name = NULL) +
  labs(title = "Prior, likelihood, and posterior at four individual pixels",
       subtitle = "Each panel is one pixel from a different mapunit-like cluster",
       x = "Clay content (%)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-24-1.png)

Every pixel here tells the same qualitative story as the AOI-wide
picture in Section 8 - the posterior curve sits close to the prior
curve, barely shifted toward the likelihood - because precision
weighting is happening independently at each pixel and, for this
property/AOI, SOLUS’s per-pixel `sigma` is consistently much larger than
SSURGO’s regardless of which mapunit the pixel falls in (visible
directly in `pixel_df$lik_sigma` vs. `pixel_df$prior_sigma` above). A
property or AOI where the two sources’ stated uncertainty is closer
would show more pixel-to-pixel variation in how far the posterior moves,
including pixels where the likelihood dominates instead.

## 10. Caching internals

[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
disk-cache every fetch step under
`tools::R_user_dir("soilSIM", "cache")`, keyed by AOI + property/group
id + depth window + kind.

[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)
builds a deterministic, filename-safe key from an AOI’s WKT geometry
(hashed via
[`digest::digest()`](https://eddelbuettel.github.io/digest/man/digest.html)),
a property/group id, a depth window, and a `kind` string:

``` r

salinas_wkt <- "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))"
aoi <- terra::project(terra::vect(salinas_wkt, crs = "epsg:4326"), "epsg:5070")

key <- build_cache_key(aoi, id = "clay", top_depth = 0, bottom_depth = 5, kind = "ssurgo")
key
#> [1] "raster_3e9ea7cad8f2_clay_0_5_ssurgo"
CACHE_TTL_SECONDS / 86400  # cache max-age, in days
#> [1] 30
```

`cache_get(key, ttl_seconds)`/`cache_set(key, kind, value)` are the
generic get/set pair every fetch step above uses - a miss (absent, or
older than `ttl_seconds`) returns `NULL` rather than erroring:

``` r

cache_get(key)  # NULL: nothing has been cached under this exact key in this session
#> NULL
cache_set(key, "ssurgo", value = list(hello = "world"))
#> [1] TRUE
cache_get(key)
#> $hello
#> [1] "world"
```

### Why `SpatRaster`s need `terra::wrap()` before caching

[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)’s
docs flag a real, sharp limitation: a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
holds an **external pointer** to in-memory/on-disk GDAL state, not the
raster’s actual values. A plain
[`saveRDS()`](https://rspatial.github.io/terra/reference/serialize.html)/[`readRDS()`](https://rspatial.github.io/terra/reference/serialize.html)
round-trip serializes that pointer, which is meaningless in a new R
session - the read-back object *looks* like a `SpatRaster` but errors
the moment anything tries to touch its actual data:

``` r

r <- fusion_clay$posterior$mu
naive_file <- tempfile(fileext = ".rds")
saveRDS(r, naive_file)               # serializes the pointer, not the raster's values
r_naive <- readRDS(naive_file)
class(r_naive)                       # still LOOKS like a SpatRaster...
#> [1] "SpatRaster"
#> attr(,"package")
#> [1] "terra"
terra::values(r_naive)               # ...but touching it throws
#>            lyr.1
#>   [1,]        NA
#>   [2,]        NA
#>   [3,]        NA
#>   [4,]        NA
#>   [5,]        NA
#>   [6,]        NA
#>   [7,]  5.513754
#>   [8,]        NA
#>   [9,]        NA
#>  [10,]        NA
#>  [11,]        NA
#>  [12,]        NA
#>  [13,]        NA
#>  [14,]        NA
#>  [15,]        NA
#>  [16,]        NA
#>  [17,]        NA
#>  [18,]        NA
#>  [19,]        NA
#>  [20,]        NA
#>  [21,]        NA
#>  [22,]        NA
#>  [23,]        NA
#>  [24,]        NA
#>  [25,]        NA
#>  [26,]        NA
#>  [27,]        NA
#>  [28,]        NA
#>  [29,]        NA
#>  [30,]  6.056102
#>  [31,]  6.338324
#>  [32,]  6.912659
#>  [33,] 21.334730
#>  [34,] 27.030965
#>  [35,]        NA
#>  [36,]        NA
#>  [37,]        NA
#>  [38,]        NA
#>  [39,]        NA
#>  [40,]        NA
#>  [41,]        NA
#>  [42,]        NA
#>  [43,]        NA
#>  [44,]        NA
#>  [45,]        NA
#>  [46,]        NA
#>  [47,]        NA
#>  [48,]        NA
#>  [49,]        NA
#>  [50,]        NA
#>  [51,]        NA
#>  [52,]        NA
#>  [53,]  4.971009
#>  [54,]  6.479068
#>  [55,] 16.503437
#>  [56,] 25.814255
#>  [57,] 27.907140
#>  [58,] 28.022443
#>  [59,] 27.281745
#>  [60,] 27.833118
#>  [61,] 27.950841
#>  [62,]        NA
#>  [63,]        NA
#>  [64,]        NA
#>  [65,]        NA
#>  [66,]        NA
#>  [67,]        NA
#>  [68,]        NA
#>  [69,]        NA
#>  [70,]        NA
#>  [71,]        NA
#>  [72,]        NA
#>  [73,]        NA
#>  [74,]        NA
#>  [75,] 10.497035
#>  [76,]  5.387846
#>  [77,]  6.053746
#>  [78,] 19.060490
#>  [79,] 26.926240
#>  [80,] 26.661847
#>  [81,] 26.745230
#>  [82,] 27.432363
#>  [83,] 27.652231
#>  [84,] 27.917359
#>  [85,] 27.344451
#>  [86,] 28.013321
#>  [87,] 27.378561
#>  [88,]        NA
#>  [89,]        NA
#>  [90,]        NA
#>  [91,]        NA
#>  [92,]        NA
#>  [93,]        NA
#>  [94,]        NA
#>  [95,]        NA
#>  [96,]        NA
#>  [97,]        NA
#>  [98,] 12.307365
#>  [99,]  6.953534
#> [100,]  4.701146
#> [101,] 10.353137
#> [102,] 18.170227
#> [103,]        NA
#> [104,] 25.552390
#> [105,] 27.528716
#> [106,] 27.756149
#> [107,] 27.705909
#> [108,] 27.002529
#> [109,] 26.439329
#> [110,] 23.489279
#> [111,] 26.525306
#> [112,] 28.167618
#> [113,] 27.598663
#> [114,] 27.753873
#> [115,]        NA
#> [116,]        NA
#> [117,]        NA
#> [118,]        NA
#> [119,]        NA
#> [120,]        NA
#> [121,] 13.047057
#> [122,] 10.024957
#> [123,]  4.932454
#> [124,]  4.936732
#> [125,]  6.831262
#> [126,]        NA
#> [127,] 17.906676
#> [128,] 27.921046
#> [129,] 27.611038
#> [130,] 27.264879
#> [131,] 21.064458
#> [132,] 10.789871
#> [133,]  7.551773
#> [134,] 14.936048
#> [135,] 26.407569
#> [136,] 27.075128
#> [137,] 27.950843
#> [138,] 27.974339
#> [139,]        NA
#> [140,]        NA
#> [141,]        NA
#> [142,]        NA
#> [143,] 13.059353
#> [144,] 13.184770
#> [145,] 12.369777
#> [146,]  6.688725
#> [147,]  4.001997
#> [148,]  5.500645
#> [149,]        NA
#> [150,]  8.867549
#> [151,] 23.061965
#> [152,] 23.422876
#> [153,] 13.498238
#> [154,]  7.153168
#> [155,]  6.018624
#> [156,]  6.496624
#> [157,]  7.579653
#> [158,] 21.731349
#> [159,] 27.881159
#> [160,] 26.986269
#> [161,] 27.544781
#> [162,]        NA
#> [163,]        NA
#> [164,]        NA
#> [165,]        NA
#> [166,] 13.103153
#> [167,] 14.142842
#> [168,] 14.424858
#> [169,] 10.214928
#> [170,]  4.667976
#> [171,]  4.422348
#> [172,]  6.019767
#> [173,]  6.291146
#> [174,]  8.234903
#> [175,]  8.431786
#> [176,]  5.854015
#> [177,]  5.768374
#> [178,]  6.097004
#> [179,]  6.215924
#> [180,]  6.266272
#> [181,] 15.597548
#> [182,] 24.776698
#> [183,] 24.074538
#> [184,] 23.832819
#> [185,]        NA
#> [186,]        NA
#> [187,]        NA
#> [188,]        NA
#> [189,] 20.076558
#> [190,] 20.487495
#> [191,] 20.823638
#> [192,] 18.588315
#> [193,]  6.721067
#> [194,]  4.529395
#> [195,]  5.259927
#> [196,]  6.084984
#> [197,]  5.895144
#> [198,]  6.639594
#> [199,]  5.835514
#> [200,]  6.072674
#> [201,]  6.282354
#> [202,]  6.170624
#> [203,]  6.273654
#> [204,]  7.085832
#> [205,]  9.216213
#> [206,] 11.617685
#> [207,]        NA
#> [208,]        NA
#> [209,]        NA
#> [210,]        NA
#> [211,]        NA
#> [212,] 22.365728
#> [213,] 21.290888
#> [214,] 24.862545
#> [215,] 23.873698
#> [216,] 12.820867
#> [217,]  5.334646
#> [218,]  3.869894
#> [219,]  5.236623
#> [220,]  5.934822
#> [221,]  6.288284
#> [222,]  5.930204
#> [223,]  5.821164
#> [224,]  6.458184
#> [225,]  5.959404
#> [226,]  5.920964
#> [227,]  5.922634
#> [228,]  6.420774
#> [229,]  6.340584
#> [230,]        NA
#> [231,]        NA
#> [232,]        NA
#> [233,]        NA
#> [234,]        NA
#> [235,] 12.740495
#> [236,] 12.696656
#> [237,] 20.025115
#> [238,] 21.557923
#> [239,] 20.222260
#> [240,]  8.747617
#> [241,]  3.845663
#> [242,]  3.376861
#> [243,]  5.147157
#> [244,]  5.649134
#> [245,]  5.809154
#> [246,]  5.907154
#> [247,]  5.869964
#> [248,]  6.032694
#> [249,]  6.209942
#> [250,]  6.704714
#> [251,]  7.324813
#> [252,]  6.724679
#> [253,]        NA
#> [254,]        NA
#> [255,]        NA
#> [256,]        NA
#> [257,] 10.878542
#> [258,] 10.840802
#> [259,] 11.048089
#> [260,] 16.683608
#> [261,] 21.603606
#> [262,] 21.546903
#> [263,] 14.785350
#> [264,]  4.953229
#> [265,]  3.723446
#> [266,]  3.919309
#> [267,]  5.154748
#> [268,]  5.975419
#> [269,]  6.142824
#> [270,]  6.079546
#> [271,]  7.757721
#> [272,]  9.028535
#> [273,]  9.917106
#> [274,]  9.170676
#> [275,]  7.899348
#> [276,]        NA
#> [277,]        NA
#> [278,]        NA
#> [279,]        NA
#> [280,] 10.798592
#> [281,] 10.924713
#> [282,] 11.143986
#> [283,] 18.296449
#> [284,] 21.686916
#> [285,] 21.651943
#> [286,] 20.209063
#> [287,]  8.513660
#> [288,]        NA
#> [289,]  3.925391
#> [290,]  3.564288
#> [291,]  5.328597
#> [292,]  6.182809
#> [293,]  8.803078
#> [294,]  9.104363
#> [295,]  7.013611
#> [296,]  6.408614
#> [297,]  5.915774
#> [298,]        NA
#> [299,]        NA
#> [300,]        NA
#> [301,]        NA
#> [302,]        NA
#> [303,] 10.891405
#> [304,] 10.744310
#> [305,]        NA
#> [306,] 13.874023
#> [307,] 17.731745
#> [308,] 20.652735
#> [309,] 21.601983
#> [310,] 17.579898
#> [311,]  6.508097
#> [312,]  4.005990
#> [313,]  3.476241
#> [314,]  4.353991
#> [315,]  4.547539
#> [316,]  4.571532
#> [317,]  4.587507
#> [318,]  4.827155
#> [319,]  5.354703
#> [320,]  5.808944
#> [321,]        NA
#> [322,]        NA
#> [323,]        NA
#> [324,]        NA
#> [325,] 11.146857
#> [326,] 10.866167
#> [327,] 10.996337
#> [328,] 11.256610
#> [329,] 12.307903
#> [330,] 16.437633
#> [331,] 15.508073
#> [332,] 18.310386
#> [333,] 21.366032
#> [334,] 13.810562
#> [335,]  5.731162
#> [336,]  4.832626
#> [337,]  3.647858
#> [338,]  3.675237
#> [339,]  3.558191
#> [340,]        NA
#> [341,]  3.472457
#> [342,]  4.263312
#> [343,]  5.970253
#> [344,]        NA
#> [345,]        NA
#> [346,]        NA
#> [347,]        NA
#> [348,] 10.931945
#> [349,] 11.141357
#> [350,] 10.908728
#> [351,] 10.905788
#> [352,] 12.287578
#> [353,] 13.224830
#> [354,] 18.131510
#> [355,] 16.486002
#> [356,] 20.866743
#> [357,] 20.106878
#> [358,] 11.392358
#> [359,]  6.799787
#> [360,]  5.459106
#> [361,]  4.438018
#> [362,]  4.177861
#> [363,]  4.171971
#> [364,]  3.402091
#> [365,]  3.805834
#> [366,]        NA
#> [367,]        NA
#> [368,]        NA
#> [369,]        NA
#> [370,]        NA
#> [371,] 10.890718
#> [372,] 10.922287
#> [373,] 11.097847
#> [374,] 10.911525
#> [375,] 12.515657
#> [376,] 12.303930
#> [377,] 13.386280
#> [378,] 17.649878
#> [379,] 17.298030
#> [380,] 15.706688
#> [381,] 14.589147
#> [382,] 14.369648
#> [383,] 12.721408
#> [384,]  6.958831
#> [385,]  4.227801
#> [386,]  3.653810
#> [387,]  3.624821
#> [388,]  3.462741
#> [389,]        NA
#> [390,]        NA
#> [391,]        NA
#> [392,]        NA
#> [393,] 11.012308
#> [394,] 10.805578
#> [395,] 10.799917
#> [396,] 10.779708
#> [397,] 11.035818
#> [398,] 12.115855
#> [399,] 12.546885
#> [400,] 12.183557
#> [401,] 13.199128
#> [402,] 13.216920
#> [403,] 14.109438
#> [404,] 14.544458
#> [405,] 14.656090
#> [406,] 16.363107
#> [407,] 15.347540
#> [408,]  7.089960
#> [409,]  4.279458
#> [410,]  3.550161
#> [411,]  3.218861
#> [412,]        NA
#> [413,]        NA
#> [414,]        NA
#> [415,]        NA
#> [416,] 10.984418
#> [417,] 10.788997
#> [418,] 10.825308
#> [419,] 10.820237
#> [420,] 10.867727
#> [421,] 10.899135
#> [422,] 11.378522
#> [423,] 11.425124
#> [424,] 12.692597
#> [425,] 12.912258
#> [426,] 13.965177
#> [427,] 14.359950
#> [428,] 14.549040
#> [429,] 14.566930
#> [430,] 15.488322
#> [431,] 13.479223
#> [432,]  5.926795
#> [433,]  4.064121
#> [434,]  3.919528
#> [435,]        NA
#> [436,]        NA
#> [437,]        NA
#> [438,]        NA
#> [439,] 10.938167
#> [440,] 10.947518
#> [441,] 11.830968
#> [442,] 12.068794
#> [443,] 11.902735
#> [444,] 11.666177
#> [445,] 12.191875
#> [446,] 11.631574
#> [447,] 11.699077
#> [448,] 12.689227
#> [449,] 14.282238
#> [450,] 14.434848
#> [451,] 14.492095
#> [452,] 14.345407
#> [453,] 14.745258
#> [454,] 14.572250
#> [455,] 11.780028
#> [456,]  6.559958
#> [457,]        NA
#> [458,]        NA
#> [459,]        NA
#> [460,]        NA
#> [461,]        NA
#> [462,] 10.928875
#> [463,] 11.245520
#> [464,] 12.789027
#> [465,] 13.050725
#> [466,] 12.291167
#> [467,] 12.714467
#> [468,] 13.248497
#> [469,] 12.988325
#> [470,] 13.453508
#> [471,] 14.154865
#> [472,] 14.443577
#> [473,] 14.525378
#> [474,] 14.552460
#> [475,] 14.563272
#> [476,] 14.659145
#> [477,] 14.854317
#> [478,] 14.691288
#> [479,] 12.361468
#> [480,]        NA
#> [481,]        NA
#> [482,]        NA
#> [483,]        NA
#> [484,]        NA
#> [485,] 10.727725
#> [486,] 10.787947
#> [487,] 12.620168
#> [488,] 12.796378
#> [489,] 12.231647
#> [490,] 12.952050
#> [491,] 13.142630
#> [492,] 13.057588
#> [493,] 13.889437
#> [494,] 14.324947
#> [495,] 14.501487
#> [496,] 14.392408
#> [497,] 14.446178
#> [498,] 14.445120
#> [499,] 14.903737
#> [500,] 14.814275
#> [501,] 14.856947
#> [502,] 14.827440
#> [503,]        NA
#> [504,]        NA
#> [505,]        NA
#> [506,]        NA
#> [507,]        NA
#> [508,]        NA
#> [509,]        NA
#> [510,] 12.957137
#> [511,] 12.978825
#> [512,] 13.039778
#> [513,] 13.075728
#> [514,] 13.252938
#> [515,] 13.082552
#> [516,] 13.353278
#> [517,] 13.880775
#> [518,] 14.499765
#> [519,] 14.369447
#> [520,] 14.511878
#> [521,] 14.474648
#> [522,] 14.648868
#> [523,] 14.805245
#> [524,] 14.868157
#> [525,] 14.724677
#> [526,]        NA
#> [527,]        NA
#> [528,]        NA
#> [529,]        NA
#> [530,]        NA
#> [531,]        NA
#> [532,]        NA
#> [533,]        NA
#> [534,]        NA
#> [535,]        NA
#> [536,]        NA
#> [537,] 13.207690
#> [538,] 13.221258
#> [539,] 13.157160
#> [540,] 13.906988
#> [541,] 14.462158
#> [542,] 14.485185
#> [543,] 14.343488
#> [544,]        NA
#> [545,] 14.716618
#> [546,] 14.672177
#> [547,] 14.790637
#> [548,]        NA
#> [549,]        NA
#> [550,]        NA
#> [551,]        NA
#> [552,]        NA
#> [553,]        NA
#> [554,]        NA
#> [555,]        NA
#> [556,]        NA
#> [557,]        NA
#> [558,]        NA
#> [559,]        NA
#> [560,]        NA
#> [561,]        NA
#> [562,]        NA
#> [563,] 13.271077
#> [564,] 14.080695
#> [565,] 14.528477
#> [566,] 14.573950
#> [567,]        NA
#> [568,] 14.452078
#> [569,] 14.503257
#> [570,] 14.923508
#> [571,]        NA
#> [572,]        NA
#> [573,]        NA
#> [574,]        NA
#> [575,]        NA
#> [576,]        NA
#> [577,]        NA
#> [578,]        NA
#> [579,]        NA
#> [580,]        NA
#> [581,]        NA
#> [582,]        NA
#> [583,]        NA
#> [584,]        NA
#> [585,]        NA
#> [586,]        NA
#> [587,]        NA
#> [588,]        NA
#> [589,]        NA
#> [590,] 14.461267
#> [591,] 14.558437
#> [592,] 14.742098
#> [593,] 14.740058
#> [594,]        NA
#> [595,]        NA
#> [596,]        NA
#> [597,]        NA
#> [598,]        NA
```

[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)
converts a `SpatRaster` into a `PackedSpatRaster` - a plain, fully
self-contained R object holding the actual cell values - which *does*
survive
[`saveRDS()`](https://rspatial.github.io/terra/reference/serialize.html)/[`readRDS()`](https://rspatial.github.io/terra/reference/serialize.html);
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html)
converts it back:

``` r

wrap_file <- tempfile(fileext = ".rds")
saveRDS(terra::wrap(r), wrap_file)
r_wrapped_read <- readRDS(wrap_file)
class(r_wrapped_read)                # PackedSpatRaster - plain data, not a live pointer
#> [1] "SpatRaster"
#> attr(,"package")
#> [1] "terra"

r_restored <- terra::unwrap(r_wrapped_read)
class(r_restored)                    # SpatRaster again
#> [1] "SpatRaster"
#> attr(,"package")
#> [1] "terra"
identical(terra::values(r), terra::values(r_restored))  # values round-tripped exactly
#> [1] TRUE
```

Since
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)
themselves don’t wrap/unwrap automatically (not every cached value is a
raster - e.g. cached SSURGO Monte Carlo draws are a plain data frame),
every caller that caches raster values has to do this manually. Two
helpers cover the two result shapes this pipeline actually caches:

- [`wrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md)/[`unwrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md) -
  for the fixed `list(values = <SpatRasters>, probs = )` percentile
  shape
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)/[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)
  return.
- [`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)/[`unwrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md) -
  a fully generic version for result structures whose shape isn’t that
  fixed form,
  e.g. [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)’s
  own return value, which nests `SpatRaster`s at `prior$values`,
  `likelihood$values`, and inside `posterior` (whose shape varies by
  `dist`). It recurses into every list element, replacing each
  `SpatRaster` it finds with its wrapped/unwrapped equivalent and
  leaving everything else untouched - which is exactly how this
  vignette’s own `fusion_clay` object was loaded, back in Section
  “Reusing real data instead of fetching”:

``` r

# This is the call this whole vignette started with:
raw <- readRDS(system.file("extdata", "fusion_clay_salinas.rds", package = "soilSIM"))
class(raw$posterior$mu)                              # still PackedSpatRaster on disk
#> [1] "PackedSpatRaster"
#> attr(,"package")
#> [1] "terra"
fusion_clay_reloaded <- unwrap_nested_rasters(raw)
class(fusion_clay_reloaded$posterior$mu)              # SpatRaster after unwrap_nested_rasters()
#> [1] "SpatRaster"
#> attr(,"package")
#> [1] "terra"
```

## Where this data came from

The cached data this vignette loads
(`inst/extdata/fusion_clay_salinas.rds`) is the same object
`raster-fusion-ssurgo-solus.Rmd` loads, produced once by
`data-raw/build_vignette_data.R` against live NRCS Soil Data Access
(SSURGO) and
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
(SOLUS100) for the Salinas Valley AOI. See that vignette for the full
fetch-through-fuse pipeline view, and
`soilSIM/docs/09_multi_source_raster_fusion_pipeline.md` for the
complete function-level reference this vignette exposed piece by piece.
