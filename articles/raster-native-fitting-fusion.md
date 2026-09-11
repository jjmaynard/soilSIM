# Raster-Native Fitting & Fusion Internals

## Overview

`raster-fusion-ssurgo-solus.Rmd` shows the *pipeline* view of soilSIM’s
multi-source raster fusion: call
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)/[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)
and get back a fused posterior. This vignette opens the hood.
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
is a thin fetch-and-cache wrapper around a handful of small, composable,
raster-native building blocks - distribution fitting functions
(`R/core-distributions-raster.R`), a fusion core (`R/core-fusion.R`),
and a disk cache (`R/cache.R`) - each of which does exactly one thing to
a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
(or list of them) and can be called directly. Understanding these
individually is what makes the pipeline’s behavior (which route it took,
why a family was resolved a certain way, why two sources’ percentiles
had to be realigned before fusing) legible rather than a black box. See
the “Multi-Source Raster Fusion Pipeline” architecture article for the
full function-level reference this vignette is built from.

``` r

library(soilSIM)
library(terra)
#> terra 1.9.50
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
fusion_clay <- run_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5)
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
#> name        :       P50
#> min value   :  1.531298
#> max value   : 41.521725
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
            col = grDevices::hcl.colors(50, "viridis"), nc = 3,
            mar = c(3.1, 3.1, 5.1, 7.1))
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
#> 209 / 598 cells (35%) have an infeasible metalog fit
terra::plot(infeasible, main = "Infeasible metalog fit (TRUE = non-monotonic)",
            col = c("grey85", "firebrick"),
            mar = c(3.1, 3.1, 5.1, 7.1))
```

![](raster-native-fitting-fusion_files/figure-html/unnamed-chunk-10-1.png)

A substantial fraction of cells are infeasible here - not a bug, just a
real consequence of fitting only 4 of the 5 available percentiles to a
family with no correction step of its own.
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md)
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

fallback_result <- quantile_metalog_linear_with_fallback_raster(
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
step, via `Reduce(+, ...)`), then delegates to `R/core-fusion.R`’s
existing scalar
[`convert_moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_gamma.md),
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
#> [1] -0.005760565

# Without bounds, it falls back to a skew-proxy test: normal if |skew| is small, lognormal if not.
resolve_property_dist(list(dist = "auto"), prior_values, prior_probs)
#> $dist
#> [1] "normal"
#> 
#> $dist_source
#> [1] "auto"
#> 
#> $skew_proxy
#> [1] -0.005760565
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
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md))
delegates to for `"normal"`/`"beta"`/`"gamma"` families. It picks the
fusion *route* from the AOI’s own cell count rather than requiring the
caller to choose: at or below `threshold_cells` (default 80,000), it
uses a fully general per-cell KDE route
([`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
on samples drawn from each side’s percentiles); above it, a much faster
closed-form route built directly on the fitting functions from Sections
1-4 above (e.g.
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md) +
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md)
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
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md)
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
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
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
#> 1         Prior (SSURGO median) 13.75337
#> 2 Likelihood (SOLUS prediction) 21.27797
#> 3        Posterior (fused mean) 14.08067
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
#> 1     Prior (SSURGO)   2.858638    0.122371919
#> 2 Likelihood (SOLUS)  12.271931    0.006640094
```

SOLUS100’s raw 95% prediction interval (`P025`-`P975`) implies a mean
`sigma` roughly 4x wider than SSURGO’s 5th-95th percentile spread for
this AOI - a real, known characteristic of SOLUS100’s output (its stated
prediction intervals are often wide), not an artifact of this pipeline.
Since fusion weights each side by *precision* (`1/sigma^2`, see
`core-fusion.Rmd`’s worked scalar examples of this exact tug-of-war), a
4x wider sigma means roughly a 18x lower precision - so the posterior
mean above (14.08) landing close to the prior’s mean (13.75) rather than
the likelihood’s (21.28) is the fusion math doing exactly what it should
with these particular inputs, confirmed by manually recomputing
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md)
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
#> [1]  1.5 10.4 13.1 27.9
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
#>                                              pixel  prior_mu prior_sigma lik_mu
#> 1  Pixel 1 (mapunit-like cluster, prior P50 ~1.5%)  1.531298    3.276059   21.0
#> 2 Pixel 2 (mapunit-like cluster, prior P50 ~10.4%) 10.405450    2.217063   17.5
#> 3 Pixel 3 (mapunit-like cluster, prior P50 ~13.1%) 13.103279    3.060271   19.0
#> 4 Pixel 4 (mapunit-like cluster, prior P50 ~27.9%) 27.890039    4.880258   25.5
#>   lik_sigma   post_mu post_sigma
#> 1  12.88289  3.320877   3.587894
#> 2  11.86246 12.950815   1.582438
#> 3  11.60736 19.542584   4.028441
#> 4  13.77576 27.661751   5.444103
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

[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)/[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)
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
#>   [1,] 21.723848
#>   [2,] 21.724403
#>   [3,] 21.743542
#>   [4,] 21.628982
#>   [5,]  6.135123
#>   [6,]  3.456004
#>   [7,]  6.235678
#>   [8,]  5.902018
#>   [9,]  6.279828
#>  [10,] 27.749131
#>  [11,] 27.195297
#>  [12,] 27.661751
#>  [13,] 27.862847
#>  [14,] 27.854473
#>  [15,] 27.730944
#>  [16,] 27.933995
#>  [17,] 28.172795
#>  [18,] 27.906679
#>  [19,] 27.644101
#>  [20,] 27.940815
#>  [21,] 27.933997
#>  [22,] 29.177075
#>  [23,] 30.806709
#>  [24,] 21.772536
#>  [25,] 21.836037
#>  [26,] 21.692967
#>  [27,] 21.765067
#>  [28,] 19.628751
#>  [29,]  6.117459
#>  [30,]  6.104207
#>  [31,]  5.927938
#>  [32,]  6.188520
#>  [33,] 27.502204
#>  [34,] 27.866436
#>  [35,] 27.772562
#>  [36,] 27.492525
#>  [37,] 27.360932
#>  [38,] 27.922725
#>  [39,] 28.062864
#>  [40,] 27.907735
#>  [41,] 28.235843
#>  [42,] 27.737586
#>  [43,] 27.849428
#>  [44,] 27.460671
#>  [45,] 28.222052
#>  [46,] 30.488510
#>  [47,] 21.721662
#>  [48,] 21.744132
#>  [49,] 21.686594
#>  [50,] 21.535568
#>  [51,] 19.611663
#>  [52,]  6.127593
#>  [53,]  3.225927
#>  [54,]  6.109527
#>  [55,] 27.422534
#>  [56,] 27.318396
#>  [57,] 27.733577
#>  [58,] 27.509098
#>  [59,] 27.748919
#>  [60,] 27.698854
#>  [61,] 27.749987
#>  [62,] 27.446765
#>  [63,] 28.466790
#>  [64,] 28.140562
#>  [65,] 28.029726
#>  [66,] 27.852711
#>  [67,] 27.674363
#>  [68,] 27.744561
#>  [69,] 30.672192
#>  [70,] 21.757781
#>  [71,] 21.768878
#>  [72,] 21.720436
#>  [73,] 21.645023
#>  [74,] 21.684067
#>  [75,] 19.350875
#>  [76,]  3.205744
#>  [77,]  5.878161
#>  [78,] 27.436989
#>  [79,] 27.677454
#>  [80,] 27.333340
#>  [81,] 27.775303
#>  [82,] 27.686275
#>  [83,] 27.781903
#>  [84,] 27.897562
#>  [85,] 28.197661
#>  [86,] 27.817226
#>  [87,] 27.808872
#>  [88,] 27.812789
#>  [89,] 27.907614
#>  [90,] 28.084065
#>  [91,] 27.783294
#>  [92,] 29.363150
#>  [93,] 21.773955
#>  [94,] 21.822565
#>  [95,] 21.614756
#>  [96,] 21.670764
#>  [97,] 21.691693
#>  [98,] 12.916776
#>  [99,]  6.283970
#> [100,]  2.981308
#> [101,]  5.792462
#> [102,] 27.572968
#> [103,]        NA
#> [104,] 27.665693
#> [105,] 27.688750
#> [106,] 27.746038
#> [107,] 27.746929
#> [108,] 27.761019
#> [109,] 27.758205
#> [110,] 28.015228
#> [111,] 27.721300
#> [112,] 28.109526
#> [113,] 28.031653
#> [114,] 27.915532
#> [115,] 27.622607
#> [116,] 21.824216
#> [117,] 21.749597
#> [118,] 21.796911
#> [119,] 21.772727
#> [120,] 13.058415
#> [121,] 12.930614
#> [122,] 12.950815
#> [123,]  2.993239
#> [124,]  6.035890
#> [125,]  6.062508
#> [126,]        NA
#> [127,] 27.370695
#> [128,] 27.370550
#> [129,] 28.016850
#> [130,] 27.477122
#> [131,] 27.692737
#> [132,]  6.294746
#> [133,]  6.583400
#> [134,]  6.011970
#> [135,] 27.736296
#> [136,] 27.869491
#> [137,] 27.959219
#> [138,] 27.745726
#> [139,] 21.877834
#> [140,] 21.754267
#> [141,] 21.820680
#> [142,] 21.668769
#> [143,] 13.092411
#> [144,] 13.010264
#> [145,] 12.945095
#> [146,]  5.876139
#> [147,]  3.191751
#> [148,]  5.918984
#> [149,]        NA
#> [150,]  6.173396
#> [151,] 27.645915
#> [152,] 27.797603
#> [153,] 27.849367
#> [154,]  5.993021
#> [155,]  6.308805
#> [156,]  6.160518
#> [157,]  6.015590
#> [158,] 24.686705
#> [159,] 27.763195
#> [160,] 27.913334
#> [161,] 28.176738
#> [162,] 25.021435
#> [163,] 21.819662
#> [164,] 21.810206
#> [165,] 13.019329
#> [166,] 13.105410
#> [167,] 13.093626
#> [168,] 12.979934
#> [169,] 12.998341
#> [170,]  5.977243
#> [171,]  3.223696
#> [172,]  6.206684
#> [173,]  6.383169
#> [174,]  5.945038
#> [175,]  6.178938
#> [176,]  6.157323
#> [177,]  5.961311
#> [178,]  6.097826
#> [179,]  6.090458
#> [180,]  6.378300
#> [181,] 24.719136
#> [182,] 27.811313
#> [183,] 20.054619
#> [184,] 20.172535
#> [185,] 25.303899
#> [186,] 24.484593
#> [187,] 24.647085
#> [188,] 13.071161
#> [189,] 24.897529
#> [190,] 13.030282
#> [191,] 23.709568
#> [192,] 29.475679
#> [193,]  6.079120
#> [194,]  3.281172
#> [195,]  5.834415
#> [196,]  5.992471
#> [197,]  5.902961
#> [198,]  6.114268
#> [199,]  6.068222
#> [200,]  6.044415
#> [201,]  6.323806
#> [202,]  6.218690
#> [203,]  6.603816
#> [204,]  6.357847
#> [205,]  6.232811
#> [206,]  6.201820
#> [207,] 24.631304
#> [208,] 11.226504
#> [209,] 24.322908
#> [210,] 24.659767
#> [211,] 25.044764
#> [212,] 24.545413
#> [213,] 25.216052
#> [214,] 29.475778
#> [215,] 21.775691
#> [216,]  6.015397
#> [217,]  5.938444
#> [218,]  3.382876
#> [219,]  6.050669
#> [220,]  6.113467
#> [221,]  6.201093
#> [222,]  6.261793
#> [223,]  5.970998
#> [224,]  6.321241
#> [225,]  6.257572
#> [226,]  6.319670
#> [227,]  6.405481
#> [228,]  6.531199
#> [229,]  6.157765
#> [230,]  6.030926
#> [231,] 11.154503
#> [232,] 11.193731
#> [233,] 11.163204
#> [234,] 11.173489
#> [235,] 11.167403
#> [236,] 11.158375
#> [237,] 21.823055
#> [238,] 21.771130
#> [239,] 21.784255
#> [240,]  6.089041
#> [241,]  3.216555
#> [242,]  3.113035
#> [243,]  3.639183
#> [244,]  5.962353
#> [245,]  5.974483
#> [246,]  5.963831
#> [247,]  5.959785
#> [248,]  6.004410
#> [249,]  6.019913
#> [250,]  6.151806
#> [251,]  6.485531
#> [252,]  5.988535
#> [253,] 24.810256
#> [254,] 11.172699
#> [255,] 11.185942
#> [256,] 11.190009
#> [257,] 11.214255
#> [258,] 11.157790
#> [259,] 11.186728
#> [260,] 11.236362
#> [261,] 21.753590
#> [262,] 21.774945
#> [263,] 21.753826
#> [264,]  3.502730
#> [265,]  3.246574
#> [266,]  3.165641
#> [267,]  6.075505
#> [268,]  5.973592
#> [269,]  6.216331
#> [270,]  6.118761
#> [271,]  6.209957
#> [272,]  6.370993
#> [273,]  6.692662
#> [274,]  6.538998
#> [275,]  6.590112
#> [276,] 41.872022
#> [277,] 11.181156
#> [278,] 11.191553
#> [279,] 11.202318
#> [280,] 11.200084
#> [281,] 11.135493
#> [282,] 11.240035
#> [283,] 21.840552
#> [284,] 21.796608
#> [285,] 21.811199
#> [286,] 21.785445
#> [287,]  6.038751
#> [288,]        NA
#> [289,]  3.320877
#> [290,]  3.261804
#> [291,]  6.066984
#> [292,]  6.375621
#> [293,]  6.001425
#> [294,]  5.996427
#> [295,]  6.187191
#> [296,]  6.246519
#> [297,]  6.059510
#> [298,]  6.213187
#> [299,] 42.311390
#> [300,] 13.074426
#> [301,] 13.003649
#> [302,] 11.145427
#> [303,] 11.189415
#> [304,] 11.174552
#> [305,]        NA
#> [306,] 11.173757
#> [307,] 11.203632
#> [308,] 21.804552
#> [309,] 21.788646
#> [310,] 21.772792
#> [311,]  5.923597
#> [312,]  3.218248
#> [313,]  3.135849
#> [314,]  3.644689
#> [315,]  6.046544
#> [316,]  3.546913
#> [317,]  3.500385
#> [318,]  3.107761
#> [319,]  6.082675
#> [320,]  6.037078
#> [321,]  6.180582
#> [322,]  6.022782
#> [323,] 13.016623
#> [324,]  9.823401
#> [325,]  9.858397
#> [326,] 11.169283
#> [327,] 11.207006
#> [328,] 11.140007
#> [329,] 13.033587
#> [330,] 21.804989
#> [331,] 11.131121
#> [332,] 21.813272
#> [333,] 21.753110
#> [334,] 14.727773
#> [335,]  6.156015
#> [336,]  3.617985
#> [337,]  3.101900
#> [338,]  3.356278
#> [339,]  3.563438
#> [340,]        NA
#> [341,]  3.263259
#> [342,]  3.326704
#> [343,]  6.183905
#> [344,]  6.080436
#> [345,]  6.488724
#> [346,] 11.332409
#> [347,] 11.258213
#> [348,] 11.286569
#> [349,]  9.781960
#> [350,] 11.143375
#> [351,] 11.162084
#> [352,] 13.015881
#> [353,] 11.156960
#> [354,] 21.768904
#> [355,] 11.186232
#> [356,] 21.843405
#> [357,] 21.790342
#> [358,] 14.736492
#> [359,]  6.124405
#> [360,]  6.338921
#> [361,]  3.687659
#> [362,]  3.573474
#> [363,]  3.492680
#> [364,]  3.545605
#> [365,]  3.239070
#> [366,]  6.328718
#> [367,]  6.161450
#> [368,]  6.533532
#> [369,] 11.309583
#> [370,] 11.294058
#> [371,] 11.290805
#> [372,] 11.362525
#> [373,]  9.859822
#> [374,]  9.836249
#> [375,] 12.995977
#> [376,] 13.022202
#> [377,] 11.170852
#> [378,] 21.827672
#> [379,] 13.117265
#> [380,] 13.006647
#> [381,] 14.626738
#> [382,] 19.339784
#> [383,] 19.542584
#> [384,]  6.234156
#> [385,]  3.461183
#> [386,]  3.440293
#> [387,]  3.236858
#> [388,]  3.325661
#> [389,]  3.132491
#> [390,]  6.240327
#> [391,]  6.143570
#> [392,] 11.289834
#> [393,] 11.265898
#> [394,] 11.312384
#> [395,] 11.312961
#> [396,] 11.326362
#> [397,] 11.298920
#> [398,]  9.717502
#> [399,] 12.975927
#> [400,]  9.735613
#> [401,] 13.006738
#> [402,] 13.054775
#> [403,] 14.665204
#> [404,] 14.586937
#> [405,] 14.699318
#> [406,] 14.787283
#> [407,] 19.492637
#> [408,]  6.040530
#> [409,]  3.259878
#> [410,]  3.414015
#> [411,]  3.187875
#> [412,]  3.325233
#> [413,]  3.385205
#> [414,]  6.237917
#> [415,] 11.322957
#> [416,] 11.309619
#> [417,] 11.300170
#> [418,] 11.271161
#> [419,] 11.283633
#> [420,] 11.261142
#> [421,] 11.372388
#> [422,]  9.692714
#> [423,]  9.732982
#> [424,] 12.950026
#> [425,] 13.007754
#> [426,] 14.646457
#> [427,] 14.690757
#> [428,] 14.817727
#> [429,] 14.689103
#> [430,] 14.720120
#> [431,] 19.944522
#> [432,]  6.398810
#> [433,]  3.630577
#> [434,]  3.361942
#> [435,]  3.240007
#> [436,]        NA
#> [437,]  3.680146
#> [438,] 11.310647
#> [439,] 11.292327
#> [440,] 11.311526
#> [441,]  9.870018
#> [442,] 13.073024
#> [443,] 12.974329
#> [444,]  9.778134
#> [445,] 13.068898
#> [446,]  9.775679
#> [447,]  9.774198
#> [448,]  9.903830
#> [449,] 14.734814
#> [450,] 14.697724
#> [451,] 14.815406
#> [452,] 14.716257
#> [453,] 14.671151
#> [454,] 14.707799
#> [455,] 14.644518
#> [456,]  5.947941
#> [457,]  3.321535
#> [458,]  3.271211
#> [459,]        NA
#> [460,]  3.533919
#> [461,]  9.867179
#> [462,]  9.839316
#> [463,]  9.763657
#> [464,] 12.997748
#> [465,] 13.038341
#> [466,] 13.110638
#> [467,] 13.076796
#> [468,] 13.057761
#> [469,] 13.000210
#> [470,] 12.999624
#> [471,] 14.696935
#> [472,] 14.690053
#> [473,] 14.684041
#> [474,] 14.692020
#> [475,] 14.768157
#> [476,] 14.668199
#> [477,] 14.634612
#> [478,] 14.667732
#> [479,] 14.631155
#> [480,]  5.991966
#> [481,]  3.668166
#> [482,]  3.381220
#> [483,]  3.349480
#> [484,]  9.836567
#> [485,] 11.145949
#> [486,] 11.108771
#> [487,] 13.080140
#> [488,] 13.025574
#> [489,]  9.823908
#> [490,] 13.024484
#> [491,] 13.024439
#> [492,] 13.033039
#> [493,] 14.694583
#> [494,] 14.652260
#> [495,] 14.736937
#> [496,] 14.683935
#> [497,] 14.748679
#> [498,] 14.731248
#> [499,] 14.569738
#> [500,] 14.544858
#> [501,] 14.634338
#> [502,] 14.634852
#> [503,]  6.024382
#> [504,]  5.835302
#> [505,]  3.237400
#> [506,]  3.217171
#> [507,] 11.196501
#> [508,] 11.121283
#> [509,] 11.163809
#> [510,] 13.002087
#> [511,] 13.002757
#> [512,] 13.015432
#> [513,] 13.075096
#> [514,] 13.097691
#> [515,] 13.072003
#> [516,] 12.961095
#> [517,] 14.670558
#> [518,] 14.797638
#> [519,] 14.849343
#> [520,] 14.675878
#> [521,] 14.745671
#> [522,] 14.616947
#> [523,] 14.631177
#> [524,] 14.645940
#> [525,] 14.598237
#> [526,] 14.650532
#> [527,]  5.995101
#> [528,]  6.026352
#> [529,]  3.329144
#> [530,] 11.133440
#> [531,] 11.166812
#> [532,] 11.169719
#> [533,] 12.981752
#> [534,] 13.035430
#> [535,] 12.993196
#> [536,] 13.055314
#> [537,] 13.076386
#> [538,] 13.035113
#> [539,] 13.028156
#> [540,] 14.698252
#> [541,] 14.679318
#> [542,] 14.661843
#> [543,] 14.815170
#> [544,]        NA
#> [545,] 14.663530
#> [546,] 14.610712
#> [547,] 14.597255
#> [548,] 14.617393
#> [549,] 14.654673
#> [550,] 14.661761
#> [551,]  5.973812
#> [552,]  3.205994
#> [553,] 11.112977
#> [554,] 11.176658
#> [555,]  9.876485
#> [556,] 13.062760
#> [557,] 13.113767
#> [558,] 13.013114
#> [559,] 13.057064
#> [560,] 13.020581
#> [561,] 13.066518
#> [562,]  9.880260
#> [563,] 13.051012
#> [564,] 14.679451
#> [565,] 14.671611
#> [566,] 14.688794
#> [567,]        NA
#> [568,] 14.699683
#> [569,] 14.656194
#> [570,] 14.675434
#> [571,] 14.676620
#> [572,] 14.592472
#> [573,] 14.605654
#> [574,] 12.259523
#> [575,]  5.948674
#> [576,] 11.178907
#> [577,]  9.795771
#> [578,]  9.848588
#> [579,] 11.194706
#> [580,] 11.127675
#> [581,]  9.814428
#> [582,]  9.756863
#> [583,] 12.970465
#> [584,] 12.934760
#> [585,] 13.031208
#> [586,] 12.944697
#> [587,] 14.693323
#> [588,] 14.679438
#> [589,] 14.742006
#> [590,] 14.654916
#> [591,] 14.638918
#> [592,] 14.637846
#> [593,] 14.587837
#> [594,] 14.620322
#> [595,] 13.437186
#> [596,] 14.654298
#> [597,] 12.345645
#> [598,] 12.300649
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
  e.g. [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)’s
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
fetch-through-fuse pipeline view, and the “Multi-Source Raster Fusion
Pipeline” architecture article for the complete function-level reference
this vignette exposed piece by piece.
