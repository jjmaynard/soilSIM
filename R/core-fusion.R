#' @title Bayesian Updating: Conjugate Fusion and General Grid-KDE Fusion
#'
#' @description Scalar (plain-vector) toolkit for fusing a prior belief
#'   distribution with a likelihood (e.g. an SSURGO-derived prior against
#'   field-measured data).
#'
#'   These are standalone building blocks. The raster fusion pipeline in
#'   `core-fusion.R` supplies the equivalent primitives used in the
#'   production SSURGO x SOLUS100 workflow; the tabular Monte Carlo pipeline
#'   in `core-montecarlo.R` does not call these functions directly.
#'
#'   Three tiers, in increasing order of generality/cost:
#'   1. `fuse_normal_normal()` - exact, closed-form, Normal-only.
#'   2. `fuse_beta()`/`fuse_gamma()` - exact, closed-form, same-family-only
#'      (each family's kernel parameters simply ADD - the same principle
#'      `fuse_normal_normal()` uses for Normal's natural parameters,
#'      applied to Beta/Gamma), with a `feasible` flag and a documented
#'      moment-based fallback route when infeasible.
#'   3. `update_prior()` - fully general grid-KDE Bayes' rule (any shape,
#'      any family, even mismatched sides), at the cost of needing raw
#'      samples (not distributional parameters) on both sides.
#'
#'   `fuse_bivariate_normal()`/`fuse_texture_group_from_triplets()` extend
#'   tier 1 to the joint 3-part-composition ILR case (clay/sand/silt texture,
#'   in whichever role order `core-montecarlo.R`'s `composition_groups$texture$members`
#'   configures - see `core-distributions.R`'s ILR section header), since
#'   independent per-fraction fusion (e.g. three separate `fuse_beta()`
#'   calls) measurably breaks sum-to-100 (up to about 10 percentage points on
#'   realistic synthetic data).
#' @name bayesian_updating
NULL

# ==============================================================================
# Tier 1: Normal-Normal conjugate fusion
# ==============================================================================

#' Precision-weighted conjugate Normal-Normal posterior
#'
#' @param prior_mu,prior_sigma Prior mean/sd (numeric scalar or vector).
#' @param lik_mu,lik_sigma Likelihood mean/sd (numeric scalar or vector, same
#'   length as the prior args).
#' @return `list(mu = posterior mean, sigma = posterior sd)`.
#'
#' @details Invariants (see `tests/testthat/test-bayesian-updating.R`):
#'   posterior `sigma <= min(prior_sigma, lik_sigma)`; posterior `mu` is
#'   bounded between `prior_mu` and `lik_mu`; symmetric under swapping
#'   `(prior_mu, prior_sigma)` <-> `(lik_mu, lik_sigma)`.
#' @family bayesian-fusion
#' @export
fuse_normal_normal <- function(prior_mu, prior_sigma, lik_mu, lik_sigma) {
  prior_prec <- 1 / (prior_sigma^2)
  lik_prec   <- 1 / (lik_sigma^2)
  post_prec  <- prior_prec + lik_prec
  post_sigma <- sqrt(1 / post_prec)
  post_mu    <- (prior_mu * prior_prec + lik_mu * lik_prec) / post_prec
  list(mu = post_mu, sigma = post_sigma)
}

#' Moment-match a raw-space Lognormal's mean/sd onto its underlying Normal's mu/sigma
#'
#' Needed because several properties (e.g. SOC, CEC) are strictly positive
#' and right-skewed - fusing their raw mean/sd as if Normal lets the
#' posterior put real probability mass below zero, which is physically
#' impossible. `fuse_normal_normal()` is reused in log-space instead.
#'
#' @param mu,sigma Raw-space mean/sd.
#' @return `list(mu = log-space mu, sigma = log-space sigma)`.
#' @family bayesian-fusion
#' @export
convert_normal_to_lognormal <- function(mu, sigma) {
  sigma_log <- sqrt(log(1 + (sigma / mu)^2))
  mu_log <- log(mu) - sigma_log^2 / 2
  list(mu = mu_log, sigma = sigma_log)
}

#' Inverse of `convert_normal_to_lognormal()`
#'
#' Converts log-space Normal(mu, sigma) parameters back to the raw-space
#' Lognormal's mean/sd.
#' @param mu_log,sigma_log Log-space mean/sd.
#' @return `list(mu = raw-space mean, sigma = raw-space sd)`.
#' @family bayesian-fusion
#' @export
convert_lognormal_to_normal <- function(mu_log, sigma_log) {
  mu <- exp(mu_log + sigma_log^2 / 2)
  sigma <- sqrt((exp(sigma_log^2) - 1) * exp(2 * mu_log + sigma_log^2))
  list(mu = mu, sigma = sigma)
}

# ==============================================================================
# Tier 2: Beta-Beta and Gamma-Gamma conjugate fusion
# ==============================================================================

#' Fuse two independent Beta belief distributions via density multiplication
#'
#' `Beta(a1,b1) x Beta(a2,b2) -> Beta(a1+a2-1, b1+b2-1)`, derived from the
#' Beta kernel `x^(a-1)(1-x)^(b-1)`. Requires `a1+a2 > 1` and `b1+b2 > 1` to
#' remain a valid distribution; `feasible` flags where that fails.
#'
#' @param prior_alpha,prior_beta,lik_alpha,lik_beta Numeric scalars or vectors.
#' @return `list(alpha=, beta=, feasible=)`.
#' @family bayesian-fusion
#' @export
fuse_beta <- function(prior_alpha, prior_beta, lik_alpha, lik_beta) {
  alpha <- prior_alpha + lik_alpha - 1
  beta <- prior_beta + lik_beta - 1
  feasible <- (prior_alpha + lik_alpha > 1) & (prior_beta + lik_beta > 1)
  list(alpha = alpha, beta = beta, feasible = feasible)
}

#' Fuse two independent Gamma (shape/rate) belief distributions
#'
#' `Gamma(k1,r1) x Gamma(k2,r2) -> Gamma(k1+k2-1, r1+r2)`, derived from the
#' Gamma kernel `x^(k-1)exp(-r*x)`. Requires `k1+k2 > 1` to remain valid.
#'
#' @param prior_shape,prior_rate,lik_shape,lik_rate Numeric scalars or vectors.
#' @return `list(shape=, rate=, feasible=)`.
#' @family bayesian-fusion
#' @export
fuse_gamma <- function(prior_shape, prior_rate, lik_shape, lik_rate) {
  shape <- prior_shape + lik_shape - 1
  rate <- prior_rate + lik_rate
  feasible <- (prior_shape + lik_shape > 1)
  list(shape = shape, rate = rate, feasible = feasible)
}

#' Method-of-moments Gamma fit
#'
#' Used by the closed-form same-family fusion route's infeasible-cell
#' fallback (re-express as Normal moments, fuse as Normal, convert back).
#' @param mean,var Mean/variance to convert.
#' @return `list(shape=, rate=)`.
#' @family bayesian-fusion
#' @export
convert_moments_to_gamma <- function(mean, var) {
  list(shape = mean^2 / var, rate = mean / var)
}

#' Method-of-moments Beta fit
#'
#' Used by the closed-form same-family fusion route's infeasible-cell
#' fallback (re-express as Normal moments, fuse as Normal, convert back).
#' @param mean,var Mean/variance to convert.
#' @return `list(alpha=, beta=)`.
#' @family bayesian-fusion
#' @export
convert_moments_to_beta <- function(mean, var) {
  common <- mean * (1 - mean) / var - 1
  list(alpha = mean * common, beta = (1 - mean) * common)
}

#' Beta's mean/variance as a function of its own (alpha, beta)
#'
#' The inverse direction of `convert_moments_to_beta()` - used to re-express an
#' infeasible same-family fusion's inputs as Normal moments before falling
#' back to `fuse_normal_normal()`.
#' @param alpha,beta Beta shape parameters.
#' @return `list(mean=, var=)`.
#' @family bayesian-fusion
#' @export
convert_beta_to_moments <- function(alpha, beta) {
  list(mean = alpha / (alpha + beta), var = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
}

#' Gamma's mean/variance as a function of its own (shape, rate)
#'
#' The inverse direction of `convert_moments_to_gamma()` - used to re-express an
#' infeasible same-family fusion's inputs as Normal moments before falling
#' back to `fuse_normal_normal()`.
#' @param shape,rate Gamma shape/rate parameters.
#' @return `list(mean=, var=)`.
#' @family bayesian-fusion
#' @export
convert_gamma_to_moments <- function(shape, rate) {
  list(mean = shape / rate, var = shape / rate^2)
}

#' Fuse a prior and likelihood belief distribution of the same family
#'
#' @param prior_params,lik_params Named lists of that family's parameters:
#'   `list(mu=,sigma=)` for `"normal"`, `list(alpha=,beta=)` for `"beta"`,
#'   `list(shape=,rate=)` for `"gamma"`.
#' @param family One of `"normal"`, `"beta"`, `"gamma"`. Both sides must
#'   already be fit in this same family - there is no closed form for fusing
#'   mismatched families; use `update_prior()` for that.
#' @return The matching `fuse_*()` function's return value.
#' @family bayesian-fusion
#' @export
fuse_distribution <- function(prior_params, lik_params, family = c("normal", "beta", "gamma")) {
  family <- match.arg(family)
  switch(family,
    normal = fuse_normal_normal(prior_params$mu, prior_params$sigma, lik_params$mu, lik_params$sigma),
    beta = fuse_beta(prior_params$alpha, prior_params$beta, lik_params$alpha, lik_params$beta),
    gamma = fuse_gamma(prior_params$shape, prior_params$rate, lik_params$shape, lik_params$rate)
  )
}

# ==============================================================================
# Tier 3: fully general grid-KDE Bayesian updating
# ==============================================================================

#' Combine a prior and likelihood distribution into a posterior via grid-based KDE
#'
#' Estimates prior and likelihood densities over a shared value grid via
#' kernel density estimation, multiplies them per Bayes' rule, and samples
#' from the resulting (normalized) posterior. Fully general - no
#' distributional family assumption, no requirement that both sides match -
#' at the cost of needing raw samples (not distributional parameters) and
#' being considerably more expensive than the closed-form tiers above.
#'
#' @param prior_distribution,likelihood_distribution Numeric vectors of samples.
#' @param grid_range Optional length-2 vector giving the grid's min/max. If
#'   `NULL`, derived from the combined range of both inputs, padded by 1.
#' @param grid_resolution Step size of the evaluation grid.
#' @param n Number of posterior samples to draw.
#' @param winsorize_probs Optional length-2 vector `c(lower, upper)` (e.g. `c(0.01, 0.99)`). When
#'   supplied, each side's sample is independently clipped to its OWN `quantile()` range at these
#'   probabilities before the `stats::density(bw = "nrd0")` call. `NULL` (default) preserves this
#'   function's original behavior exactly.
#' @param posterior_probs Optional numeric vector of probabilities (0-1) at which to report
#'   posterior percentiles (e.g. `c(0.05, 0.5, 0.95)`). `NULL` (default) preserves this function's
#'   original return shape EXACTLY - a plain numeric vector of `n` posterior samples - so every
#'   existing caller (`fuse_property()`, `R/core-montecarlo.R`'s `fuse_observed_data_into_priors()`,
#'   this file's own tests) is completely unaffected. When supplied, the return shape changes to
#'   the richer list described below (see `@return`); see `@section Grid-based percentiles` for why
#'   this is the primary use case for `posterior_probs`, not `quantile()` on the resampled vector.
#' @return When `posterior_probs` is `NULL` (default): a numeric vector of `n` samples drawn from
#'   the posterior distribution (unchanged from this function's original contract). When
#'   `posterior_probs` is supplied: a list with -
#'   - `samples`: the same `n`-length resampled vector as the `NULL` case (still available for
#'     callers that want actual posterior draws, e.g. downstream Monte Carlo propagation).
#'   - `mean`, `var`: the posterior mean/variance computed EXACTLY from the discretized
#'     `posterior_prob`/`value_grid` distribution (`sum(value_grid * posterior_prob)` and its
#'     variance analogue) - zero resampling noise, strictly at least as accurate as `mean(samples)`/
#'     `var(samples)` for any `n`.
#'   - `percentiles`: a named numeric vector (names `"P5"`, `"P50"`, ... matching `posterior_probs`),
#'     read directly off `cumsum(posterior_prob)` - also zero resampling noise, accurate to within
#'     `grid_resolution`.
#'   - `value_grid`, `posterior_prob`: the underlying discretized distribution itself, for callers
#'     that want to compute something these three summaries don't cover.
#'
#' @section Grid-based percentiles have zero resampling noise:
#' This function already computes the FULL discretized posterior distribution
#' (`posterior_prob` over `value_grid`) before its final step, which historically only existed to
#' hand back a plain vector (`sample(value_grid, size = n, prob = posterior_prob, replace = TRUE)`).
#' Reading percentiles/moments directly off `posterior_prob`'s cumulative distribution - finding
#' where `cumsum(posterior_prob)` crosses each target probability - is exact relative to
#' `grid_resolution` and has no resampling noise at all, unlike `quantile(samples, probs = ...)`
#' would have. No amount of increasing `n` improves on reading percentiles directly off an
#' already-fully-known discrete distribution, so this is `posterior_probs`' primary computation
#' path, not a convenience wrapper around the resampled `samples` vector.
#'
#' @section Known limitation (raw-draws fusion vs. percentile-reconstruction fusion):
#' A dedicated adversarial test found that for
#' a skewed prior fused against a weak/wide likelihood, `soilSIM`'s `prior_fusion_method =
#' "raw_draws"` route (a genuine empirical resample as the prior side) can produce a fused posterior
#' mean measurably *less* close to the prior population's true mean than the default
#' percentile-reconstruction route, in a way that first looked like a bandwidth-selection bug. Task
#' P3.1's follow-up investigation ruled that out directly (forcing both routes onto the exact same
#' `stats::density()` bandwidth left the gap essentially unchanged) and traced the real mechanism to
#' `R/core-distributions.R`'s `sim_linear_cdf_batch()`: reconstructing a distribution from only 5
#' percentile knots is an inherent information-loss approximation for skewed shapes (confirmed
#' method-agnostic - `linear_cdf`/`spline`/`kde` reconstruction all show a comparable P75-P95
#' overshoot on the same synthetic scenario), not a defect in this function or in raw_draws. See
#' that task's write-up for the full investigation; `winsorize_probs` remains available below as a
#' harmless, independently-useful safety net for genuine extreme-outlier inputs, not as a fix for
#' this finding (there was nothing here to fix).
#' @family bayesian-fusion
#' @export
update_prior <- function(prior_distribution, likelihood_distribution, grid_range = NULL,
                             grid_resolution = 0.01, n = 1000, winsorize_probs = NULL,
                             posterior_probs = NULL) {
  if (!is.null(winsorize_probs)) {
    clip <- function(x) {
      bounds <- stats::quantile(x, probs = winsorize_probs, na.rm = TRUE, names = FALSE)
      pmin(pmax(x, bounds[1]), bounds[2])
    }
    prior_distribution <- clip(prior_distribution)
    likelihood_distribution <- clip(likelihood_distribution)
  }

  if (is.null(grid_range)) {
    grid_min <- min(c(prior_distribution, likelihood_distribution)) - 1
    grid_max <- max(c(prior_distribution, likelihood_distribution)) + 1
  } else {
    grid_min <- grid_range[1]
    grid_max <- grid_range[2]
  }

  value_grid <- seq(grid_min, grid_max, by = grid_resolution)

  prior_density <- density(prior_distribution, bw = "nrd0", from = grid_min, to = grid_max, n = length(value_grid))
  likelihood_density <- density(likelihood_distribution, bw = "nrd0", from = grid_min, to = grid_max, n = length(value_grid))

  prior_prob <- approxfun(prior_density$x, prior_density$y)(value_grid)
  likelihood_prob <- approxfun(likelihood_density$x, likelihood_density$y)(value_grid)

  prior_prob <- prior_prob / sum(prior_prob)
  likelihood_prob <- likelihood_prob / sum(likelihood_prob)

  posterior_prob <- prior_prob * likelihood_prob
  posterior_prob <- posterior_prob / sum(posterior_prob)

  samples <- sample(value_grid, size = n, prob = posterior_prob, replace = TRUE)

  if (is.null(posterior_probs)) {
    return(samples)
  }

  # Grid-based exact percentiles/moments - see @section Grid-based percentiles above for why this
  # is preferred over quantile()/mean()/var() on `samples`.
  cdf <- cumsum(posterior_prob)
  idx <- pmin(findInterval(posterior_probs, cdf) + 1L, length(value_grid))
  percentiles <- stats::setNames(value_grid[idx], paste0("P", round(posterior_probs * 100)))

  mean_exact <- sum(value_grid * posterior_prob)
  var_exact <- sum((value_grid - mean_exact)^2 * posterior_prob)

  list(
    samples = samples,
    mean = mean_exact,
    var = var_exact,
    percentiles = percentiles,
    value_grid = value_grid,
    posterior_prob = posterior_prob
  )
}

# ==============================================================================
# Joint (compositional) fusion for texture (clay/sand/silt)
# ==============================================================================

#' Fuse two independent bivariate Normal beliefs about the same 2D quantity
#'
#' Precision MATRICES add - the direct multivariate generalization of
#' `fuse_normal_normal()`'s scalar-precision addition.
#'
#' @param mu1,mu2 Length-2 mean vectors.
#' @param Sigma1,Sigma2 2x2 covariance matrices.
#' @return `list(mu = fused mean vector, Sigma = fused covariance matrix)`.
#' @family bayesian-fusion
#' @export
fuse_bivariate_normal <- function(mu1, Sigma1, mu2, Sigma2) {
  P1 <- solve(Sigma1)
  P2 <- solve(Sigma2)
  Sigma_post <- solve(P1 + P2)
  list(mu = as.numeric(Sigma_post %*% (P1 %*% mu1 + P2 %*% mu2)), Sigma = Sigma_post)
}

#' Fuse two full clay/sand/silt low-rep-high triplets jointly via ILR fusion
#'
#' Wraps `core-distributions.R`'s `estimate_ilr_moments_mc()` (once per side) +
#' `fuse_bivariate_normal()` + `sample_ilr_posterior()` for future joint
#' texture-prior updating. Independent per-fraction fusion (e.g. three
#' separate `fuse_beta()` calls) measurably breaks sum-to-100 (documented
#' upstream: up to 10.5 percentage points on realistic synthetic data) -
#' this is why the fusion happens jointly, in ILR space, instead.
#'
#' @param prior_triplets,lik_triplets Named list with elements `clay`,
#'   `sand`, `silt`, each a length-3 numeric vector `c(low, rep, high)`. These
#'   keys are positional-role placeholders matching `estimate_ilr_moments_mc()`'s
#'   `low_clay`/`low_sand`/`low_silt` parameter naming (position 1/2/3 of the
#'   ILR sequential binary partition - see `core-distributions.R`'s ILR section
#'   header), not an identity requirement; callers (e.g.
#'   `core-montecarlo.R`'s `fuse_observed_data_into_priors()`) build them from
#'   `composition_groups$texture$members` in configured order.
#' @param z_prior,z_lik Standard-normal quantile matching each side's
#'   low/high interval (e.g. `qnorm(0.95)` for a 5th/95th-percentile
#'   low/high).
#' @param n_mc Monte Carlo sample size passed to `estimate_ilr_moments_mc()`.
#' @param n_samples Number of posterior composition samples to draw.
#' @param total Composition target sum (default 100).
#' @return `list(posterior_samples = <n_samples x 3 matrix, columns clay/sand/silt,
#'   guaranteed sum-to-`total` by construction>, ilr_mu=, ilr_Sigma=)`.
#' @keywords internal
fuse_texture_group_from_triplets <- function(prior_triplets, lik_triplets, z_prior, z_lik,
                                              n_mc = 2000, n_samples = 1000, total = 100) {
  prior_moments <- estimate_ilr_moments_mc(
    low_clay = prior_triplets$clay[1], rep_clay = prior_triplets$clay[2], high_clay = prior_triplets$clay[3],
    low_sand = prior_triplets$sand[1], rep_sand = prior_triplets$sand[2], high_sand = prior_triplets$sand[3],
    low_silt = prior_triplets$silt[1], rep_silt = prior_triplets$silt[2], high_silt = prior_triplets$silt[3],
    z = z_prior, n_mc = n_mc
  )
  lik_moments <- estimate_ilr_moments_mc(
    low_clay = lik_triplets$clay[1], rep_clay = lik_triplets$clay[2], high_clay = lik_triplets$clay[3],
    low_sand = lik_triplets$sand[1], rep_sand = lik_triplets$sand[2], high_sand = lik_triplets$sand[3],
    low_silt = lik_triplets$silt[1], rep_silt = lik_triplets$silt[2], high_silt = lik_triplets$silt[3],
    z = z_lik, n_mc = n_mc
  )

  fused <- fuse_bivariate_normal(prior_moments$mu, prior_moments$Sigma, lik_moments$mu, lik_moments$Sigma)
  posterior_samples <- sample_ilr_posterior(fused$mu, fused$Sigma, n = n_samples, total = total)

  list(posterior_samples = posterior_samples, ilr_mu = fused$mu, ilr_Sigma = fused$Sigma)
}

# ==============================================================================
# Tabular dispatch entry point
# ==============================================================================

#' Fuse a prior and likelihood belief distribution, dispatching by input shape
#'
#' Unlike the raster-native reference (`fuse_adaptive()`), which dispatches
#' between the general and closed-form routes by AOI cell count - a concept
#' with no tabular analogue - this dispatches on the SHAPE of `prior`/
#' `likelihood`: atomic numeric vectors (raw samples) route to the fully
#' general `update_prior()`; named lists of family-native parameters
#' route to the closed-form `fuse_distribution()`. This removes the "which
#' heuristic" ambiguity entirely, since the caller's input shape already
#' determines which route applies.
#'
#' @param prior,likelihood Either numeric vectors of raw samples (both sides
#'   must be vectors), or named lists of family-native parameters matching
#'   `family` (both sides must be lists).
#' @param family Required when `prior`/`likelihood` are parameter lists; one
#'   of `"normal"`, `"beta"`, `"gamma"`. Ignored (and inferred as "general")
#'   when `prior`/`likelihood` are raw-sample vectors.
#' @param bounds Unused currently; kept for interface symmetry with
#'   percentile-triplet-based callers.
#' @param method Optional assertion/override: `"general"` or `"closed_form"`.
#'   Errors if it doesn't match what the input shape implies, rather than
#'   silently overriding it.
#' @param n_samples,grid_resolution Passed to `update_prior()` for the
#'   general route.
#' @return For the general route: a numeric vector of posterior samples
#'   (`update_prior()`'s native output). For the closed-form route: the
#'   family-native posterior parameter list (`fuse_distribution()`'s native output).
#'   These are NOT the same shape - documented deliberately rather than
#'   forced into a fake-uniform contract, since the caller already knows
#'   which shape it's passing in and therefore which shape it gets back.
#' @family bayesian-fusion
#' @export
fuse_property <- function(prior, likelihood, ...) UseMethod("fuse_property")

#' @rdname fuse_property
#' @export
fuse_property.default <- function(prior, likelihood, family = NULL, bounds = NULL, method = NULL,
                                   n_samples = 1000, grid_resolution = 0.01, ...) {
  prior_is_vector <- is.atomic(prior) && is.numeric(prior)
  lik_is_vector <- is.atomic(likelihood) && is.numeric(likelihood)

  if (prior_is_vector != lik_is_vector) {
    stop("fuse_property(): prior and likelihood must be the SAME shape (both raw-sample vectors, or both parameter lists).")
  }

  resolved_method <- if (prior_is_vector) "general" else "closed_form"
  if (!is.null(method) && !identical(method, resolved_method)) {
    stop(sprintf(
      "fuse_property(): method = '%s' was requested but the input shape implies '%s' - pass matching inputs or omit method.",
      method, resolved_method
    ))
  }

  if (resolved_method == "general") {
    return(update_prior(prior, likelihood, grid_resolution = grid_resolution, n = n_samples))
  }

  if (is.null(family)) {
    stop("fuse_property(): family is required when prior/likelihood are parameter lists (closed_form route).")
  }
  fuse_distribution(prior, likelihood, family = family)
}

#' Tag a list of percentile-value rasters for `fuse_property()`'s raster method
#'
#' `fuse_property()` dispatches on the class of `prior`: a plain numeric vector or plain list
#' goes to [fuse_property.default()] (scalar/vector fusion); a list tagged with this class goes to
#' [fuse_property.percentile_rasters()] (raster-native fusion, forwarding to
#' [fuse_property_adaptive()] unchanged). A plain, unclassed `list()` of percentile-value rasters
#' cannot be distinguished from a closed-form parameter list by class alone, hence the explicit tag.
#'
#' @param x A list of `terra::SpatRaster` percentile-value layers.
#' @return `x`, with class `c("percentile_rasters", class(x))`.
#' @export
new_percentile_rasters <- function(x) {
  stopifnot(is.list(x))
  class(x) <- unique(c("percentile_rasters", class(x)))
  x
}

#' @rdname fuse_property
#' @param ... Passed to the dispatched method.
#' @param prior_probs,lik_probs Percentile levels for `prior`/`likelihood` (raster method only).
#' @param property_config A property configuration list selecting the fusion distribution family
#'   (raster method only) - see [fuse_property_adaptive()].
#' @param threshold_cells Raster method only; see [fuse_property_adaptive()].
#' @export
fuse_property.percentile_rasters <- function(prior, likelihood, prior_probs, lik_probs,
                                              property_config, threshold_cells = 80000, ...) {
  fuse_property_adaptive(unclass(prior), prior_probs, unclass(likelihood), lik_probs,
                          property_config, threshold_cells = threshold_cells, ...)
}


# ============================================================================
# merged from core-fusion.R (P1 file reorg)
# ============================================================================

#' @title Raster-Native Bayesian Fusion (Prior x Likelihood, per Cell)
#'
#' @description Per-cell Bayesian fusion of a prior belief distribution
#'   (e.g. an SSURGO-derived percentile raster set) with a likelihood (e.g.
#'   an independent percentile raster set from another source), across a
#'   whole `terra::SpatRaster` AOI at once. This is the raster-native
#'   counterpart of `R/core-fusion.R`.
#'
#'   Uses `core-fusion.R`'s scalar/vector fusion functions
#'   directly - `fuse_normal_normal()`,
#'   `fuse_beta()`, `fuse_gamma()`, `convert_moments_to_gamma()`/`convert_moments_to_beta()`/
#'   `convert_beta_to_moments()`/`convert_gamma_to_moments()`, `convert_normal_to_lognormal()`/
#'   `convert_lognormal_to_normal()`, `update_prior()`, and
#'   `R/core-distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`, plus
#'   `fuse_bivariate_normal()` - are all pure elementwise arithmetic, so they
#'   work unchanged on `SpatRaster` inputs.
#'
#' @section Entry points:
#' `run_fusion()`/`run_fusion_group()` are the top-level
#' orchestrators, wiring the fetch-and-cache layer (`R/cache.R`,
#' `R/adapter-ssurgo-simulate.R`, `R/adapter-solus.R`) to the fusion core.
#' `fuse_property_adaptive()` is the lower-level entry point for callers
#' who already have their own pre-fetched `prior_value_rasters`/
#' `lik_value_rasters` and want to skip the fetch-and-cache wrapper.
#'
#' @section Configuration:
#' `group_members()`/`fuse_property_adaptive()` take a per-call config
#' list/member vector directly, reusing soilSIM's
#' `config$monte_carlo$composition_groups` convention (see
#' `R/core-distributions.R`'s `resolve_composition_groups()`).
#' @name raster_fusion
NULL

# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

#' Method-of-moments Gamma fit from a list of percentile-value rasters
#'
#' Thin raster wrapper around the existing scalar `convert_moments_to_gamma()`
#' (`core-fusion.R`) - the mean/variance computation is the only
#' genuinely raster-specific part (combining a *list* of rasters via `Reduce()`).
#' @param value_rasters List of percentile-value SpatRasters.
#' @return `list(shape = SpatRaster, rate = SpatRaster)`.
#' @family raster-fusion
#' @export
fit_gamma_mom_raster <- function(value_rasters) {
  k <- length(value_rasters)
  mean_r <- Reduce(`+`, value_rasters) / k
  var_r <- Reduce(`+`, lapply(value_rasters, function(r) (r - mean_r)^2)) / k
  convert_moments_to_gamma(mean_r, var_r)
}

#' Locate the low/median/high indices and probabilities within a percentile
#' probability vector.
#' @param probs Numeric probability vector.
#' @return `list(lo_idx=, p50_idx=, hi_idx=, p_lo=, p_hi=)`.
#' @keywords internal
percentile_index <- function(probs) {
  p50_idx <- which(probs == 0.5)
  tail_idx <- probs > 0 & probs < 1 & probs != 0.5
  p_lo <- min(probs[tail_idx]); p_hi <- max(probs[tail_idx])
  list(lo_idx = which(probs == p_lo), p50_idx = p50_idx, hi_idx = which(probs == p_hi),
       p_lo = p_lo, p_hi = p_hi)
}

#' Interior percentiles only (excludes p=0/p=1, where the metalog logit
#' transform is undefined) - matches `fit_metalog_linear_raster()`'s own restriction.
#' @param probs Numeric probability vector.
#' @return `list(idx=, probs=)`.
#' @keywords internal
percentile_interior <- function(probs) {
  idx <- which(probs > 0 & probs < 1)
  list(idx = idx, probs = probs[idx])
}

#' Re-express one side's percentile-value rasters onto the other side's
#' probability grid when they differ (e.g. one source's 5th/95th percentiles
#' vs another's 2.5th/97.5th) - every fusion route below fits both sides
#' using ONE shared `probs` vector, so passing mismatched probabilities
#' through unchanged would silently mis-fit whichever side's actual
#' percentiles don't match the assumed labels.
#'
#' Uses `quantile_linear_cdf_raster()` to reinterpolate the wider-ranging
#' side onto the narrower side's exact probabilities (the narrower side is
#' always safely interpolatable FROM the wider one, never the reverse,
#' since `quantile_linear_cdf_raster()` requires the target probability to
#' fall within the source's own range).
#' @param prior_value_rasters,lik_value_rasters Lists of percentile-value SpatRasters.
#' @param prior_probs,lik_probs Matching probability vectors.
#' @return `list(prior_value_rasters=, lik_value_rasters=, probs=)`.
#' @family raster-fusion
#' @export
align_percentile_probs <- function(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs) {
  if (identical(prior_probs, lik_probs)) {
    return(list(prior_value_rasters = prior_value_rasters, lik_value_rasters = lik_value_rasters, probs = prior_probs))
  }

  target_probs <- if (diff(range(prior_probs)) <= diff(range(lik_probs))) prior_probs else lik_probs
  reinterpolate <- function(value_rasters, probs) {
    if (identical(probs, target_probs)) return(value_rasters)
    if (min(target_probs) < min(probs) || max(target_probs) > max(probs)) {
      stop("align_percentile_probs(): target probabilities fall outside this side's own percentile range - cannot interpolate.")
    }
    lapply(target_probs, function(p) quantile_linear_cdf_raster(value_rasters, probs, p))
  }

  list(
    prior_value_rasters = reinterpolate(prior_value_rasters, prior_probs),
    lik_value_rasters = reinterpolate(lik_value_rasters, lik_probs),
    probs = target_probs
  )
}

# ---------------------------------------------------------------------------
# Size-adaptive fusion for normal/beta/gamma
# ---------------------------------------------------------------------------

#' Analytic Normal percentiles, raster-native (`mu + sigma * qnorm(p)`)
#'
#' Pure raster arithmetic - `qnorm(p)` is a single scalar per requested probability (unlike
#' `qbeta()`/`qgamma()`, whose quantile depends on the shape parameters too, not just `p`), so no
#' per-cell dispatch (`terra::app()`) is needed at all; exact, and effectively free.
#' @param mu,sigma Normal parameter SpatRasters.
#' @param posterior_probs Numeric probabilities (0-1).
#' @return Named list of percentile `SpatRaster`s (names `"P1"`, `"P50"`, ... matching `posterior_probs`).
#' @keywords internal
qnorm_percentiles_raster <- function(mu, sigma, posterior_probs) {
  z <- stats::qnorm(posterior_probs)
  stats::setNames(
    lapply(z, function(zi) mu + sigma * zi),
    paste0("P", round(posterior_probs * 100))
  )
}

#' Analytic Beta/Gamma percentiles, raster-native, via per-cell `terra::app()`
#'
#' Unlike the Normal case, `qbeta()`/`qgamma()`'s quantile is a genuinely nonlinear function of
#' both the target probability AND the shape parameters, so there is no scalar-arithmetic
#' shortcut - this dispatches per cell via `terra::app()`, but still computes every requested
#' probability for a cell in one `vapply()` call (R's `qbeta()`/`qgamma()` are already vectorized
#' over their shape-parameter arguments), not one `terra::app()` pass per probability.
#' @param param1,param2 The two shape-parameter SpatRasters (`alpha`/`beta` or `shape`/`rate`).
#' @param qfun `stats::qbeta` or `stats::qgamma` (or any `function(p, param1, param2)`).
#' @param posterior_probs Numeric probabilities (0-1).
#' @return Named list of percentile `SpatRaster`s.
#' @keywords internal
closed_form_percentiles_raster <- function(param1, param2, qfun, posterior_probs) {
  combined <- c(param1, param2)
  np <- length(posterior_probs)
  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1)
    p1 <- row_mat[, 1]; p2 <- row_mat[, 2]
    vapply(posterior_probs, function(p) qfun(p, p1, p2), numeric(length(p1)))
  }
  res <- terra::app(combined, fun)
  stats::setNames(
    lapply(seq_len(np), function(i) res[[i]]),
    paste0("P", round(posterior_probs * 100))
  )
}

#' Per-unique-mukey closed-form family fit from real Monte Carlo draws, broadcast to a raster
#'
#' Precomputes a family-native parameter fit (Normal `mu`/`sigma`, Beta `alpha`/`beta`, or Gamma
#' `shape`/`rate`) once per UNIQUE mukey directly from that mukey's real simulated draws, then
#' broadcasts the result onto every cell sharing that mukey via one vectorized categorical lookup
#' (`terra::subst()`) - not a per-cell computation, unlike `fuse_general_kde()`'s raw_draws branch
#' (which needs a per-cell `density()` call and so is meaningfully slower). Verified cheap
#': at 100,172 synthetic cells and 150
#' unique mukeys, the whole precompute+broadcast pipeline measured ~0.18s total - negligible next
#' to the existing family-specific fit costs it complements (beta's Newton-Raphson fit alone
#' measured ~20s at the same cell count), because cost here scales with mukey CARDINALITY, not
#' cell count.
#'
#' @param mukey_raster Categorical mukey SpatRaster, aligned to the target grid.
#' @param mukey_draws A `lookup_mukey_draws()` result - named list keyed by mukey (character), each
#'   element a numeric vector of real simulated draws for that mukey.
#' @param family One of "normal", "beta", "gamma", "lognormal". `"lognormal"` (added
#') fits `mu`/`sigma` directly from
#'   `log(draws)` - the same LOG-SPACE parameterization `fuse_lognormal_adaptive()`'s large-AOI
#'   branch already uses internally (`convert_normal_to_lognormal()`'s output shape), so the merge
#'   there is a drop-in replacement of the percentile-triplet-derived log-space fit, not a new
#'   shape. Non-positive draws are dropped before taking `log()` (a real lognormal draw is always
#'   positive; any non-positive value reflects upstream noise, not signal) - a mukey left with fewer
#'   than 2 positive draws falls back to `NA` like every other family here.
#' @param bounds Required for `family = "beta"` - draws are scaled to `[0,1]` before fitting, to
#'   match `fit_beta_mle_newton_raster()`'s own scaled-space convention (see that function's docs) -
#'   `qfun`/downstream percentile code already expects `alpha`/`beta` in that scaled space.
#' @return A named list of `SpatRaster`s (`mu`/`sigma`, `alpha`/`beta`, or `shape`/`rate`) - `NA`
#'   at any cell whose mukey has no `mukey_draws` entry. No fallback is applied here; the caller
#'   merges this against the percentile-based fit for such cells (matching every other raw_draws
#'   branch's per-cell degrade-to-percentile-reconstruction convention).
#' @keywords internal
mukey_draws_closed_form_fit_raster <- function(mukey_raster, mukey_draws, family, bounds = NULL) {
  param_names <- switch(family, normal = c("mu", "sigma"), beta = c("alpha", "beta"),
                         gamma = c("shape", "rate"), lognormal = c("mu", "sigma"))

  mukey_codes <- suppressWarnings(as.numeric(names(mukey_draws)))
  valid <- !is.na(mukey_codes) & vapply(mukey_draws, length, integer(1)) > 0
  mukey_codes <- mukey_codes[valid]
  draws_list <- mukey_draws[valid]

  mukey_raster_plain <- mukey_raster
  levels(mukey_raster_plain) <- NULL

  if (length(mukey_codes) == 0) {
    na_r <- mukey_raster_plain * NA_real_
    return(stats::setNames(lapply(param_names, function(p) na_r), param_names))
  }

  fits <- lapply(draws_list, function(draws) {
    switch(family,
      normal = list(mu = mean(draws), sigma = stats::sd(draws)),
      gamma = convert_moments_to_gamma(mean(draws), stats::var(draws)),
      lognormal = {
        pos <- draws[draws > 0]
        if (length(pos) < 2) list(mu = NA_real_, sigma = NA_real_)
        else list(mu = mean(log(pos)), sigma = stats::sd(log(pos)))
      },
      beta = {
        scaled <- pmin(pmax((draws - bounds[1]) / (bounds[2] - bounds[1]), 1e-6), 1 - 1e-6)
        fit <- tryCatch(fitdistrplus::fitdist(scaled, "beta"), error = function(e) NULL)
        if (is.null(fit)) list(alpha = NA_real_, beta = NA_real_)
        else list(alpha = unname(fit$estimate[1]), beta = unname(fit$estimate[2]))
      }
    )
  })

  stats::setNames(lapply(param_names, function(p) {
    to <- vapply(fits, function(f) f[[p]], numeric(1))
    # others = NA_real_ is required: terra::subst()'s default behavior for a mukey code with no
    # match in `from` is to leave the cell's ORIGINAL value unchanged (the raw mukey code number
    # itself), not NA - silently corrupting fuse_closed_form()'s merge_raw_fit() fallback (which
    # relies on is.na() to detect "no raw-draws coverage, use the percentile fit instead") into
    # treating a mukey code like 602 as a fitted mu/alpha/shape of 602. Without this, an uncovered cell's Normal-Normal update would see
    # an absurd prior (mu=602, sigma=602) instead of falling back to the percentile-based fit.
    terra::subst(mukey_raster_plain, from = mukey_codes, to = to, others = NA_real_)
  }), param_names)
}

#' Per-unique-mukey empirical percentiles from real Monte Carlo draws, broadcast to a raster
#'
#' Companion to `mukey_draws_closed_form_fit_raster()` for routes that consume percentile VALUES
#' directly rather than a family parameter fit - `fuse_metalog_adapter()`, whose metalog fit
#' is an exact linear interpolation through fixed percentile knots (`solve()` on a basis matrix), not
#' a moment/density fit. There is no "metalog fit to raw draws" analogous to
#' `mukey_draws_closed_form_fit_raster()`'s normal/beta/gamma/lognormal cases - the real extension
#' for a percentile-interpolation method is to interpolate through REAL empirical percentiles
#' (`stats::quantile()` on the actual draws) instead of the percentile-reconstruction values. Same
#' per-unique-mukey precompute + `terra::subst()` broadcast pattern as its sibling - cheap, scales
#' with mukey cardinality, not cell count.
#'
#' @param mukey_raster Categorical mukey SpatRaster, aligned to the target grid.
#' @param mukey_draws A `lookup_mukey_draws()` result - named list keyed by mukey (character), each
#'   element a numeric vector of real simulated draws for that mukey.
#' @param probs Numeric probabilities (0-1) at which to compute each mukey's empirical quantile -
#'   pass the same knot probabilities the caller's `fit_metalog_linear_raster()` call will use.
#' @return A named list of `SpatRaster`s, one per `probs` entry (named `paste0("P", round(probs *
#'   100))`, matching this file's other percentile-list conventions) - `NA` at any cell whose mukey
#'   has no `mukey_draws` entry. No fallback is applied here; the caller merges this against the
#'   percentile-reconstruction values for such cells, same convention as every other raw_draws branch.
#' @keywords internal
mukey_draws_percentiles_raster <- function(mukey_raster, mukey_draws, probs) {
  mukey_codes <- suppressWarnings(as.numeric(names(mukey_draws)))
  valid <- !is.na(mukey_codes) & vapply(mukey_draws, length, integer(1)) > 0
  mukey_codes <- mukey_codes[valid]
  draws_list <- mukey_draws[valid]

  mukey_raster_plain <- mukey_raster
  levels(mukey_raster_plain) <- NULL
  layer_names <- paste0("P", round(probs * 100))

  if (length(mukey_codes) == 0) {
    na_r <- mukey_raster_plain * NA_real_
    return(stats::setNames(lapply(layer_names, function(nm) na_r), layer_names))
  }

  quantile_mat <- vapply(draws_list, function(draws) as.numeric(stats::quantile(draws, probs, names = FALSE)),
                          numeric(length(probs)))
  # quantile_mat: probs (rows) x mukeys (cols) when length(probs) > 1; vapply drops to a plain
  # vector when length(probs) == 1, so re-matrix defensively either way.
  quantile_mat <- matrix(quantile_mat, nrow = length(probs), ncol = length(draws_list))

  stats::setNames(lapply(seq_along(probs), function(i) {
    to <- quantile_mat[i, ]
    # others = NA_real_ - see mukey_draws_closed_form_fit_raster()'s identical note; without it an
    # uncovered mukey's cell would silently take on its own mukey code as a "percentile value".
    terra::subst(mukey_raster_plain, from = mukey_codes, to = to, others = NA_real_)
  }), layer_names)
}

#' Closed-form same-family fusion path for `fuse_adaptive()`, with a
#' per-cell fallback to Normal-moment fusion (re-expressed back into the
#' requested family) where the naive same-family fusion is infeasible.
#' @param prior_value_rasters,lik_value_rasters Lists of percentile-value SpatRasters.
#' @param percentile_probs Numeric probabilities matching the value-raster list order.
#' @param family One of "normal", "beta", "gamma".
#' @param bounds Required for `family = "beta"`.
#' @param posterior_probs Numeric probabilities (0-1) at which to report posterior percentiles.
#'   `NULL` (default) resolves to [FUSE_POSTERIOR_DEFAULT_PROBS]. Computed analytically
#'   (`qnorm`/`qbeta`/`qgamma` at the fused family parameters) - exact regardless of how extreme
#'   the requested probability is, unlike `fuse_general_kde()`'s sample/grid-based route.
#' @param mukey_raster,mukey_draws Optional - opts into `prior_fusion_method = "raw_draws"` for
#'   this (closed-form) route: fits the
#'   PRIOR side's family parameters directly from each cell's mukey's real Monte Carlo draws (via
#'   `mukey_draws_closed_form_fit_raster()`) instead of `fit_normal_raster()`/
#'   `fit_beta_mle_newton_raster()`/`fit_gamma_mom_raster()`'s percentile-triplet formulas - cells
#'   whose mukey has no draws entry fall back to the percentile-based fit for that cell only, same
#'   convention as `fuse_general_kde()`'s raw_draws branch. The LIKELIHOOD (SOLUS) side is
#'   unaffected either way - no raw-draws equivalent exists for it. `NULL` (default) for either
#'   parameter preserves the original percentile-only behavior exactly.
#' @return `list(posterior=, route_detail=, n_fallback_cells=)`. `posterior` carries a new
#'   `percentiles` field (named list of `SpatRaster`s) alongside the existing family-native
#'   parameters.
#' @keywords internal
fuse_closed_form <- function(prior_value_rasters, lik_value_rasters, percentile_probs, family, bounds,
                              posterior_probs = NULL, mukey_raster = NULL, mukey_draws = NULL) {
  effective_posterior_probs <- if (is.null(posterior_probs)) FUSE_POSTERIOR_DEFAULT_PROBS else posterior_probs
  use_raw_draws <- !is.null(mukey_raster) && !is.null(mukey_draws)

  # Merges a raw-draws-based fit onto the percentile-based fit, cell-by-cell - a cell whose mukey
  # has no mukey_draws entry (raw_fit == NA there) keeps its percentile-based value instead.
  merge_raw_fit <- function(percentile_fit, raw_fit) {
    stats::setNames(lapply(names(percentile_fit), function(p) {
      terra::ifel(is.na(raw_fit[[p]]), percentile_fit[[p]], raw_fit[[p]])
    }), names(percentile_fit))
  }

  if (family == "normal") {
    idx <- percentile_index(percentile_probs)
    fit_prior <- fit_normal_raster(prior_value_rasters[[idx$lo_idx]], prior_value_rasters[[idx$p50_idx]], prior_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
    if (use_raw_draws) {
      fit_prior <- merge_raw_fit(fit_prior, mukey_draws_closed_form_fit_raster(mukey_raster, mukey_draws, "normal"))
    }
    fit_lik <- fit_normal_raster(lik_value_rasters[[idx$lo_idx]], lik_value_rasters[[idx$p50_idx]], lik_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
    posterior <- fuse_normal_normal(fit_prior$mu, fit_prior$sigma, fit_lik$mu, fit_lik$sigma)
    posterior$percentiles <- qnorm_percentiles_raster(posterior$mu, posterior$sigma, effective_posterior_probs)
    return(list(posterior = posterior, route_detail = NULL, n_fallback_cells = 0))
  }

  if (family == "beta") {
    fit_prior <- fit_beta_mle_newton_raster(prior_value_rasters, bounds)
    if (use_raw_draws) {
      fit_prior <- merge_raw_fit(fit_prior, mukey_draws_closed_form_fit_raster(mukey_raster, mukey_draws, "beta", bounds = bounds))
    }
    fit_lik <- fit_beta_mle_newton_raster(lik_value_rasters, bounds)
    fused <- fuse_beta(fit_prior$alpha, fit_prior$beta, fit_lik$alpha, fit_lik$beta)

    prior_m <- convert_beta_to_moments(fit_prior$alpha, fit_prior$beta)
    lik_m <- convert_beta_to_moments(fit_lik$alpha, fit_lik$beta)
    span <- bounds[2] - bounds[1]
    fallback_normal <- fuse_normal_normal(
      bounds[1] + prior_m$mean * span, sqrt(prior_m$var) * span,
      bounds[1] + lik_m$mean * span, sqrt(lik_m$var) * span
    )
    fallback_mean <- (fallback_normal$mu - bounds[1]) / span
    fallback_var <- (fallback_normal$sigma / span)^2
    fallback_beta <- convert_moments_to_beta(fallback_mean, fallback_var)

    alpha_final <- terra::ifel(fused$feasible, fused$alpha, fallback_beta$alpha)
    beta_final <- terra::ifel(fused$feasible, fused$beta, fallback_beta$beta)
    route_detail <- terra::ifel(fused$feasible, 1, 0)
    n_fallback <- sum(!terra::values(fused$feasible), na.rm = TRUE)

    # alpha_final/beta_final fit Beta(alpha,beta) on the BOUNDS-SCALED [0,1] space (see
    # fit_beta_mle_newton_raster()) - qbeta() output must be rescaled back to raw units.
    scaled_percentiles <- closed_form_percentiles_raster(alpha_final, beta_final, stats::qbeta, effective_posterior_probs)
    percentiles <- lapply(scaled_percentiles, function(r) bounds[1] + r * span)

    posterior <- list(alpha = alpha_final, beta = beta_final, percentiles = percentiles)
    return(list(posterior = posterior, route_detail = route_detail, n_fallback_cells = n_fallback))
  }

  if (family == "gamma") {
    fit_prior <- fit_gamma_mom_raster(prior_value_rasters)
    if (use_raw_draws) {
      fit_prior <- merge_raw_fit(fit_prior, mukey_draws_closed_form_fit_raster(mukey_raster, mukey_draws, "gamma"))
    }
    fit_lik <- fit_gamma_mom_raster(lik_value_rasters)
    fused <- fuse_gamma(fit_prior$shape, fit_prior$rate, fit_lik$shape, fit_lik$rate)

    prior_m <- convert_gamma_to_moments(fit_prior$shape, fit_prior$rate)
    lik_m <- convert_gamma_to_moments(fit_lik$shape, fit_lik$rate)
    fallback_normal <- fuse_normal_normal(prior_m$mean, sqrt(prior_m$var), lik_m$mean, sqrt(lik_m$var))
    fallback_gamma <- convert_moments_to_gamma(fallback_normal$mu, fallback_normal$sigma^2)

    shape_final <- terra::ifel(fused$feasible, fused$shape, fallback_gamma$shape)
    rate_final <- terra::ifel(fused$feasible, fused$rate, fallback_gamma$rate)
    route_detail <- terra::ifel(fused$feasible, 1, 0)
    n_fallback <- sum(!terra::values(fused$feasible), na.rm = TRUE)

    percentiles <- closed_form_percentiles_raster(shape_final, rate_final, stats::qgamma, effective_posterior_probs)
    posterior <- list(shape = shape_final, rate = rate_final, percentiles = percentiles)
    return(list(posterior = posterior, route_detail = route_detail, n_fallback_cells = n_fallback))
  }
}

#' Default `update_prior()` grid resolution used by `fuse_general_kde()` (per-cell KDE fusion
#' route), coarser than `update_prior()`'s own standalone default of `0.01`.
#'
#' @section Why 0.1, not 0.01:
#' `update_prior()`'s `stats::density()` calls (kernel density estimation over
#' `seq(grid_min, grid_max, by = grid_resolution)`) are the dominant cost of
#' `fuse_general_kde()`'s per-cell loop - `Rprof()` profiling on a synthetic 10,000-cell raster
#' attributed 65% of total wall-clock time to `density()`
#' (`dnorm`/`fft` internals) at the `0.01` default, vs. 12% for the per-cell percentile-sampling
#' step. At `0.01`, a typical soil-property range (e.g. 20-50%) produces a ~3,000-point evaluation
#' grid per cell, per side. Coarsening to `0.1` cuts that grid ~10x and measured **~2.4x** faster
#' wall-clock on a real `fuse_adaptive()` call (23.86s -> 9.91s at 2,024 cells), while an accuracy
#' sweep across four representative percentile scenarios (narrow/wide/skewed distributions) found
#' the posterior mean/variance error at `0.1` stays under 0.02%/0.12% respectively vs. the `0.01`
#' reference - far inside the ~5-10% variability typical of field-measured soil properties (the
#' same bar the analogous cubic-spline-to-linear-interpolation tradeoff used elsewhere in this
#' project's Python sibling report was validated against). Going coarser (`0.25`+) is where error
#' starts compounding fast (variance error reaches double digits to 470% by `2.0`) - `0.1` is the
#' sweet spot, not an arbitrary round number.
#'
#' `update_prior()`'s own standalone default is deliberately left at `0.01` - this constant
#' only overrides the grid resolution `fuse_general_kde()` requests, so any other direct caller of
#' `update_prior()` is unaffected.
#' @keywords internal
FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION <- 0.1

#' Default `winsorize_probs` `fuse_general_kde()` passes to `update_prior()` on its `raw_draws`
#' branch only.
#'
#' @section Why this exists:
#' `update_prior()`'s `stats::density(bw = "nrd0")` bandwidth choice is sensitive to outliers a
#' raw-draws resample can genuinely contain - confirmed via a dedicated adversarial test: a
#' right-skewed synthetic prior's raw resample (realized range ~4-195) produced a measurably wider
#' KDE bandwidth, and a *less* accurate fused posterior mean on average across 8 seeds, than the
#' percentile-reconstruction route's resample (confined to ~7-64, the same data's P5-P95 range).
#' Clipping each side to its own P01-P99 range before density estimation directly counters that
#' outlier-driven inflation while retaining the bulk of the real distribution's shape - only the
#' most extreme 1% on each tail is affected, which is exactly the region a KDE-based fused mean is
#' most vulnerable to.
#'
#' `update_prior()`'s own standalone default is deliberately left `NULL` (unclipped) - this
#' constant only affects `fuse_general_kde()`'s `raw_draws` branch, which is the one path
#' confirmed to need it; the default percentile-reconstruction branch is already structurally
#' bounded near its own P5-P95 input range and is left untouched.
#' @keywords internal
FUSE_GENERAL_KDE_RAW_DRAWS_WINSORIZE_PROBS <- c(0.01, 0.99)

#' Default posterior percentile set computed by `fuse_general_kde()` (and, as later routes adopt
#' it, every other fusion route).
#'
#' @section Why 9 points, not the 5-point SSURGO/SOLUS input convention:
#' The 5-point set (`c(0.05, 0.25, 0.5, 0.75, 0.95)`) is what SSURGO/SOLUS happen to provide as raw
#' *input* percentile data - a real external constraint on the prior/likelihood side. The
#' *posterior* has no such constraint, since this package computes it, and 5 points under-resolves
#' exactly the tail behavior risk/cost-assessment use cases need (exceedance probabilities,
#' asymmetric credible intervals). This set covers a 98% interval (P01/P99, standard for tail-risk
#' work), a 90% interval (P05/P95), and an 80% interval (P10/P90), alongside the existing
#' median/IQR.
#'
#' @section Reliability note:
#' For the closed-form analytic routes (`qnorm`/`qbeta`/`qgamma`), extra percentiles are free
#' and exact regardless of how extreme they are. For
#' `fuse_general_kde()`'s sample/grid-based route, percentiles are read directly off
#' `update_prior()`'s discretized `posterior_prob` (see that function's own "Grid-based
#' percentiles" section) - exact relative to `grid_resolution` and NOT resampling-noise-limited,
#' but still subject to how well the underlying KDE itself approximates the tails from a finite
#' `n_samples`-sized input resample; the extreme percentiles here (P01/P99) are inherently noisier
#' than the inner ones (P10-P90) for that reason, not because of resampling.
#'
#' A caller needing even more extreme tail resolution (e.g. P0.1/P99.9) can pass its own
#' `posterior_probs` instead of relying on this default.
#' @keywords internal
FUSE_POSTERIOR_DEFAULT_PROBS <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)

#' Fully general fusion path for `fuse_adaptive()`: per cell, draw samples
#' from both sides' percentiles (via `simulate_from_percentiles(method=
#' "linear_cdf")`), fuse via `update_prior()`, and moment-match the
#' posterior samples back into the requested family so the output contract
#' matches the closed-form route.
#' @inheritParams fuse_closed_form
#' @param n_samples Samples drawn per side, per cell.
#' @param grid_resolution Passed to `update_prior()`. `NULL` (the caller-facing default from
#'   `fuse_adaptive()`) resolves to [FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION], not
#'   `update_prior()`'s own standalone default of `0.01` - see that constant's docs for why.
#' @param mukey_raster,mukey_draws Optional - opts into `prior_fusion_method = "raw_draws"` (see
#'). `mukey_raster` must already be
#'   aligned to the same grid as `prior_value_rasters`/`lik_value_rasters` (nearest-neighbor
#'   resampled, since it's categorical - bilinear would fabricate nonsensical mukey codes).
#'   `mukey_draws` is a `lookup_mukey_draws()` result. When both are supplied, the PRIOR side's
#'   samples for a cell are drawn from that cell's mukey's real Monte Carlo draws (computed
#'   once per unique mukey, not once per cell) instead of `sim_linear_cdf_batch()`'s
#'   percentile-reconstructed approximation - a genuine fidelity improvement, not just a speed one,
#'   since the real draws carry information the 5-percentile summary lost. The LIKELIHOOD (SOLUS)
#'   side is unaffected either way, since it has no raw-draws equivalent (see that function's own
#'   docs). `NULL` (default) for either parameter preserves the original percentile-reconstruction
#'   behavior for the prior side exactly - bit-identical output to before this parameter existed.
#' @param posterior_probs Numeric vector of probabilities (0-1) at which to report posterior
#'   percentiles, computed via `update_prior()`'s grid-based exact percentile machinery (see
#'   that function's docs). `NULL` (default) resolves to [FUSE_POSTERIOR_DEFAULT_PROBS] - unlike
#'   `update_prior()`'s own `NULL` default (which means "no percentiles, preserve the original
#'   plain-vector return"), this is an internal orchestration parameter with no external
#'   plain-vector consumers, so its `NULL` means "use the standard rich default" instead. Also used
#'   to compute `mu`/`sigma` (via the exact grid-based mean/var, not `mean()`/`var()` on a resampled
#'   vector - a strictly more accurate replacement for the previous computation, at no extra cost).
#' @keywords internal
fuse_general_kde <- function(prior_value_rasters, lik_value_rasters, percentile_probs,
                              family, bounds, n_samples, grid_resolution,
                              mukey_raster = NULL, mukey_draws = NULL, posterior_probs = NULL) {
  percentile_cols <- paste0("P", round(percentile_probs * 100))
  k <- length(percentile_probs)
  prior_stack <- terra::rast(prior_value_rasters); names(prior_stack) <- paste0("prior_", percentile_cols)
  lik_stack <- terra::rast(lik_value_rasters); names(lik_stack) <- paste0("lik_", percentile_cols)
  combined <- c(prior_stack, lik_stack)

  effective_grid_resolution <- if (is.null(grid_resolution)) FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION else grid_resolution
  effective_posterior_probs <- if (is.null(posterior_probs)) FUSE_POSTERIOR_DEFAULT_PROBS else posterior_probs
  np <- length(effective_posterior_probs)
  posterior_col_names <- paste0("P", round(effective_posterior_probs * 100))

  use_raw_draws <- !is.null(mukey_raster) && !is.null(mukey_draws)
  if (use_raw_draws) {
    # Precompute, once per unique mukey (not once per cell), a fixed n_samples-sized resample of
    # that mukey's real Monte Carlo draws. Reused via a per-cell mukey lookup below, rather
    # than reconstructing an approximation from percentiles for every cell that happens to share
    # the same mukey.
    mukey_prior_samples <- lapply(mukey_draws, function(vals) {
      if (length(vals) == 0) return(NULL)
      sample(vals, n_samples, replace = TRUE)
    })
    # Strip factor status before stacking - terra::c() silently coerces a categorical layer to
    # numeric anyway (with a warning), and the raw numeric codes are exactly the mukey values this
    # function needs to look up in mukey_prior_samples/mukey_draws; dropping the levels explicitly
    # gets the same numeric codes without the warning.
    mukey_raster_plain <- mukey_raster
    levels(mukey_raster_plain) <- NULL
    combined <- c(combined, mukey_raster_plain)
  }

  # terra::app() calls fun with a plain matrix (cells x layers) for a
  # chunk, and makes an initial probe call with a bare vector (not a
  # matrix) to infer the output layer count - fun must normalize both shapes.
  #
  # Cells with too many NA percentile
  # values (e.g. at raster edges after terra::resample()/align_percentile_probs() alignment
  # between two differently-extented grids) make extract_percentile_pairs() hard-stop() with
  # "Need at least 2 valid quantile columns." terra::app() is not per-cell-error-tolerant by
  # default, so one such cell would abort the entire raster operation - wrapped per-cell so a bad
  # cell degrades to NA output instead of crashing the whole fusion.
  #
  # Fast path (below) batches the percentile-sampling step (sim_linear_cdf_batch()) across every
  # cell in the chunk at once for cells with no missing/-1-sentinel percentile values (the common
  # case) - Rprof() profiling attributed ~9% of fuse_general_kde()'s total wall-clock to
  # extract_percentile_pairs()'s per-cell dispatch overhead alone (
  # Tier 4). update_prior()'s density() call still runs per cell (no vectorized form exists in
  # base R), and any cell with a missing value falls back to the original per-cell
  # simulate_from_percentiles() path unchanged, preserving the exact NA-degradation contract above.
  #
  # The raw_draws path (mukey_prior_samples lookup) is deliberately NOT merged into the
  # fast/slow-path vectorized split above - it's a separate, simpler row-by-row branch so the
  # already-tuned default path (still the only path any existing caller reaches, since
  # mukey_raster/mukey_draws default to NULL) is untouched. `update_prior()`'s density() call
  # still runs per cell either way (no vectorized form exists in base R) - but, the
  # LIKELIHOOD-side sampling (sim_linear_cdf_batch()) is now batched once across every
  # lik_no_na cell in the chunk up front, the same technique the fast path above already proves
  # out, instead of one row-by-row call per cell inside the loop
  #. The PRIOR-side fallback sampling (for a
  # cell whose mukey has no draws coverage) stays row-by-row - it's the rarer path and out of this
  # task's scope, though the same batching technique would apply there too.
  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1)
    n_cells <- nrow(row_mat)
    prior_mat <- row_mat[, 1:k, drop = FALSE]
    lik_mat <- row_mat[, (k + 1):(2 * k), drop = FALSE]
    mukey_col <- if (use_raw_draws) row_mat[, 2 * k + 1] else NULL

    # Columns: mean, var, then one column per `effective_posterior_probs` entry (in that order).
    result <- matrix(NA_real_, nrow = n_cells, ncol = 2 + np)

    if (use_raw_draws) {
      # Row-by-row: look up each cell's mukey's precomputed prior samples; cells with no mukey
      # coverage (NA) or an unknown/empty mukey fall back to the percentile-reconstruction route
      # for that cell only, preserving the existing degrade-to-something-reasonable behavior.
      lik_no_na <- stats::complete.cases(lik_mat) & rowSums(lik_mat == -1, na.rm = TRUE) == 0
      # Batch the likelihood-side sampling once across every lik_no_na cell in the chunk,
      # instead of one sim_linear_cdf_batch() call per cell inside the loop below - lik_row_pos[i]
      # is cell i's row within lik_samples_mat (only meaningful where lik_no_na[i] is TRUE).
      lik_samples_mat <- NULL
      lik_row_pos <- cumsum(lik_no_na)
      if (any(lik_no_na)) {
        lik_samples_mat <- tryCatch(
          sim_linear_cdf_batch(percentile_probs, lik_mat[which(lik_no_na), , drop = FALSE], n_samples),
          error = function(e) NULL
        )
      }
      for (i in seq_len(n_cells)) {
        mk <- mukey_col[i]
        prior_samples <- if (!is.na(mk)) mukey_prior_samples[[as.character(mk)]] else NULL
        if (is.null(prior_samples)) {
          # Fallback: no real draws for this cell's mukey - reconstruct from its own percentiles,
          # same as the default (non-raw_draws) path would.
          row <- prior_mat[i, ]
          if (anyNA(row) || any(row == -1, na.rm = TRUE)) next
          prior_samples <- tryCatch(
            sim_linear_cdf_batch(percentile_probs, matrix(row, nrow = 1), n_samples)[1, ],
            error = function(e) NULL
          )
          if (is.null(prior_samples)) next
        }
        if (!lik_no_na[i] || is.null(lik_samples_mat)) next
        lik_samples <- lik_samples_mat[lik_row_pos[i], ]
        post <- tryCatch(
          update_prior(prior_samples, lik_samples, grid_resolution = effective_grid_resolution,
                           winsorize_probs = FUSE_GENERAL_KDE_RAW_DRAWS_WINSORIZE_PROBS,
                           posterior_probs = effective_posterior_probs),
          error = function(e) NULL
        )
        if (!is.null(post)) result[i, ] <- c(post$mean, post$var, post$percentiles)
      }
      return(result)
    }

    no_na <- stats::complete.cases(prior_mat) & stats::complete.cases(lik_mat)
    no_sentinel <- rowSums(prior_mat == -1, na.rm = TRUE) == 0 & rowSums(lik_mat == -1, na.rm = TRUE) == 0
    fast_path <- no_na & no_sentinel

    if (any(fast_path)) {
      fast_idx <- which(fast_path)
      prior_samples_mat <- sim_linear_cdf_batch(percentile_probs, prior_mat[fast_idx, , drop = FALSE], n_samples)
      lik_samples_mat <- sim_linear_cdf_batch(percentile_probs, lik_mat[fast_idx, , drop = FALSE], n_samples)
      for (j in seq_along(fast_idx)) {
        post <- tryCatch(
          update_prior(prior_samples_mat[j, ], lik_samples_mat[j, ], grid_resolution = effective_grid_resolution,
                           posterior_probs = effective_posterior_probs),
          error = function(e) NULL
        )
        if (!is.null(post)) result[fast_idx[j], ] <- c(post$mean, post$var, post$percentiles)
      }
    }

    if (any(!fast_path)) {
      for (i in which(!fast_path)) {
        row <- row_mat[i, ]
        res <- tryCatch({
          prior_row <- as.data.frame(stats::setNames(as.list(row[1:k]), percentile_cols))
          lik_row <- as.data.frame(stats::setNames(as.list(row[(k + 1):(2 * k)]), percentile_cols))
          prior_samples <- simulate_from_percentiles(prior_row, method = "linear_cdf", percentile_cols = percentile_cols, n = n_samples)
          lik_samples <- simulate_from_percentiles(lik_row, method = "linear_cdf", percentile_cols = percentile_cols, n = n_samples)
          post <- update_prior(prior_samples, lik_samples, grid_resolution = effective_grid_resolution,
                                   posterior_probs = effective_posterior_probs)
          c(post$mean, post$var, post$percentiles)
        }, error = function(e) rep(NA_real_, 2 + np))
        result[i, ] <- res
      }
    }
    result
  }
  moments_r <- terra::app(combined, fun)
  mean_r <- moments_r[[1]]; var_r <- moments_r[[2]]
  percentile_rs <- stats::setNames(
    lapply(seq_len(np), function(i) moments_r[[2 + i]]),
    posterior_col_names
  )

  posterior <- switch(family,
    normal = list(mu = mean_r, sigma = sqrt(var_r)),
    beta = { span <- bounds[2] - bounds[1]; convert_moments_to_beta((mean_r - bounds[1]) / span, var_r / span^2) },
    gamma = convert_moments_to_gamma(mean_r, var_r)
  )
  posterior$percentiles <- percentile_rs
  list(posterior = posterior, route_detail = NULL, n_fallback_cells = 0)
}

#' Fuse a prior and likelihood belief distribution, choosing the fusion
#' route by AOI size rather than requiring the caller to pick.
#'
#' @param prior_value_rasters,lik_value_rasters Named lists of percentile-value
#'   `SpatRaster`s (e.g. `list(P0=, P5=, P50=, P95=, P100=)`), same layer
#'   count/order on both sides.
#' @param percentile_probs Numeric probabilities matching the value-raster
#'   list order (e.g. `c(0, 0.05, 0.5, 0.95, 1)`).
#' @param family One of "normal", "beta", "gamma" - the family both sides
#'   are fit in for the closed-form route, and the family the general
#'   route's output is moment-matched back into for a consistent output
#'   contract across routes.
#' @param bounds Required for `family = "beta"` (`c(lower, upper)`).
#' @param threshold_cells AOI cell count at or below which the general
#'   `update_prior()` route is used; above it, the closed-form route.
#' @param n_samples Samples drawn per side for the general route's KDE fusion.
#' @param grid_resolution Passed to `update_prior()` for the general route. `NULL` (default)
#'   resolves to [FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION] (`0.1`) - see that constant's docs for
#'   the profiling/accuracy justification for why this differs from `update_prior()`'s own
#'   standalone default of `0.01`.
#' @param verbose If TRUE (default), print the chosen route and why.
#' @param mukey_raster,mukey_draws Optional - opts into `prior_fusion_method = "raw_draws"` on
#'   whichever route runs (`fuse_general_kde()` for the general route -;
#'   `fuse_closed_form()` for the closed-form route -, added at negligible cost since
#'   that route's raw_draws fit is a per-mukey precompute + `terra::subst()` broadcast, not a
#'   per-cell computation - see `mukey_draws_closed_form_fit_raster()`'s docs). `NULL` (default)
#'   preserves original behavior exactly on both routes.
#' @param posterior_probs Passed through to whichever route runs (`fuse_general_kde()` for the
#'   general route, `fuse_closed_form()` for the closed-form route - see [FUSE_POSTERIOR_DEFAULT_PROBS]).
#'   `NULL` (default) resolves to the same rich default on both routes. The general route computes
#'   percentiles from `update_prior()`'s discretized posterior (exact relative to
#'   `grid_resolution`, but tail percentiles are limited by how well the KDE approximates the tails
#'   from a finite sample); the closed-form route computes them analytically (`qnorm`/`qbeta`/
#'   `qgamma` at the fused family parameters - exact regardless of how extreme the probability).
#' @return A list:
#'   - `posterior`: family-native parameter list (`mu`/`sigma` for "normal",
#'     `alpha`/`beta` for "beta", `shape`/`rate` for "gamma") - always this
#'     shape regardless of which route ran. Also carries `percentiles`: a named list of
#'     `SpatRaster`s (`"P01"`..`"P99"` by default), one per `posterior_probs` entry, on both
#'     routes.
#'   - `route`: one of `"bayesian_update_general"`, `"closed_form_normal"`,
#'     `"closed_form_beta"`, `"closed_form_gamma"`.
#'   - `route_detail`: for `family %in% c("beta","gamma")` on the
#'     closed-form route, a `SpatRaster` of 1 (native same-family fusion) /
#'     0 (per-cell fallback). `NULL` for "normal" and for the general route.
#'   - `n_fallback_cells`: count of cells that used the fallback.
#'   - `diagnostics`: `list(ncell=, threshold_cells=, family=, elapsed_sec=)`.
#' @family raster-fusion
#' @export
fuse_adaptive <- function(prior_value_rasters, lik_value_rasters, percentile_probs,
                           family = c("normal", "beta", "gamma"), bounds = NULL,
                           threshold_cells = 80000, n_samples = 500,
                           grid_resolution = NULL, verbose = TRUE,
                           mukey_raster = NULL, mukey_draws = NULL, posterior_probs = NULL) {
  family <- match.arg(family)
  if (family == "beta" && is.null(bounds)) stop("bounds = c(lower, upper) is required for family = 'beta'.")

  ncell <- terra::ncell(prior_value_rasters[[1]])
  use_general <- ncell <= threshold_cells
  route <- if (use_general) "bayesian_update_general" else paste0("closed_form_", family)

  if (verbose) {
    # %.0f / format(): `ncell` and `threshold_cells` are doubles, and `fuse_lognormal_adaptive()`
    # passes threshold_cells = Inf to force the general route - sprintf("%d", Inf) errors
    # ("invalid format '%d' ... use %f/%e/%g/%a"), which blocked every dist = "auto"/lognormal
    # run_fusion() call whenever verbose = TRUE.
    cat(sprintf(
      "[fuse_adaptive] AOI: %.0f cells (threshold: %s) -> route: %s\n",
      ncell, format(threshold_cells), route
    ))
  }

  t0 <- proc.time()[["elapsed"]]
  result <- if (use_general) {
    fuse_general_kde(prior_value_rasters, lik_value_rasters, percentile_probs,
                      family, bounds, n_samples, grid_resolution,
                      mukey_raster = mukey_raster, mukey_draws = mukey_draws,
                      posterior_probs = posterior_probs)
  } else {
    fuse_closed_form(prior_value_rasters, lik_value_rasters, percentile_probs, family, bounds,
                      posterior_probs = posterior_probs,
                      mukey_raster = mukey_raster, mukey_draws = mukey_draws)
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  if (verbose && !is.null(result$n_fallback_cells) && result$n_fallback_cells > 0) {
    cat(sprintf(
      "[fuse_adaptive] %d/%d cells (%.1f%%) required fallback to Normal-moment fusion (infeasible %s parameters)\n",
      result$n_fallback_cells, ncell, 100 * result$n_fallback_cells / ncell, family
    ))
  }

  list(
    posterior = result$posterior,
    route = route,
    route_detail = result$route_detail,
    n_fallback_cells = result$n_fallback_cells %||% 0,
    diagnostics = list(ncell = ncell, threshold_cells = threshold_cells, family = family, elapsed_sec = elapsed)
  )
}

# ---------------------------------------------------------------------------
# Lognormal and metalog adapters (fuse_adaptive() has no native route for these)
# ---------------------------------------------------------------------------

#' Lognormal adapter for `fuse_property_adaptive()`: size-adaptive general/closed-form route on
#' log-transformed values.
#' @keywords internal
#' @param posterior_probs Numeric probabilities (0-1) at which to report posterior percentiles.
#'   `NULL` (default) resolves to [FUSE_POSTERIOR_DEFAULT_PROBS]. Small-AOI branch: forwarded to
#'   `fuse_adaptive()`'s general route unchanged (see that function's docs). Large-AOI closed-form
#'   branch: computed via `qlnorm()` on the LOG-SPACE fused parameters (`posterior_log$mu`/`sigma`,
#'   captured before the moment-match-back-to-raw-space step below) - more faithful to the assumed
#'   lognormal family than deriving percentiles from the raw-space Normal-moment approximation
#'   (`posterior$mu`/`sigma`) would be.
#' @param mukey_raster,mukey_draws Optional raw-draws inputs (see `fuse_adaptive()`'s docs). Small-AOI
#'   branch: forwarded to `fuse_adaptive()`'s general route unchanged. Large-AOI closed-form branch
#'   Raw-draws inputs: each mukey's real draws are fit
#'   directly in log-space (`mukey_draws_closed_form_fit_raster(..., family = "lognormal")`) and
#'   merged over the percentile-triplet-derived log-space fit wherever a mukey has draws coverage -
#'   the same merge pattern used for normal/beta/gamma on the closed-form route.
fuse_lognormal_adaptive <- function(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs,
                                     ncell, threshold_cells, n_samples = 500, grid_resolution = NULL, verbose = TRUE,
                                     posterior_probs = NULL, mukey_raster = NULL, mukey_draws = NULL) {
  if (ncell <= threshold_cells) {
    # update_prior()'s general route makes no distributional assumption
    # at all, so it needs no log-space detour - force fuse_adaptive() into
    # its general route (threshold_cells=Inf) directly on raw values;
    # family="normal" only selects the (mu,sigma) OUTPUT shape.
    r <- fuse_adaptive(prior_value_rasters, lik_value_rasters, prior_probs,
                        family = "normal", threshold_cells = Inf,
                        n_samples = n_samples, grid_resolution = grid_resolution, verbose = verbose,
                        posterior_probs = posterior_probs,
                        mukey_raster = mukey_raster, mukey_draws = mukey_draws)
    return(list(posterior = r$posterior, route = r$route, route_detail = r$route_detail, n_fallback_cells = r$n_fallback_cells))
  }

  if (verbose) {
    cat(sprintf("[fuse_property_adaptive] AOI: %.0f cells (threshold: %s) -> route: closed_form_lognormal\n", ncell, format(threshold_cells)))
  }
  effective_posterior_probs <- if (is.null(posterior_probs)) FUSE_POSTERIOR_DEFAULT_PROBS else posterior_probs
  idx <- percentile_index(prior_probs)
  fit_prior <- fit_normal_raster(prior_value_rasters[[idx$lo_idx]], prior_value_rasters[[idx$p50_idx]], prior_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
  fit_lik <- fit_normal_raster(lik_value_rasters[[idx$lo_idx]], lik_value_rasters[[idx$p50_idx]], lik_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
  prior_log <- convert_normal_to_lognormal(fit_prior$mu, fit_prior$sigma)
  lik_log <- convert_normal_to_lognormal(fit_lik$mu, fit_lik$sigma)
  # Raw-draws extension: fit each mukey's real draws' log() directly (a more faithful
  # log-space fit than the percentile-triplet-derived approximation above), broadcast via
  # terra::subst(), and merge over prior_log wherever a mukey has draws coverage - the exact
  # merge_raw_fit() pattern fuse_closed_form() already uses for normal/beta/gamma. The
  # LIKELIHOOD (SOLUS) side has no raw-draws equivalent, same as every other route.
  if (!is.null(mukey_raster) && !is.null(mukey_draws)) {
    raw_prior_log <- mukey_draws_closed_form_fit_raster(mukey_raster, mukey_draws, "lognormal")
    prior_log <- list(
      mu = terra::ifel(is.na(raw_prior_log$mu), prior_log$mu, raw_prior_log$mu),
      sigma = terra::ifel(is.na(raw_prior_log$sigma), prior_log$sigma, raw_prior_log$sigma)
    )
  }
  posterior_log <- fuse_normal_normal(prior_log$mu, prior_log$sigma, lik_log$mu, lik_log$sigma)
  posterior <- convert_lognormal_to_normal(posterior_log$mu, posterior_log$sigma)
  posterior$percentiles <- closed_form_percentiles_raster(posterior_log$mu, posterior_log$sigma, stats::qlnorm, effective_posterior_probs)
  list(posterior = posterior, route = "closed_form_lognormal", route_detail = NULL, n_fallback_cells = 0)
}

#' Compute a metalog fit's mean/sd via numerical quadrature over its own raw
#' quantile function (`mean = integral of Q(p) dp`, `var = integral of
#' Q(p)^2 dp - mean^2`).
#'
#' @section Known limitation:
#' This is new glue code with no prior validated version (unlike the rest of
#' this file's math, which was validated against `fitdistrplus`/
#' closed-form references) - spot-check against real data before trusting it
#' in production, per the original source bundle's own caveat.
#'
#' Uses `quantile_metalog_linear_raster()` (valid for any `p` in `(0,1)`)
#' rather than the `_with_fallback()` wrapper for the quadrature itself:
#' `terra::ifel()` isn't lazy, so if the fallback's
#' `quantile_linear_cdf_raster()` were evaluated at every quadrature grid
#' point, it would error whenever the grid extends outside `full_probs`'s
#' own range. Instead, the quadrature integrates the raw metalog quantile
#' function over the full range, and INFEASIBLE cells fall back to
#' `fit_normal_raster()`'s closed-form moments computed directly from the
#' same percentile triplet.
#' @param fit,infeasible_r,bounds,boundedness Same as
#'   `check_metalog_feasibility_raster()`'s arguments/output for this `fit`.
#' @param full_value_rasters,full_probs The FULL percentile set (including
#'   p=0/p=1 if available) - used only for the infeasible-cell Normal fallback.
#' @param p_grid Probability grid for the quadrature.
#' @keywords internal
metalog_moments_raster <- function(fit, infeasible_r, full_value_rasters, full_probs, bounds, boundedness,
                                    p_grid = seq(0.001, 0.999, by = 0.005)) {
  n <- length(p_grid)
  sum_q <- NULL; sum_q2 <- NULL
  for (p in p_grid) {
    q_r <- quantile_metalog_linear_raster(fit, p, bounds, boundedness)
    sum_q <- if (is.null(sum_q)) q_r else sum_q + q_r
    sum_q2 <- if (is.null(sum_q2)) q_r^2 else sum_q2 + q_r^2
  }
  mean_r <- sum_q / n
  var_r <- sum_q2 / n - mean_r^2

  idx <- percentile_index(full_probs)
  fallback_fit <- fit_normal_raster(full_value_rasters[[idx$lo_idx]], full_value_rasters[[idx$p50_idx]], full_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)

  list(
    mu = terra::ifel(infeasible_r, fallback_fit$mu, mean_r),
    sigma = terra::ifel(infeasible_r, fallback_fit$sigma, sqrt(var_r))
  )
}

#' Metalog adapter for `fuse_property_adaptive()`: closed-form route via metalog-fitted
#' Normal-equivalent moments, used regardless of AOI size.
#' @keywords internal
#' @param posterior_probs Numeric probabilities (0-1) at which to report posterior percentiles.
#'   `NULL` (default) resolves to [FUSE_POSTERIOR_DEFAULT_PROBS]. Computed via
#'   `qnorm_percentiles_raster()` on the final fused `(mu, sigma)` - this route already reduces
#'   both sides to Normal-equivalent moments (`metalog_moments_raster()`'s quadrature) before the
#'   final `fuse_normal_normal()` step, so there is no richer distributional shape to draw
#'   percentiles from beyond that Normal approximation; carries forward this route's existing
#'   "less-validated glue code" caveat (see `metalog_moments_raster()`'s own docs) to its
#'   percentile output too.
#' @param mukey_raster,mukey_draws Optional raw-draws inputs. Metalog has no moment/density fit analogous to
#'   `fuse_closed_form()`'s normal/beta/gamma/lognormal routes (its fit is an exact linear
#'   interpolation through fixed percentile knots) - the real extension is interpolating through
#'   each mukey's REAL empirical percentiles (`mukey_draws_percentiles_raster()`) instead of the
#'   percentile-reconstruction values, wherever a mukey has draws coverage. `NULL` (default, or
#'   absorbed via `...` if passed positionally-incompatible) preserves the original
#'   percentile-reconstruction behavior exactly.
#' @param ... Absorbed and ignored - lets this function sit behind `fuse_property_adaptive()`'s
#'   `...` passthrough (which also forwards `n_samples`/`grid_resolution` meant for other routes)
#'   without erroring on arguments this route doesn't use.
fuse_metalog_adapter <- function(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs, property_config,
                                  verbose = TRUE, posterior_probs = NULL, mukey_raster = NULL, mukey_draws = NULL, ...) {
  interior <- percentile_interior(prior_probs)
  bounds <- property_config$bounds
  boundedness <- property_config$boundedness %||% "b"

  if (verbose) {
    cat("[fuse_property_adaptive] -> route: closed_form_metalog_adapter (always used regardless of AOI size)\n")
  }

  effective_posterior_probs <- if (is.null(posterior_probs)) FUSE_POSTERIOR_DEFAULT_PROBS else posterior_probs

  prior_interior_rasters <- prior_value_rasters[interior$idx]
  if (!is.null(mukey_raster) && !is.null(mukey_draws)) {
    raw_prior_percentiles <- mukey_draws_percentiles_raster(mukey_raster, mukey_draws, interior$probs)
    prior_interior_rasters <- Map(function(default_r, raw_r) terra::ifel(is.na(raw_r), default_r, raw_r),
                                   prior_interior_rasters, raw_prior_percentiles)
  }

  fit_prior <- fit_metalog_linear_raster(prior_interior_rasters, interior$probs, bounds, boundedness)
  fit_lik <- fit_metalog_linear_raster(lik_value_rasters[interior$idx], interior$probs, bounds, boundedness)

  infeasible_prior <- check_metalog_feasibility_raster(fit_prior, bounds, boundedness)
  infeasible_lik <- check_metalog_feasibility_raster(fit_lik, bounds, boundedness)

  prior_m <- metalog_moments_raster(fit_prior, infeasible_prior, prior_value_rasters, prior_probs, bounds, boundedness)
  lik_m <- metalog_moments_raster(fit_lik, infeasible_lik, lik_value_rasters, lik_probs, bounds, boundedness)

  posterior <- fuse_normal_normal(prior_m$mu, prior_m$sigma, lik_m$mu, lik_m$sigma)
  posterior$percentiles <- qnorm_percentiles_raster(posterior$mu, posterior$sigma, effective_posterior_probs)
  n_fallback <- sum(terra::values(infeasible_prior) | terra::values(infeasible_lik), na.rm = TRUE)

  list(
    posterior = posterior,
    route = "closed_form_metalog_adapter",
    route_detail = list(prior_infeasible = infeasible_prior, lik_infeasible = infeasible_lik),
    n_fallback_cells = n_fallback
  )
}

# ---------------------------------------------------------------------------
# dist = "auto": data-driven refinement, computed once per AOI-run
# ---------------------------------------------------------------------------

#' Resolve a property's distribution family, refining `dist = "auto"` from
#' the AOI's own percentile skew rather than a fixed per-property label.
#'
#' Computes a skewness proxy directly from the AOI's own percentile triplet:
#' `((high-median)-(median-low))/(high-low)`, aggregated once across the AOI
#' (default `median`) - NOT per-cell, which would create map discontinuity
#' artifacts at family-boundary cells.
#'
#' Deliberately narrow: only ever resolves to "beta" (if `bounds` is
#' configured), "normal", or "lognormal" - never "metalog" or "gamma",
#' since distinguishing those from a 3-point proxy isn't reliable. NEVER
#' overrides an explicit (non-"auto") `dist` value.
#'
#' @param property_config List with `dist` and optionally `bounds`/`auto_skew_threshold`.
#' @param prior_value_rasters,prior_probs Prior percentile-value rasters/probabilities.
#' @param aggregate_fun Function to aggregate the AOI-wide skew proxy (default `stats::median`).
#' @return `list(dist=, dist_source= one of "config"/"auto", skew_proxy=
#'   NA_real_ unless dist_source=="auto")`.
#' @family raster-fusion
#' @export
resolve_property_dist <- function(property_config, prior_value_rasters, prior_probs, aggregate_fun = stats::median) {
  if (!identical(property_config$dist, "auto")) {
    return(list(dist = property_config$dist, dist_source = "config", skew_proxy = NA_real_))
  }

  idx <- percentile_index(prior_probs)
  agg_low <- aggregate_fun(terra::values(prior_value_rasters[[idx$lo_idx]]), na.rm = TRUE)
  agg_mid <- aggregate_fun(terra::values(prior_value_rasters[[idx$p50_idx]]), na.rm = TRUE)
  agg_high <- aggregate_fun(terra::values(prior_value_rasters[[idx$hi_idx]]), na.rm = TRUE)
  skew_proxy <- ((agg_high - agg_mid) - (agg_mid - agg_low)) / (agg_high - agg_low)

  threshold <- property_config$auto_skew_threshold %||% 0.15
  dist <- if (!is.null(property_config$bounds)) "beta" else if (abs(skew_proxy) < threshold) "normal" else "lognormal"
  list(dist = dist, dist_source = "auto", skew_proxy = skew_proxy)
}

# ---------------------------------------------------------------------------
# Top-level entry point (non-compositional properties)
# ---------------------------------------------------------------------------

#' Fuse one property's prior and likelihood, dispatching to the right
#' family/route automatically from `property_config`.
#'
#' The raster-native counterpart of `core-fusion.R`'s `fuse_property()`,
#' and this toolkit's top-level entry point for non-compositional properties
#' (see `fuse_texture_group()` for the compositional/texture case) - the
#' raster analogue of `run_fusion()` minus the SSURGO/SOLUS
#' fetch-and-cache wrapper (see this file's `@section Deliberately out of
#' scope:` above).
#'
#' @param prior_value_rasters,lik_value_rasters Named lists of percentile-value SpatRasters.
#' @param prior_probs,lik_probs Matching probability vectors (reconciled via
#'   `align_percentile_probs()` if they differ).
#' @param property_config List with `dist` (one of "normal"/"beta"/"gamma"/
#'   "lognormal"/"metalog"/"auto"), and optionally `bounds`/`boundedness`/
#'   `auto_skew_threshold`.
#' @param threshold_cells AOI cell count at or below which `fuse_adaptive()`
#'   uses its general route.
#' @param ... Additional arguments passed to `fuse_adaptive()`/`fuse_lognormal_adaptive()`/
#'   `fuse_metalog_adapter()` (e.g. `n_samples`, `grid_resolution`, `verbose`, `posterior_probs` -
#'   see [FUSE_POSTERIOR_DEFAULT_PROBS]; reaches every route - `fuse_metalog_adapter()`
#'   absorbs and ignores the arguments it doesn't use, e.g. `n_samples`). Also how
#'   `mukey_raster`/`mukey_draws` (opting into `prior_fusion_method = "raw_draws"`, see
#'   `fuse_general_kde()`'s docs) reach the underlying routes: both the general-KDE and
#'   closed-form routes support it for `dist %in% c("normal","beta","gamma")`, and it
#'   reaches `fuse_lognormal_adaptive()`'s small-AOI general branch AND large-AOI closed-form
#'   branch the same way. `fuse_metalog_adapter()` also has a raw-draws fit
#'   (empirical per-mukey percentiles feed the same interpolation machinery) - but
#'   `resolve_want_raw_draws()`'s default still excludes `dist = "metalog"`, mirroring the
#'   texture-group route's treatment: a fit existing is not the same as its cost/benefit being
#'   benchmarked, and that default decision is kept separate. Explicitly passing `mukey_raster`/`mukey_draws`
#'   for `dist = "metalog"` still works.
#' @return `list(posterior=, route=, route_detail=, n_fallback_cells=, dist=,
#'   dist_source=, skew_proxy=)`.
#' @keywords internal
fuse_property_adaptive <- function(prior_value_rasters, prior_probs,
                                    lik_value_rasters, lik_probs,
                                    property_config, threshold_cells = 80000, ...) {
  aligned <- align_percentile_probs(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs)
  prior_value_rasters <- aligned$prior_value_rasters
  lik_value_rasters <- aligned$lik_value_rasters
  prior_probs <- lik_probs <- aligned$probs

  resolved <- resolve_property_dist(property_config, prior_value_rasters, prior_probs)
  dist <- resolved$dist
  ncell <- terra::ncell(prior_value_rasters[[1]])

  result <- if (dist %in% c("normal", "beta", "gamma")) {
    r <- fuse_adaptive(prior_value_rasters, lik_value_rasters, prior_probs,
                        family = dist, bounds = property_config$bounds,
                        threshold_cells = threshold_cells, ...)
    list(posterior = r$posterior, route = r$route, route_detail = r$route_detail, n_fallback_cells = r$n_fallback_cells)
  } else if (dist == "lognormal") {
    fuse_lognormal_adaptive(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs, ncell, threshold_cells, ...)
  } else if (dist == "metalog") {
    fuse_metalog_adapter(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs, property_config, ...)
  } else {
    stop(sprintf("Unknown dist '%s'.", dist))
  }

  c(result, list(dist = dist, dist_source = resolved$dist_source, skew_proxy = resolved$skew_proxy))
}

# ---------------------------------------------------------------------------
# Compositional (texture) fusion
# ---------------------------------------------------------------------------

#' Look up a compositional group's member property ids, in configured order
#'
#' Reuses soilSIM's own `config$monte_carlo$composition_groups` convention
#' (see `R/core-distributions.R`'s `resolve_composition_groups()`) rather than
#' the original source bundle's separate `PROPERTIES`/`config.R` global
#' registry - order is trusted from the config, not hardcoded to literal
#' "clay"/"sand"/"silt" names, exactly like the already-ported
#' `fuse_texture_group_from_triplets()` (`core-fusion.R`) already
#' documents ("positional-role placeholders... not an identity requirement").
#'
#' @param group Composition group name (e.g. `"texture"`).
#' @param composition_groups A `config$monte_carlo$composition_groups`-shaped
#'   list, e.g. `list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))`.
#' @return Character vector of member property ids, in configured order.
#' @keywords internal
group_members <- function(group, composition_groups) {
  group_def <- composition_groups[[group]]
  if (is.null(group_def)) {
    stop(sprintf("group_members(): no composition group '%s' found.", group))
  }
  members <- group_def$members
  if (length(members) < 2) {
    stop(sprintf("group_members(): composition group '%s' needs at least 2 members.", group))
  }
  members
}

#' Vectorized per-chunk core of `fuse_texture_group()`
#'
#' Computes, for every cell (row) in `row_mat` at once, the same
#' prior/likelihood ILR moments + bivariate-normal fusion that
#' `estimate_ilr_moments_mc()` + `fuse_bivariate_normal()` compute one cell at
#' a time - see `fuse_texture_group()`'s PERF comment at its call site.
#'
#' RNG draws are generated as a single `rnorm()` call ordered so the
#' underlying standard-normal stream is consumed in exactly the same
#' cell-major, then prior-before-lik, then clay/sand/silt, then
#' draw-1..n_mc order the original per-cell loop consumed it in (rnorm()'s
#' `mean`/`sd` args only affine-transform each already-drawn standard normal;
#' they don't change how many values are drawn or in what order) - so this is
#' a pure vectorization, not a behavior change. `fuse_bivariate_normal()`'s
#' 2x2 `solve()`-based fusion is replaced with the closed-form 2x2 matrix
#' inverse formula (`solve()` isn't vectorizable across cells), which can
#' differ from `solve()` by floating-point rounding (~1e-10 relative) but is
#' otherwise the identical linear algebra.
#'
#' @param row_mat A matrix with one row per cell and the 18
#'   `<id>_<prior|lik>_<lo|p50|hi>` columns `fuse_texture_group()` builds (plus one trailing
#'   `mukey_code` column when `mukey_prior_texture_samples` is supplied).
#' @param clay_id,sand_id,silt_id The three members' `id`s, in ILR order.
#' @param prior_z,lik_z Standard-normal quantiles for the prior/likelihood
#'   low-high intervals (as in `estimate_ilr_moments_mc()`'s `z` argument).
#' @param n_mc Monte Carlo sample size per side (matches
#'   `estimate_ilr_moments_mc()`'s default).
#' @param max_cells_per_subchunk Upper bound on cells processed by one
#'   `fuse_texture_group_batch_core()` call - each cell needs `6 * n_mc`
#'   doubles of intermediate MC-draw storage (`n_mc=2000` -> ~96KB/cell), so
#'   `terra::app()`'s own chunk size (which only accounts for the raw,
#'   non-MC-expanded raster layers) can hand `fuse_texture_group()` chunks
#'   large enough to exhaust memory - confirmed empirically: a 50,000-cell
#'   chunk failed with "cannot allocate vector of size 4.5 Gb". Sub-chunking
#'   here, in cell order, doesn't change the RNG stream order (still
#'   cell-major) so results are unaffected.
#' @param mukey_prior_texture_samples Optional - opts into `prior_fusion_method = "raw_draws"`
#'. A named list, keyed by mukey (as
#'   character), each element an `n_mc x 3` matrix (columns `clay_total`/`sand_total`/
#'   `silt_total`) of real joint texture draws for that mukey - see
#'   `lookup_mukey_texture_draws()`. When supplied, `row_mat`'s last column must be the per-cell
#'   mukey code. `NULL` (default) preserves original behavior exactly.
#' @return A `nrow(row_mat) x 5` matrix, columns `mu1, mu2, S11, S12, S22`.
#' @keywords internal
fuse_texture_group_batch <- function(row_mat, clay_id, sand_id, silt_id, prior_z, lik_z, n_mc = 2000,
                                      max_cells_per_subchunk = 2000,
                                      mukey_prior_texture_samples = NULL) {
  ncell_total <- nrow(row_mat)
  if (ncell_total <= max_cells_per_subchunk) {
    return(fuse_texture_group_batch_core(row_mat, clay_id, sand_id, silt_id, prior_z, lik_z, n_mc,
                                          mukey_prior_texture_samples = mukey_prior_texture_samples))
  }
  out <- matrix(NA_real_, nrow = ncell_total, ncol = 5, dimnames = list(NULL, c("mu1", "mu2", "S11", "S12", "S22")))
  for (start in seq(1, ncell_total, by = max_cells_per_subchunk)) {
    end <- min(start + max_cells_per_subchunk - 1, ncell_total)
    out[start:end, ] <- fuse_texture_group_batch_core(
      row_mat[start:end, , drop = FALSE], clay_id, sand_id, silt_id, prior_z, lik_z, n_mc,
      mukey_prior_texture_samples = mukey_prior_texture_samples
    )
  }
  out
}

#' @rdname fuse_texture_group_batch
#' @keywords internal
fuse_texture_group_batch_core <- function(row_mat, clay_id, sand_id, silt_id, prior_z, lik_z, n_mc = 2000,
                                           mukey_prior_texture_samples = NULL) {
  use_raw_draws <- !is.null(mukey_prior_texture_samples)
  mukey_col <- if (use_raw_draws) row_mat[, ncol(row_mat)] else NULL
  if (use_raw_draws) row_mat <- row_mat[, -ncol(row_mat), drop = FALSE]

  ncell <- nrow(row_mat)
  col <- function(nm) row_mat[, nm]

  side_mean_sd <- function(side, z) {
    list(
      clay = list(mean = col(paste0(clay_id, "_", side, "_p50")),
                  sd = (col(paste0(clay_id, "_", side, "_hi")) - col(paste0(clay_id, "_", side, "_lo"))) / (2 * z)),
      sand = list(mean = col(paste0(sand_id, "_", side, "_p50")),
                  sd = (col(paste0(sand_id, "_", side, "_hi")) - col(paste0(sand_id, "_", side, "_lo"))) / (2 * z)),
      silt = list(mean = col(paste0(silt_id, "_", side, "_p50")),
                  sd = (col(paste0(silt_id, "_", side, "_hi")) - col(paste0(silt_id, "_", side, "_lo"))) / (2 * z))
    )
  }
  lik <- side_mean_sd("lik", lik_z)

  # Cell-major block: n_mc rows x ncell cols, column c filled with v[c] (rep(v, each=n_mc),
  # column-major-filled, puts v[c] in every row of column c without matrix()'s byrow=TRUE - which
  # internally transposes and is measurably slower for large matrices) - i.e. column c is v[c]
  # repeated n_mc times, matching one cell's n_mc consecutive draws. Kept in this n_mc(row) x
  # ncell(col) orientation throughout (never transposed) since R matrices are column-major, so a
  # cell's n_mc draws are already contiguous - transposing to ncell x n_mc just to use
  # rowMeans()/rowSums() (as an earlier version of this function did) re-shuffles the whole
  # allocation for no benefit; colMeans()/colSums() over this orientation do the same reduction
  # without it.
  block <- function(v) matrix(rep(v, each = n_mc), nrow = n_mc, ncol = ncell)

  if (use_raw_draws) {
    # PRIOR side: look up each cell's mukey in the precomputed (once-per-unique-mukey) real
    # joint-texture-draws sample instead of reconstructing independent per-fraction normals from
    # percentiles - preserves the real clay/sand/silt cross-correlation those normals discard.
    # Cells with no real-draws entry for their mukey (NA mukey, or a mukey missing from the
    # lookup) fall back to the percentile-reconstruction approach for that cell only.
    prior_clay_s <- matrix(NA_real_, nrow = n_mc, ncol = ncell)
    prior_sand_s <- matrix(NA_real_, nrow = n_mc, ncol = ncell)
    prior_silt_s <- matrix(NA_real_, nrow = n_mc, ncol = ncell)
    fallback_cells <- logical(ncell)
    for (j in seq_len(ncell)) {
      mk <- mukey_col[j]
      samp <- if (!is.na(mk)) mukey_prior_texture_samples[[as.character(mk)]] else NULL
      if (is.null(samp)) {
        fallback_cells[j] <- TRUE
        next
      }
      prior_clay_s[, j] <- samp[, "clay_total"]
      prior_sand_s[, j] <- samp[, "sand_total"]
      prior_silt_s[, j] <- samp[, "silt_total"]
    }
    if (any(fallback_cells)) {
      prior_fb <- side_mean_sd("prior", prior_z)
      for (j in which(fallback_cells)) {
        prior_clay_s[, j] <- stats::rnorm(n_mc, prior_fb$clay$mean[j], prior_fb$clay$sd[j])
        prior_sand_s[, j] <- stats::rnorm(n_mc, prior_fb$sand$mean[j], prior_fb$sand$sd[j])
        prior_silt_s[, j] <- stats::rnorm(n_mc, prior_fb$silt$mean[j], prior_fb$silt$sd[j])
      }
    }
    prior_clay_s <- pmax(prior_clay_s, 0.01)
    prior_sand_s <- pmax(prior_sand_s, 0.01)
    prior_silt_s <- pmax(prior_silt_s, 0.01)

    lik_ordered_props <- list(lik$clay, lik$sand, lik$silt)
    lik_mean_stack <- do.call(rbind, lapply(lik_ordered_props, function(p) block(p$mean)))
    lik_sd_stack <- do.call(rbind, lapply(lik_ordered_props, function(p) block(p$sd)))
    lik_draws <- stats::rnorm(length(lik_mean_stack), mean = as.vector(lik_mean_stack), sd = as.vector(lik_sd_stack))
    lik_draws_mat <- matrix(lik_draws, nrow = 3 * n_mc, ncol = ncell)
    lik_side_block <- function(i) pmax(lik_draws_mat[((i - 1) * n_mc + 1):(i * n_mc), , drop = FALSE], 0.01)

    side_block_list <- list(prior_clay_s, prior_sand_s, prior_silt_s,
                             lik_side_block(1), lik_side_block(2), lik_side_block(3))
    side_block <- function(i) side_block_list[[i]]
  } else {
    prior <- side_mean_sd("prior", prior_z)
    ordered_props <- list(prior$clay, prior$sand, prior$silt, lik$clay, lik$sand, lik$silt)
    mean_stack <- do.call(rbind, lapply(ordered_props, function(p) block(p$mean)))
    sd_stack <- do.call(rbind, lapply(ordered_props, function(p) block(p$sd)))

    draws <- stats::rnorm(length(mean_stack), mean = as.vector(mean_stack), sd = as.vector(sd_stack))
    draws_mat <- matrix(draws, nrow = 6 * n_mc, ncol = ncell)
    side_block <- function(i) pmax(draws_mat[((i - 1) * n_mc + 1):(i * n_mc), , drop = FALSE], 0.01)
  }

  moments <- function(clay_s, sand_s, silt_s) {
    total <- clay_s + sand_s + silt_s
    z_res <- ilr_forward(as.vector(clay_s / total * 100), as.vector(sand_s / total * 100), as.vector(silt_s / total * 100))
    z1 <- matrix(z_res[, "z1"], nrow = n_mc, ncol = ncell)
    z2 <- matrix(z_res[, "z2"], nrow = n_mc, ncol = ncell)
    mu1 <- colMeans(z1); mu2 <- colMeans(z2)
    d1 <- z1 - rep(mu1, each = n_mc); d2 <- z2 - rep(mu2, each = n_mc)
    denom <- n_mc - 1
    list(mu1 = mu1, mu2 = mu2,
         S11 = colSums(d1 * d1) / denom, S12 = colSums(d1 * d2) / denom, S22 = colSums(d2 * d2) / denom)
  }
  prior_m <- moments(side_block(1), side_block(2), side_block(3))
  lik_m <- moments(side_block(4), side_block(5), side_block(6))

  # Vectorized closed-form 2x2 fuse_bivariate_normal(): precision = inverse(covariance),
  # posterior precision = sum of precisions, posterior mean = posterior_Sigma %*% sum(precision
  # %*% mu). 2x2 inverse of [[a,b],[b,d]] is (1/(ad-b^2)) * [[d,-b],[-b,a]].
  det1 <- prior_m$S11 * prior_m$S22 - prior_m$S12^2
  P1_11 <- prior_m$S22 / det1; P1_12 <- -prior_m$S12 / det1; P1_22 <- prior_m$S11 / det1
  det2 <- lik_m$S11 * lik_m$S22 - lik_m$S12^2
  P2_11 <- lik_m$S22 / det2; P2_12 <- -lik_m$S12 / det2; P2_22 <- lik_m$S11 / det2

  Pp_11 <- P1_11 + P2_11; Pp_12 <- P1_12 + P2_12; Pp_22 <- P1_22 + P2_22
  detp <- Pp_11 * Pp_22 - Pp_12^2
  S11 <- Pp_22 / detp; S12 <- -Pp_12 / detp; S22 <- Pp_11 / detp

  T1 <- P1_11 * prior_m$mu1 + P1_12 * prior_m$mu2 + P2_11 * lik_m$mu1 + P2_12 * lik_m$mu2
  T2 <- P1_12 * prior_m$mu1 + P1_22 * prior_m$mu2 + P2_12 * lik_m$mu1 + P2_22 * lik_m$mu2
  mu1 <- S11 * T1 + S12 * T2
  mu2 <- S12 * T1 + S22 * T2

  cbind(mu1 = mu1, mu2 = mu2, S11 = S11, S12 = S12, S22 = S22)
}

#' Raster-native ILR-posterior percentile sampler for `fuse_texture_group()`
#'
#' Draws `n_mc` bivariate-Normal ILR-space samples per cell from the fused `(mu1, mu2, S11, S12,
#' S22)` posterior - the same distribution `R/core-distributions.R`'s `sample_ilr_posterior()` draws
#' from in the scalar case, here vectorized raster-wide via `terra::app()` - applies each cell's
#' own 2x2 Cholesky factor (closed-form: `u11 = sqrt(S11)`, `u12 = S12/u11`, `u22 = sqrt(S22 -
#' u12^2)`, matching base R's `chol()` upper-triangular convention exactly), `ilr_inverse()`s the
#' result to `(clay, sand, silt)` triples, then computes `quantile()` across draws per cell per
#' fraction. Unlike `update_prior()`'s grid-based routes, this posterior has no discretized-grid
#' shortcut - it's a genuine bivariate Monte Carlo sampler, so `n_mc` is the only lever on
#' percentile reliability here.
#'
#' @section Sum-to-100 caveat (expected, not a bug):
#' The point estimate (`posterior$value`, via `ilr_inverse()` on the fused mean) sums to 100 by
#' construction for every cell. Per-fraction PERCENTILES do **not** generally sum to 100 -
#' percentiles of correlated marginals don't add linearly (e.g. clay's P95 and sand's P95 can both
#' be simultaneously "high" draws from different tail regions of the joint distribution, so their
#' sum can exceed - or their P5s' sum can fall short of - 100). This is an expected mathematical
#' property of taking marginal percentiles of a compositional joint distribution, not a defect.
#' @param ilr_mu_r 2-layer SpatRaster (`mu1`, `mu2`) - `fuse_texture_group_batch()`'s `ilr_mu` output.
#' @param ilr_Sigma_r 3-layer SpatRaster (`S11`, `S12`, `S22`) - that function's `ilr_Sigma` output.
#' @param posterior_probs Numeric probabilities (0-1) at which to report posterior percentiles.
#' @param n_mc Monte Carlo draw count per cell. Default 2000 matches this file's existing
#'   `texture_n_mc` convention (the separate MC sample size `fuse_texture_group_batch_core()` uses
#'   to ESTIMATE the `(mu1,mu2,S11,S12,S22)` moments in the first place) - a coincidentally shared
#'   default, not the same computation; this parameter controls a distinct, later sampling step
#'   that draws FROM the already-fused posterior distribution.
#' @param max_cells_per_subchunk Same memory rationale as `fuse_texture_group_batch()`'s own
#'   parameter - each cell needs `2*n_mc` standard-normal draws plus a `3*n_mc`-row ILR-inverse
#'   composition matrix, so large `terra::app()` chunks can exhaust memory.
#' @return `list(clay=, sand=, silt=)`, each a named list of percentile `SpatRaster`s (names
#'   `"P1"`, `"P50"`, ... matching `posterior_probs`).
#' @keywords internal
texture_group_percentiles_raster <- function(ilr_mu_r, ilr_Sigma_r, posterior_probs, n_mc = 2000,
                                              max_cells_per_subchunk = 2000) {
  np <- length(posterior_probs)
  percentile_names <- paste0("P", round(posterior_probs * 100))
  combined <- c(ilr_mu_r, ilr_Sigma_r)

  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1)
    ncell_total <- nrow(row_mat)
    result <- matrix(NA_real_, nrow = ncell_total, ncol = 3 * np)

    for (start in seq(1, ncell_total, by = max_cells_per_subchunk)) {
      end <- min(start + max_cells_per_subchunk - 1, ncell_total)
      chunk <- row_mat[start:end, , drop = FALSE]
      ncell <- nrow(chunk)
      mu1 <- chunk[, 1]; mu2 <- chunk[, 2]
      S11 <- chunk[, 3]; S12 <- chunk[, 4]; S22 <- chunk[, 5]

      u11 <- sqrt(S11); u12 <- S12 / u11; u22 <- sqrt(pmax(S22 - u12^2, 0))

      e1 <- matrix(stats::rnorm(n_mc * ncell), nrow = n_mc, ncol = ncell)
      e2 <- matrix(stats::rnorm(n_mc * ncell), nrow = n_mc, ncol = ncell)
      z1 <- e1 * rep(u11, each = n_mc) + rep(mu1, each = n_mc)
      z2 <- e1 * rep(u12, each = n_mc) + e2 * rep(u22, each = n_mc) + rep(mu2, each = n_mc)

      comp <- ilr_inverse(as.vector(z1), as.vector(z2))
      clay_s <- matrix(comp[, "clay"], nrow = n_mc, ncol = ncell)
      sand_s <- matrix(comp[, "sand"], nrow = n_mc, ncol = ncell)
      silt_s <- matrix(comp[, "silt"], nrow = n_mc, ncol = ncell)

      # Cells outside real SSURGO/SOLUS coverage (a documented, common occurrence - e.g..8's
      # live-AOI validation found 401/598 NA cells for one real AOI) arrive here with NA
      # mu1/mu2/S11/S12/S22, which propagates to an all-NA draws column for that cell. quantile()'s
      # default na.rm = FALSE hard-errors on an all-NA input, which - unlike every other per-cell
      # fallback in this file (fuse_general_kde()'s tryCatch-wrapped loop, etc.) - would crash the
      # ENTIRE terra::app() call for the whole raster over one missing-coverage cell. Found live via
      #'s real-AOI validation (exactly the class of
      # bug synthetic unit tests, which never included an NA cell, couldn't catch). Degrade that one
      # cell to NA percentiles instead, matching this file's established NA-cell convention.
      col_quantile <- function(col) {
        if (!any(is.finite(col))) return(rep(NA_real_, np))
        stats::quantile(col, probs = posterior_probs, names = FALSE, na.rm = TRUE)
      }
      # apply(..., FUN = col_quantile) returns an np x ncell matrix (one column per cell) -
      # transpose each fraction's block into the ncell x np shape `result` expects.
      clay_q <- t(apply(clay_s, 2, col_quantile))
      sand_q <- t(apply(sand_s, 2, col_quantile))
      silt_q <- t(apply(silt_s, 2, col_quantile))
      result[start:end, ] <- cbind(clay_q, sand_q, silt_q)
    }
    result
  }

  res <- terra::app(combined, fun)
  fraction_percentiles <- function(offset) {
    stats::setNames(lapply(seq_len(np), function(i) res[[offset + i]]), percentile_names)
  }
  list(clay = fraction_percentiles(0), sand = fraction_percentiles(np), silt = fraction_percentiles(2 * np))
}

#' Fuse a compositional group's members JOINTLY via ILR fusion
#' (`R/core-distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`,
#' `R/core-fusion.R`'s `fuse_bivariate_normal()`), rather than
#' independently via `fuse_beta()` per member - independent fusion
#' measurably breaks sum-to-100 (up to 10.5 percentage points on realistic
#' synthetic data). The raster counterpart of the
#' already-ported scalar `fuse_texture_group_from_triplets()`.
#'
#' @param fetched A list (one entry per group member, in `group_members()`'s
#'   order) of `list(id=, prior=<named list of ALIGNED percentile-value
#'   rasters>, prior_probs=, lik=<value rasters>, lik_probs=)`.
#' @param mukey_raster,mukey_texture_draws Optional - opts into `prior_fusion_method =
#'   "raw_draws"` (see `fuse_texture_group_batch_core()`'s docs). `mukey_raster` must already be
#'   aligned to the same grid as `fetched`'s percentile rasters (nearest-neighbor resampled - it's
#'   categorical). `mukey_texture_draws` is a `lookup_mukey_texture_draws()` result. `NULL`
#'   (default) for either uses the percentile-reconstruction path.
#' @param posterior_probs Numeric probabilities (0-1) at which to report each fraction's posterior
#'   percentiles. `NULL` (default) resolves to [FUSE_POSTERIOR_DEFAULT_PROBS]. Computed via
#'   `texture_group_percentiles_raster()` - see that function's `@section Sum-to-100 caveat` for
#'   why, unlike `value`, these do NOT generally sum to 100 across the group's three members.
#' @param posterior_n_mc Monte Carlo draw count per cell for the percentile sampler (default 2000).
#'   See `texture_group_percentiles_raster()`'s docs for why this is a separate lever from the
#'   moment-estimation MC sample size used elsewhere in this function.
#' @return A named list keyed by each member's `id`:
#'   `list(posterior = list(value = <inverse-ILR point-estimate raster for
#'   that fraction>, ilr_mu = <2-layer raster, shared across the group>,
#'   ilr_Sigma = <3-layer raster (S11,S12,S22), shared>, percentiles = <named list of percentile
#'   SpatRasters for THIS member's own fraction - see the sum-to-100 caveat above>), dist =
#'   "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL,
#'   n_fallback_cells = 0)`. To draw posterior samples for a specific
#'   fraction/cell directly, extract that cell's `ilr_mu`/`ilr_Sigma` and pass to
#'   `sample_ilr_posterior()`.
#' @keywords internal
fuse_texture_group <- function(fetched, mukey_raster = NULL, mukey_texture_draws = NULL,
                                posterior_probs = NULL, posterior_n_mc = 2000) {
  stopifnot(length(fetched) == 3)
  ids <- vapply(fetched, function(f) f$id, character(1))

  prior_idx <- percentile_index(fetched[[1]]$prior_probs)
  lik_idx <- percentile_index(fetched[[1]]$lik_probs)
  prior_z <- stats::qnorm(prior_idx$p_hi)
  lik_z <- stats::qnorm(lik_idx$p_hi)

  stack_list <- list()
  for (f in fetched) {
    stack_list[[paste0(f$id, "_prior_lo")]] <- f$prior[[prior_idx$lo_idx]]
    stack_list[[paste0(f$id, "_prior_p50")]] <- f$prior[[prior_idx$p50_idx]]
    stack_list[[paste0(f$id, "_prior_hi")]] <- f$prior[[prior_idx$hi_idx]]
    stack_list[[paste0(f$id, "_lik_lo")]] <- f$lik[[lik_idx$lo_idx]]
    stack_list[[paste0(f$id, "_lik_p50")]] <- f$lik[[lik_idx$p50_idx]]
    stack_list[[paste0(f$id, "_lik_hi")]] <- f$lik[[lik_idx$hi_idx]]
  }
  combined <- terra::rast(stack_list)
  col_names <- names(stack_list)
  clay_id <- ids[1]; sand_id <- ids[2]; silt_id <- ids[3]

  # n_mc matches fuse_texture_group_batch()/fuse_texture_group_batch_core()'s own default - kept
  # as a local constant here (not a new exposed parameter) since fuse_texture_group() itself never
  # exposed n_mc before this either.
  texture_n_mc <- 2000

  use_raw_draws <- !is.null(mukey_raster) && !is.null(mukey_texture_draws)
  mukey_prior_texture_samples <- NULL
  if (use_raw_draws) {
    # Precompute, once per unique mukey (not once per cell), a fixed texture_n_mc-row resample of
    # that mukey's real JOINT (clay, sand, silt) draws - preserves real cross-fraction correlation
    # a per-fraction percentile reconstruction can't. Row-sampled (not per-column independently),
    # so each resampled row stays a genuine (clay, sand, silt) triple from one real simulated
    # cokey/replicate.
    mukey_prior_texture_samples <- lapply(mukey_texture_draws, function(mat) {
      if (is.null(mat) || nrow(mat) == 0) return(NULL)
      mat[sample(nrow(mat), texture_n_mc, replace = TRUE), , drop = FALSE]
    })
    mukey_raster_plain <- mukey_raster
    levels(mukey_raster_plain) <- NULL
    combined <- c(combined, mukey_raster_plain)
    col_names <- c(col_names, "mukey_code")
  }

  # estimate_ilr_moments_mc() (3x rnorm(n_mc=2000) + an ILR
  # transform + cov()) twice per raster cell via terra::app()+apply(), an R-level scalar
  # closure dispatched once per cell - the dominant cost for any non-trivial AOI (see
  #). fuse_texture_group_batch() below computes the
  # same thing for a whole chunk of cells at once via vectorized matrix arithmetic, with the
  # rnorm() draws ordered to consume the RNG stream in EXACTLY the same
  # cell-then-(prior/lik)-then-(clay/sand/silt)-then-mc-draw order the old per-cell loop did,
  # so results are bit-for-bit identical (see test-raster-fusion.R's regression test), not just
  # statistically similar. (This guarantee only applies to the default, non-raw_draws path.)
  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1, dimnames = list(NULL, col_names))
    fuse_texture_group_batch(row_mat, clay_id, sand_id, silt_id, prior_z, lik_z,
                              n_mc = texture_n_mc,
                              mukey_prior_texture_samples = mukey_prior_texture_samples)
  }
  result_r <- terra::app(combined, fun)
  names(result_r) <- c("mu1", "mu2", "S11", "S12", "S22")

  mu1_vals <- terra::values(result_r[["mu1"]])[, 1]
  mu2_vals <- terra::values(result_r[["mu2"]])[, 1]
  point_est <- ilr_inverse(mu1_vals, mu2_vals)

  effective_posterior_probs <- if (is.null(posterior_probs)) FUSE_POSTERIOR_DEFAULT_PROBS else posterior_probs
  fraction_percentiles <- texture_group_percentiles_raster(
    result_r[[c("mu1", "mu2")]], result_r[[c("S11", "S12", "S22")]],
    effective_posterior_probs, n_mc = posterior_n_mc
  )
  # Positional roles (fraction_percentiles$clay/sand/silt) match ilr_inverse()'s fixed column
  # order, which in turn matches how clay_id/sand_id/silt_id (= ids[1]/ids[2]/ids[3]) were assigned
  # above - so member i's own percentiles are fraction_percentiles at that same position.
  member_percentiles <- list(fraction_percentiles$clay, fraction_percentiles$sand, fraction_percentiles$silt)

  template <- fetched[[1]]$prior[[1]]
  make_layer <- function(v) { r <- template; terra::values(r) <- v; r }

  stats::setNames(lapply(seq_along(ids), function(i) {
    list(
      posterior = list(
        value = make_layer(point_est[, i]),
        ilr_mu = result_r[[c("mu1", "mu2")]],
        ilr_Sigma = result_r[[c("S11", "S12", "S22")]],
        percentiles = member_percentiles[[i]]
      ),
      dist = "texture_ilr",
      route = "closed_form_ilr_group",
      route_detail = NULL,
      n_fallback_cells = 0
    )
  }), ids)
}

# ---------------------------------------------------------------------------
# Top-level orchestrators: fetch (cached) SSURGO prior + SOLUS likelihood, align, fuse
# ---------------------------------------------------------------------------

#' Resolve whether a `run_fusion()` call should fuse against real per-mukey Monte Carlo
#' draws (`prior_fusion_method = "raw_draws"`) or the percentile-reconstructed approximation.
#'
#' An explicit `prior_fusion_method` value ("raw_draws" or "percentile") always wins; this
#' function only decides what happens when it is unset (`NULL`, the common case).
#'
#' Unset defaults to raw_draws for any `dist` whose fusion route has a raw-draws fit:
#' `"normal"`/`"beta"`/`"gamma"` (both AOI-size routes) and `"lognormal"` (both its small-AOI
#' general branch and large-AOI closed-form branch). `"auto"` is included, since
#' `resolve_property_dist()`'s auto-selection only ever resolves to `"beta"`/`"normal"`/
#' `"lognormal"`, never `"metalog"`. An explicit `dist = "metalog"` uses the percentile
#' reconstruction: `fuse_metalog_adapter()` has an empirical-percentile raw-draws variant, but
#' its cost/benefit is not benchmarked, so it is kept off by default.
#'
#' On the general-KDE route raw_draws costs about 16-25% more than the percentile path; on the
#' closed-form route it adds about 0.18s on top of the fit cost at 100,172 cells / 150 mukeys.
#'
#' @param prior_fusion_method `property_config$prior_fusion_method` - `NULL`/unset resolves per
#'   `dist` (see above); `"raw_draws"`/`"percentile"` are explicit overrides, always honored.
#' @param dist `property_config$dist` - the value as configured, BEFORE `resolve_property_dist()`
#'   resolves `"auto"` to a concrete family (this function runs earlier, before the mukey draws
#'   simulation that resolution would otherwise duplicate the cost of triggering).
#' @return `TRUE`/`FALSE`.
#' @keywords internal
resolve_want_raw_draws <- function(prior_fusion_method, dist) {
  if (identical(prior_fusion_method, "raw_draws")) return(TRUE)
  if (identical(prior_fusion_method, "percentile")) return(FALSE)
  isTRUE(dist %in% c("auto", "normal", "beta", "gamma", "lognormal"))
}

#' Align a Prior and Likelihood onto One Shared Grid
#'
#' Every fusion route requires cell-aligned prior/likelihood rasters, so `run_fusion()`,
#' `run_fusion_group()`, and `run_fusion_multiproperty()` each resample one side onto the
#' other's grid before fusing. Factored out here so the resampling *direction* is one
#' source-agnostic decision point (prior vs. likelihood, not "SSURGO" vs. "SOLUS") rather than
#' several copies of the same `terra::resample()` call.
#'
#' @param prior_values,lik_values Named lists of percentile `SpatRaster`s (a `prior$values` /
#'   `likelihood$values`-shaped list, as `run_fusion()` etc. hold them).
#' @param resampling `"down"` (default, matches all prior behavior exactly) resamples the prior
#'   DOWN onto the likelihood's grid via bilinear interpolation - the fused output sits at the
#'   likelihood's own resolution, which for the SSURGO/SOLUS100 pipeline means SSURGO's finer
#'   native grid (~30 m) is blurred/averaged into SOLUS100's coarser 100 m cells.
#'
#'   `"up"` instead resamples the likelihood UP onto the prior's grid via nearest-neighbor
#'   ("block replication" - every fine cell takes its enclosing coarse cell's exact value, with
#'   no invented intermediate values) - the fused output sits at the prior's own (finer)
#'   resolution instead, preserving the prior's real spatial detail (e.g. SSURGO map-unit
#'   boundaries) rather than discarding it. Nearest-neighbor, not bilinear, is deliberate here:
#'   bilinear would interpolate BETWEEN two independent likelihood-cell predictions, implying a
#'   smooth underlying gradient the likelihood source never actually asserted - it does not add
#'   real resolution to the likelihood, only fabricates it. Every fine cell's likelihood value
#'   under `"up"` is traceable back to one real published likelihood cell value.
#' @return `list(prior_values=, lik_values=, reference_grid=)` - `prior_values`/`lik_values` are
#'   both now on `reference_grid` (whichever input's own grid the other side was resampled onto),
#'   returned for reuse aligning e.g. a categorical mukey raster onto the same grid.
#' @keywords internal
align_prior_likelihood <- function(prior_values, lik_values, resampling = c("down", "up")) {
  resampling <- match.arg(resampling)
  if (resampling == "down") {
    reference_grid <- lik_values[[1]]
    list(prior_values = lapply(prior_values, terra::resample, y = reference_grid, method = "bilinear"),
         lik_values = lik_values,
         reference_grid = reference_grid)
  } else {
    reference_grid <- prior_values[[1]]
    list(prior_values = prior_values,
         lik_values = lapply(lik_values, terra::resample, y = reference_grid, method = "near"),
         reference_grid = reference_grid)
  }
}

#' Fuse an Already-Fetched Prior and Likelihood (Stage 1 fusion tail)
#'
#' The portion of [run_fusion()] that runs once its SSURGO prior and SOLUS likelihood
#' percentile lists are in hand: resample the prior onto the SOLUS grid, optionally attach the real
#' per-mukey Monte Carlo draws for raw-draws fusion, fuse via [fuse_property_adaptive()], and
#' assemble the standard `run_fusion()` return list. Factored out so `run_fusion()`
#' and `run_fusion_multiproperty()` share exactly one non-compositional fusion code path.
#'
#' @param property_config As in [run_fusion()].
#' @param prior,solus The pre-alignment `list(values=, probs=)` percentile lists for the SSURGO
#'   prior and SOLUS likelihood.
#' @param draws Optional in-memory draws data frame (already simulated by the caller) for the
#'   raw-draws fusion path; `NULL` skips it.
#' @param mukey_raster_native Optional native-grid mukey raster matching `draws` - required for the
#'   raw-draws path.
#' @param verbose Passed through to [fuse_property_adaptive()] (default `TRUE`, matching the
#'   original inline behavior).
#' @param resampling Passed through to [align_prior_likelihood()] - `"down"` (default) resamples
#'   the SSURGO prior onto SOLUS100's coarser grid (today's behavior); `"up"` instead resamples
#'   SOLUS100 onto SSURGO's finer native grid, preserving SSURGO's map-unit boundary resolution.
#'   See [align_prior_likelihood()]'s docs for the full discussion.
#' @return `run_fusion()`'s return list
#'   (`list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=,
#'   route_detail=, n_fallback_cells=)`).
#' @keywords internal
stage1_fuse_from_prior_solus <- function(property_config, prior, solus,
                                          draws = NULL, mukey_raster_native = NULL,
                                          verbose = TRUE, resampling = c("down", "up")) {
  resampling <- match.arg(resampling)
  ssurgo_property_id <- if (!is.null(property_config$solus_variable)) property_config$solus_variable else property_config$id
  want_raw_draws <- resolve_want_raw_draws(property_config$prior_fusion_method, property_config$dist)

  # SSURGO prior (~30m, from mukey.wcs()) and SOLUS likelihood (100m) are on different grids -
  # align them onto one shared grid before combining, since every fusion route requires
  # cell-aligned rasters. See align_prior_likelihood()'s docs for the "down"/"up" tradeoff.
  aligned <- align_prior_likelihood(prior$values, solus$values, resampling = resampling)
  prior_aligned <- aligned$prior_values
  solus_aligned_values <- aligned$lik_values
  reference_grid <- aligned$reference_grid

  # Opt-in raw-draws fusion (see run_fusion()'s @section): `draws` is already in memory from
  # the caller whenever want_raw_draws is TRUE, so this never triggers a further simulation. Align
  # the mukey raster onto the same reference grid via NEAREST-NEIGHBOR resampling (not bilinear,
  # which would fabricate nonsensical interpolated mukey codes between categories) - a no-op when
  # `resampling = "up"`, since the reference grid is already mukey_raster_native's own grid.
  extra_fusion_args <- list()
  if (want_raw_draws && !is.null(draws) && !is.null(mukey_raster_native)) {
    mukey_draws <- lookup_mukey_draws(draws, ssurgo_property_id)
    if (!is.null(mukey_draws)) {
      extra_fusion_args$mukey_raster <- terra::resample(mukey_raster_native, reference_grid, method = "near")
      extra_fusion_args$mukey_draws <- mukey_draws
    }
  }

  fused <- do.call(fuse_property_adaptive, c(
    list(
      prior_value_rasters = prior_aligned, prior_probs = prior$probs,
      lik_value_rasters = solus_aligned_values, lik_probs = solus$probs,
      property_config = property_config, verbose = verbose
    ),
    extra_fusion_args
  ))

  list(
    prior = list(values = prior_aligned, probs = prior$probs),
    likelihood = list(values = solus_aligned_values, probs = solus$probs),
    posterior = fused$posterior,
    dist = fused$dist, dist_source = fused$dist_source, skew_proxy = fused$skew_proxy,
    route = fused$route, route_detail = fused$route_detail, n_fallback_cells = fused$n_fallback_cells
  )
}

#' Run Stage 1 Fusion for One Property/Depth over an AOI
#'
#' Fetches the (cached) SSURGO prior (`R/adapter-ssurgo-simulate.R`'s `fetch_ssurgo_percentiles()`) and
#' SOLUS likelihood (`R/adapter-solus.R`'s `fetch_solus_percentiles()`) percentiles, aligns the
#' SSURGO grid onto the SOLUS grid via `terra::resample()`, then fuses via
#' `fuse_property_adaptive()`.
#'
#' Compositional properties (`property_config$composition_group` set) are detected up front and
#' delegated to `run_fusion_group()` instead, since they must be fetched/fit jointly.
#'
#' @param aoi_vect A `terra::SpatVector` AOI (projected, e.g. EPSG:5070).
#' @param property_config A list with `id` (used as `fetch_ssurgo_percentiles()`'s property id and
#'   as the cache key), `solus_variable` (a `soilDB::fetchSOLUS()`-recognized variable name), and
#'   `dist`/`bounds`/`boundedness`/`auto_skew_threshold`/`composition_group` as needed by
#'   `fuse_property_adaptive()`.
#' @param top_depth,bottom_depth Numeric depth bounds in cm.
#' @param composition_groups,property_configs Only needed when `property_config$composition_group`
#'   is set - passed through to `run_fusion_group()` (see its docs).
#' @param parallel,n_cores Passed through to `simulate_ssurgo_mapunit_draws()`'s
#'   `parallel`/`n_cores` (via `fetch_ssurgo_percentiles()` or `run_fusion_group()`,
#'   whichever path this call takes) - the SSURGO adapter's per-cokey depth-trend GP fitting is
#'   this pipeline's dominant cost for AOIs with many cokeys. Default `parallel = FALSE` matches
#'   prior behavior exactly.
#' @param seed Optional integer for opt-in determinism. `NULL` (default) keeps the current
#'   stochastic behavior. When set, `set.seed(seed)` runs once at the top of this call and `seed`
#'   is forwarded to `simulate_ssurgo_mapunit_draws()`, so the whole call (SSURGO simulation,
#'   depth-trend GP fit, and the posterior-sampling fusion step) is reproducible given identical
#'   upstream SSURGO/SOLUS data. See `simulate_ssurgo_mapunit_draws()`'s `seed` docs for the
#'   conditional-determinism caveats.
#' @param resampling `"down"` (default, matches all prior behavior exactly) resamples the SSURGO
#'   prior onto SOLUS100's coarser 100 m grid before fusing. `"up"` instead resamples SOLUS100
#'   onto SSURGO's finer native grid (nearest-neighbor block replication, not interpolation),
#'   so the fused output preserves SSURGO's real map-unit boundary resolution instead of blurring
#'   it away - at the cost of ~(100/native_res)^2 more output cells and no gain in SOLUS100's own
#'   information content (its real resolving power stays 100 m either way). See
#'   [align_prior_likelihood()]'s docs for the full discussion.
#' @return `list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=,
#'   route_detail=, n_fallback_cells=)`, or `NULL` if the SSURGO or SOLUS side failed.
#'   `posterior`'s shape depends on `dist` - see `fuse_property_adaptive()`'s docs. As of
#', `posterior` also carries a
#'   `percentiles` field (a named list of `SpatRaster`s, `"P01"`..`"P99"` by default - see
#'   [FUSE_POSTERIOR_DEFAULT_PROBS]) regardless of `dist`, for every non-compositional property -
#'   this function doesn't expose its own `posterior_probs` parameter (matching the existing
#'   convention that `n_samples`/`grid_resolution` also aren't threaded this far down; the
#'   underlying `fuse_*()` route's own default applies).
#' @section SSURGO simulation scope:
#' The SSURGO prior simulation is restricted to just this call's own property (plus the full
#' sand/silt/clay draw for any texture member) via `simulate_ssurgo_mapunit_draws()`'s
#' `requested_properties` - the per-cokey depth-trend GP fit is the pipeline's dominant cost and
#' scales with the property count. Consequences: (1) the prior is *statistically equivalent* to,
#' not bit-identical to, a full-property simulation's matching column - the pipeline is unseeded by
#' default (a fresh draw either way), and even with `seed` set, a subset consumes the RNG stream
#' differently from a full run; (2) a map unit whose components carry no data for *this* property
#' does not gain from those components' other properties being present, so a single-property
#' prior - and hence the fused posterior - can be `NA` for a few cells a full-property run would
#' have covered. Those cells genuinely lack data for the property in question.
#' @section Fusion fidelity - `prior_fusion_method`:
#' Controls whether fusion runs against each mukey's real per-cell Monte Carlo draws
#' (`"raw_draws"`) or a percentile-reconstructed approximation (`"percentile"`) - see
#' `fuse_general_kde()`'s docs. Leaving `property_config$prior_fusion_method` unset defaults to
#' `"raw_draws"` for `dist %in% c("auto","normal","beta","gamma","lognormal")` - see
#' `resolve_want_raw_draws()` for the exact rule. An explicit `dist = "metalog"` uses
#' `"percentile"`. Set `property_config$prior_fusion_method = "percentile"` explicitly to force
#' the percentile path for any `dist` (e.g. to avoid the live simulation dependency);
#' `"raw_draws"` forces it on even for `dist = "metalog"`.
#'
#' The real draws are held in memory only for the duration of this one call, never disk-cached
#' (see `simulate_ssurgo_mapunit_draws()`'s docs for why) - so raw_draws fusion (whether defaulted
#' or explicit) forces a fresh simulation even when this AOI/property/depth-window's `"ssurgo"`
#' percentile cache is already warm (the draws that produced that cached percentile summary weren't
#' kept). That simulation still runs at most once per call regardless of cache state - it is shared
#' with the percentile-cache-population step when both are needed in the same call.
#' @keywords internal
run_fusion <- function(aoi_vect, property_config, top_depth, bottom_depth,
                               composition_groups = NULL, property_configs = NULL,
                               parallel = FALSE, n_cores = NULL, seed = NULL,
                               resampling = c("down", "up")) {
  resampling <- match.arg(resampling)
  # Opt-in determinism: seed once up front so the whole call - including a fusion-sampling step
  # that runs even on a warm percentile cache - is reproducible. simulate_ssurgo_mapunit_draws()
  # re-seeds to the same value for its own scope (and the parallel depth-trend path); nothing
  # between here and there consumes the RNG.
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(property_config$composition_group)) {
    group_result <- run_fusion_group(
      aoi_vect, property_config$composition_group, composition_groups, property_configs,
      top_depth, bottom_depth, parallel = parallel, n_cores = n_cores, seed = seed,
      resampling = resampling
    )
    return(if (is.null(group_result)) NULL else new_soilSIM_fusion(group_result[[property_config$id]]))
  }

  ssurgo_property_id <- if (!is.null(property_config$solus_variable)) property_config$solus_variable else property_config$id
  want_raw_draws <- resolve_want_raw_draws(property_config$prior_fusion_method, property_config$dist)

  ssurgo_key <- build_cache_key(aoi_vect, property_config$id, top_depth, bottom_depth, "ssurgo")
  prior <- cache_get_valid_percentiles(ssurgo_key)

  # The mukey grid/draws are computed AT MOST ONCE and reused in-memory for both the percentile
  # cache (below) and the raw-draws lookup (further down) - needed whenever the percentile cache
  # is cold (have to simulate anyway) OR raw-draws fusion was explicitly requested (needs the real
  # draws regardless of whether percentiles are already cached). This avoids literally re-running
  # the expensive simulation twice in one call, which an earlier version of this function did.
  # simulate_ssurgo_mapunit_draws() itself is deliberately NOT disk-cached - see its own docs.
  mukey_raster_native <- NULL
  draws <- NULL
  if (is.null(prior) || want_raw_draws) {
    mukey_raster_native <- fetch_ssurgo_mukey_raster(aoi_vect)
    if (is.null(mukey_raster_native)) return(NULL)
    # Simulate only the one property this call fuses (plus texture coupling) - the per-cokey GP
    # depth-trend fit is the pipeline's dominant cost and scales with the property count. Only
    # ever law-equivalent (not bit-identical) to a full simulation's matching column, and only
    # under the default joint_copula vertical method - see simulate_ssurgo_mapunit_draws()'s
    # `requested_properties` docs. run_fusion() never forwards `config`, so it always gets
    # joint_copula. A cokey with no data for this one property no longer rides along on the
    # others, so a single-property prior can be NA for a mukey a full-property run would have
    # covered - documented in this function's @return.
    draws <- simulate_ssurgo_mapunit_draws(aoi_vect, top_depth, bottom_depth,
                                            parallel = parallel, n_cores = n_cores,
                                            mukey_raster = mukey_raster_native,
                                            requested_properties = ssurgo_property_id, seed = seed)
    if (is.null(draws)) return(NULL)
  }

  if (is.null(prior)) {
    prior <- percentiles_from_draws(mukey_raster_native, draws, ssurgo_property_id)
    if (is.null(prior)) return(NULL)
    cache_set(ssurgo_key, "ssurgo", wrap_percentile_list(prior))
  }

  solus_key <- build_cache_key(aoi_vect, property_config$id, top_depth, bottom_depth, "solus")
  solus <- cache_get_valid_percentiles(solus_key)
  if (is.null(solus)) {
    solus <- fetch_solus_percentiles(aoi_vect, property_config$solus_variable, top_depth, bottom_depth)
    if (is.null(solus)) return(NULL)
    cache_set(solus_key, "solus", wrap_percentile_list(solus))
  }

  # `draws`/`mukey_raster_native` are already in memory from above whenever want_raw_draws is TRUE
  # (part of the is.null(prior) || want_raw_draws condition), so the shared fusion tail never
  # triggers a further simulation.
  new_soilSIM_fusion(stage1_fuse_from_prior_solus(
    property_config, prior, solus,
    draws = draws, mukey_raster_native = mukey_raster_native, resampling = resampling
  ))
}

#' Fuse an Already-Fetched Texture Group (Stage 1 group-fusion tail)
#'
#' The portion of [run_fusion_group()] that runs once every member's aligned SSURGO prior
#' and SOLUS likelihood are assembled into `fetched`: optionally prepare the real joint per-mukey
#' texture draws for raw-draws fusion, then fuse the group via [fuse_texture_group()]. Factored out
#' so `run_fusion_group()` and `run_fusion_multiproperty()` share one group-fusion code path.
#'
#' @param fetched A list (one entry per group member, in `group_members()` order) of
#'   `list(id=, prior=<aligned percentile-value rasters>, prior_probs=, lik=<value rasters>,
#'   lik_probs=)`.
#' @param want_raw_draws Whether any member requested `prior_fusion_method = "raw_draws"` - the
#'   raw-draws texture path only runs when this is `TRUE` (a cold percentile cache alone can leave
#'   `shared_draws` non-`NULL` without raw-draws fusion being requested).
#' @param shared_draws,shared_mukey_raster The one-per-group in-memory draws data frame and its
#'   native-grid mukey raster (already computed by the caller); both required for the raw-draws
#'   path.
#' @param reference_grid The grid every `fetched` member's prior/likelihood was already aligned
#'   onto (i.e. `align_prior_likelihood()`'s `reference_grid` from aligning the first member) -
#'   the mukey raster is aligned onto this same grid, regardless of which side (prior or
#'   likelihood) it came from. Defaults to `fetched[[1]]$lik[[1]]` for backward compatibility
#'   with callers that haven't been updated to pass it explicitly (matches the original
#'   `resampling = "down"`-only behavior).
#' @return `fuse_texture_group()`'s return value - a named list keyed by member id.
#' @keywords internal
stage1_fuse_texture_group_from_fetched <- function(fetched, want_raw_draws = FALSE,
                                                    shared_draws = NULL, shared_mukey_raster = NULL,
                                                    reference_grid = fetched[[1]]$lik[[1]]) {
  # Align the mukey raster onto the same reference grid every member's prior/likelihood was
  # already aligned onto, via NEAREST-NEIGHBOR resampling (not bilinear, which would fabricate
  # nonsensical interpolated mukey codes between categories) - a no-op when that grid is already
  # the mukey raster's own native grid (the "up" resampling direction).
  mukey_texture_draws <- NULL
  mukey_raster_aligned <- NULL
  if (isTRUE(want_raw_draws) && !is.null(shared_draws) && !is.null(shared_mukey_raster)) {
    mukey_texture_draws <- lookup_mukey_texture_draws(shared_draws)
    if (!is.null(mukey_texture_draws)) {
      mukey_raster_aligned <- terra::resample(shared_mukey_raster, reference_grid, method = "near")
    }
  }

  fuse_texture_group(fetched, mukey_raster = mukey_raster_aligned,
                     mukey_texture_draws = mukey_texture_draws)
}

#' Run Stage 1 Fusion for a Whole Compositional Group Jointly
#'
#' Fuses a compositional group's members (currently only `"texture"`: clay/sand/silt) jointly via
#' `fuse_texture_group()`, instead of one independent `run_fusion()` call per member -
#' independent fusion measurably breaks sum-to-100. Each member's SSURGO/SOLUS percentiles are
#' still fetched/cached per-property (so a later request for just one member's raw percentiles
#' still hits cache), but the fusion itself runs once for the whole group, under a `"texture_group"`
#' cache kind (`"texture_group_raw_draws"` when `prior_fusion_method = "raw_draws"` was requested -
#' see the section below for why these must be distinct kinds), with each member's resulting
#' posterior also seeded into a per-property `"posterior"`/`"posterior_raw_draws"` cache kind - so
#' three sequential `run_fusion()` calls for clay, then sand, then silt trigger the joint
#' fetch+fusion exactly once, not three times.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param group A composition group name (e.g. `"texture"`).
#' @param composition_groups A `config$monte_carlo$composition_groups`-shaped list (see
#'   `group_members()`) giving this group's member ids in order.
#' @param property_configs A named list, keyed by each member id `composition_groups` lists,
#'   of per-member config lists (each needs at least `id`/`solus_variable`) - the per-call
#'   replacement for the original bundle's global `PROPERTIES` registry.
#' @param top_depth,bottom_depth Numeric depth bounds in cm.
#' @param parallel,n_cores Passed through to `simulate_ssurgo_mapunit_draws()`'s
#'   `parallel`/`n_cores` for the shared draws computation below (only relevant when at least one
#'   member isn't already disk-cached). Default `parallel = FALSE` matches prior behavior exactly.
#' @param seed Optional integer for opt-in determinism, forwarded to
#'   `simulate_ssurgo_mapunit_draws()` (and `set.seed()` once up front). `NULL` (default) =
#'   current stochastic behavior.
#' @param resampling Passed through to [align_prior_likelihood()] for every member's
#'   prior/likelihood alignment - `"down"` (default) matches all prior behavior; `"up"`
#'   preserves SSURGO's native grid instead. See [align_prior_likelihood()]'s docs.
#' @return The full group result: a named list keyed by member id, each element
#'   `list(posterior = list(value=, ilr_mu=, ilr_Sigma=, percentiles=), dist = "texture_ilr", route =
#'   "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)` (see
#'   `fuse_texture_group()`'s docs - NOT the uniform `(mu,sigma)`/`(alpha,beta)` contract
#'   non-compositional properties get) - or `NULL` if any member's SSURGO/SOLUS fetch failed.
#'   `run_fusion()`'s own dispatch slices this down to the single requested member
#'   (`group_result[[property_config$id]]`) to keep its own per-property return contract
#'   consistent regardless of `dist`. `percentiles` is THIS member's own fraction's
#'   percentiles (a named list of `SpatRaster`s) - see `fuse_texture_group()`'s `@section
#'   Sum-to-100 caveat`: unlike `value`, these do NOT generally sum to 100 across the group's three
#'   members.
#' @section Higher-fidelity fusion (opt-in):
#' Set `prior_fusion_method = "raw_draws"` on ANY member's `property_configs` entry to fuse the
#' whole group against the real joint (clay, sand, silt) Monte Carlo draws instead of independent
#' percentile-reconstructed marginals - see `fuse_texture_group()`'s docs. It is checked across
#' all members (not just one), since the group is always fused jointly regardless of which
#' member's config carries the setting. This forces the shared simulation to run even when every
#' member's own `"ssurgo"` percentile cache is already warm (the draws that produced those cached
#' percentiles are not kept - see `simulate_ssurgo_mapunit_draws()`'s docs). Unset on every
#' member (default) uses the percentile-reconstruction path; unlike `run_fusion()`, the
#' joint texture-group route does not default to raw draws.
#' @keywords internal
run_fusion_group <- function(aoi_vect, group, composition_groups, property_configs,
                                     top_depth, bottom_depth, parallel = FALSE, n_cores = NULL,
                                     seed = NULL, resampling = c("down", "up")) {
  resampling <- match.arg(resampling)
  if (!is.null(seed)) set.seed(seed)
  member_ids <- group_members(group, composition_groups)
  members <- property_configs[member_ids]
  want_raw_draws <- any(vapply(members, function(m) identical(m$prior_fusion_method, "raw_draws"), logical(1)))

  # The "kind" strings below encode want_raw_draws (e.g. "texture_group" vs
  # "texture_group_raw_draws") so a cached FUSED result from one prior_fusion_method is never
  # silently served back to a caller requesting the other: without
  # this, switching prior_fusion_method for the same AOI/depth-window returned a stale,
  # method-mismatched cached group result (the group-level cache stores the fused POSTERIOR, unlike
  # run_fusion()'s own "ssurgo" cache, which only stores the pre-fusion PRIOR - re-fused
  # fresh every call - so this staleness risk is specific to the group path).
  group_kind <- if (want_raw_draws) "texture_group_raw_draws" else "texture_group"
  member_kind <- if (want_raw_draws) "posterior_raw_draws" else "posterior"

  group_key <- build_cache_key(aoi_vect, group, top_depth, bottom_depth, group_kind)
  # cache_set()'s own docs flag that SpatRasters must be wrap()/unwrap()ed around a cache
  # round-trip - group_result nests them inside each member's $posterior$value/ilr_mu/ilr_Sigma,
  # not the simple list(values=, probs=) shape wrap_percentile_list() handles, so the generic
  # wrap_nested_rasters()/unwrap_nested_rasters() pair is needed here instead. A malformed hit
  # (see cache_get_valid_percentiles()'s docs for why this matters) is discarded rather than
  # trusted, via cache_get_valid_group_result()'s shape check.
  group_result <- cache_get_valid_group_result(group_key, member_ids)

  if (is.null(group_result)) {
    # Each member's own SSURGO percentile cache is still checked first (so a later single-
    # property request, or a partially-cached group, still hits cache per member) - but if any
    # member is missing, run the shared mukey-raster fetch + Monte Carlo simulation ONCE here
    # rather than once per member. simulate_cokey_generalized() already simulates every
    # recognized property (clay/sand/silt/db/ph/...) jointly per cokey in one pass, so calling
    # fetch_ssurgo_percentiles() independently per texture member would re-run that expensive
    # simulation (the dominant cost of the whole fusion pipeline) 3 times for a 3-member group,
    # to extract 3 columns a single simulation pass already produces together.
    # percentiles_from_draws() (R/adapter-ssurgo-simulate.R) is the shared quantile/
    # rasterize step, factored out of fetch_ssurgo_percentiles() for exactly this reuse.
    ssurgo_keys <- lapply(members, function(m) build_cache_key(aoi_vect, m$id, top_depth, bottom_depth, "ssurgo"))
    ssurgo_cached <- lapply(ssurgo_keys, cache_get_valid_percentiles)

    # Restrict the shared simulation to this group's members - for texture that normalizes to
    # "just the sand/silt/clay draw" (both ILR axes), skipping the other ~6 properties' GP fits.
    member_props <- unique(vapply(members, function(m) {
      if (!is.null(m$solus_variable)) m$solus_variable else m$id
    }, character(1)))

    shared_mukey_raster <- NULL
    shared_draws <- NULL
    if (want_raw_draws || any(vapply(ssurgo_cached, is.null, logical(1)))) {
      shared_mukey_raster <- fetch_ssurgo_mukey_raster(aoi_vect)
      if (!is.null(shared_mukey_raster)) {
        shared_draws <- simulate_ssurgo_mapunit_draws(aoi_vect, top_depth, bottom_depth,
                                                        parallel = parallel, n_cores = n_cores,
                                                        mukey_raster = shared_mukey_raster,
                                                        requested_properties = member_props,
                                                        seed = seed)
      }
    }

    fetched <- lapply(seq_along(members), function(idx) {
      m <- members[[idx]]

      prior <- ssurgo_cached[[idx]]
      if (is.null(prior)) {
        prior <- if (is.null(shared_mukey_raster) || is.null(shared_draws)) {
          NULL
        } else {
          ssurgo_property_id <- if (!is.null(m$solus_variable)) m$solus_variable else m$id
          percentiles_from_draws(shared_mukey_raster, shared_draws, ssurgo_property_id)
        }
        if (!is.null(prior)) cache_set(ssurgo_keys[[idx]], "ssurgo", wrap_percentile_list(prior))
      }

      solus_key <- build_cache_key(aoi_vect, m$id, top_depth, bottom_depth, "solus")
      solus <- cache_get_valid_percentiles(solus_key)
      if (is.null(solus)) {
        solus <- fetch_solus_percentiles(aoi_vect, m$solus_variable, top_depth, bottom_depth)
        if (!is.null(solus)) cache_set(solus_key, "solus", wrap_percentile_list(solus))
      }

      if (is.null(prior) || is.null(solus)) return(NULL)
      aligned <- align_prior_likelihood(prior$values, solus$values, resampling = resampling)
      list(id = m$id, prior = aligned$prior_values, prior_probs = prior$probs,
           lik = aligned$lik_values, lik_probs = solus$probs, reference_grid = aligned$reference_grid)
    })
    if (any(vapply(fetched, is.null, logical(1)))) return(NULL)

    # `shared_draws`/`shared_mukey_raster` are already in memory from above whenever want_raw_draws
    # is TRUE, so the shared texture-group fusion tail never triggers a further simulation.
    group_result <- stage1_fuse_texture_group_from_fetched(
      fetched, want_raw_draws = want_raw_draws,
      shared_draws = shared_draws, shared_mukey_raster = shared_mukey_raster,
      reference_grid = fetched[[1]]$reference_grid
    )
    cache_set(group_key, group_kind, wrap_nested_rasters(group_result))

    for (id in member_ids) {
      member_key <- build_cache_key(aoi_vect, id, top_depth, bottom_depth, member_kind)
      cache_set(member_key, member_kind, wrap_nested_rasters(group_result[[id]]))
    }
  }

  group_result
}

#' Run Stage 1 Fusion for Many Properties over an AOI in One Simulation Pass
#'
#' The multi-property / multi-depth-window counterpart of [run_fusion()]: runs the
#' expensive SSURGO Monte Carlo simulation **once** - covering every requested property and every
#' depth window via [simulate_ssurgo_mapunit_draws()]'s `depth_windows` argument - then fuses each
#' property/window against SOLUS, seeding the same per-property `"ssurgo"`/`"solus"`/`"posterior"`
#' disk caches [run_fusion()] reads. A later single-property `run_fusion()` call for
#' the same AOI/property/window therefore hits cache instead of re-simulating.
#'
#' Calling `run_fusion()` in a loop over N properties x M windows runs the (dominant-cost)
#' per-cokey simulation N*M times, each discarding all but one of the ~9 jointly-simulated
#' properties; this runs it once. Every leaf shares that one draw set, so - unlike independent
#' `run_fusion()` calls - their `NA` masks are mutually consistent.
#'
#' @param aoi_vect A `terra::SpatVector` AOI (projected, e.g. EPSG:5070).
#' @param property_configs A **named** list, keyed by each config's own `id`, of `property_config`
#'   lists exactly as [run_fusion()] takes them (each needs `id`/`solus_variable`, optionally
#'   `dist`/`bounds`/`composition_group`/`prior_fusion_method`/...). Configs carrying a
#'   `composition_group` are fused jointly per [run_fusion_group()]; **every** member of any
#'   referenced group must have its own entry here.
#' @param depth_windows A non-empty list of `c(top, bottom)` numeric pairs (cm, `bottom > top`).
#' @param composition_groups A `config$monte_carlo$composition_groups`-shaped list (see
#'   [group_members()]); required if any config sets `composition_group`.
#' @param n_mc,parallel,n_cores Passed to [simulate_ssurgo_mapunit_draws()].
#' @param seed Optional integer for opt-in determinism. `NULL` (default) = current stochastic
#'   behavior. When set, `set.seed(seed)` runs once up front and `seed` is forwarded to the one
#'   shared [simulate_ssurgo_mapunit_draws()] call, so the whole batch is reproducible given
#'   identical upstream SSURGO/SOLUS data.
#' @param verbose Passed to [fuse_property_adaptive()] - default `FALSE` (a multi-leaf batch is
#'   otherwise noisy with per-route messages).
#' @param simplify If `TRUE`, each leaf is reduced to `list(percentiles = <named SpatRasters>)`
#'   (what per-pixel re-marginalization consumes); default `FALSE` returns the full
#'   [run_fusion()]-shaped list per leaf.
#' @param resampling Passed through to [align_prior_likelihood()] for every leaf's
#'   prior/likelihood alignment - `"down"` (default) matches all prior behavior; `"up"`
#'   preserves SSURGO's native grid instead. See [align_prior_likelihood()]'s docs.
#' @return A nested named list `result[[config_id]][["<top>-<bottom>"]]`. Each leaf is a
#'   [run_fusion()] return list (or, for group members, the `texture_ilr`-shaped list from
#'   [run_fusion_group()]), or `NULL` on that leaf's own failure; the whole call returns
#'   `NULL` only if the shared simulation itself fails.
#' @seealso [run_fusion()], [run_fusion_group()], [extract_mukey_joint_ensemble()]
#' @keywords internal
run_fusion_multiproperty <- function(aoi_vect, property_configs, depth_windows,
                                     composition_groups = NULL,
                                     n_mc = 1000, parallel = FALSE, n_cores = NULL,
                                     seed = NULL, verbose = FALSE, simplify = FALSE,
                                     resampling = c("down", "up")) {
  resampling <- match.arg(resampling)
  if (!is.null(seed)) set.seed(seed)

  ## --- validate ---------------------------------------------------------------
  if (!is.list(property_configs) || length(property_configs) == 0 ||
      is.null(names(property_configs)) || any(!nzchar(names(property_configs)))) {
    stop("run_fusion_multiproperty(): `property_configs` must be a non-empty named list keyed by config id.")
  }
  if (!all(vapply(seq_along(property_configs),
                  function(i) identical(property_configs[[i]]$id, names(property_configs)[i]),
                  logical(1)))) {
    stop("run_fusion_multiproperty(): each `property_configs` entry's $id must match its list name.")
  }
  ok_pair <- function(w) length(w) == 2 && is.numeric(w) && is.finite(w[[1]]) &&
    is.finite(w[[2]]) && w[[2]] > w[[1]]
  if (!is.list(depth_windows) || length(depth_windows) == 0 ||
      !all(vapply(depth_windows, ok_pair, logical(1)))) {
    stop("run_fusion_multiproperty(): `depth_windows` must be a non-empty list of c(top, bottom) pairs with bottom > top.")
  }
  window_names <- vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), character(1))

  ## --- partition standalone vs compositional-group configs -------------------
  is_group <- vapply(property_configs, function(cfg) !is.null(cfg$composition_group), logical(1))
  group_names <- unique(vapply(property_configs[is_group],
                               function(cfg) cfg$composition_group, character(1)))
  if (length(group_names) && is.null(composition_groups)) {
    stop("run_fusion_multiproperty(): a config sets `composition_group` but `composition_groups` is NULL.")
  }
  for (g in group_names) {
    miss <- setdiff(group_members(g, composition_groups), names(property_configs))
    if (length(miss)) {
      stop(sprintf("run_fusion_multiproperty(): composition group '%s' member(s) missing from property_configs: %s",
                   g, paste(miss, collapse = ", ")))
    }
  }
  standalone_ids <- names(property_configs)[!is_group]

  ## --- decide whether the shared simulation is needed at all -----------------
  cfg_id_prop <- function(cfg) if (!is.null(cfg$solus_variable)) cfg$solus_variable else cfg$id
  need_sim <- any(vapply(property_configs,
                         function(cfg) resolve_want_raw_draws(cfg$prior_fusion_method, cfg$dist),
                         logical(1)))
  if (!need_sim) {
    for (w in depth_windows) {
      for (cfg in property_configs) {
        k <- build_cache_key(aoi_vect, cfg$id, w[[1]], w[[2]], "ssurgo")
        if (is.null(cache_get_valid_percentiles(k))) { need_sim <- TRUE; break }
      }
      if (need_sim) break
    }
  }

  ## --- one simulation, every window -----------------------------------------
  mukey_raster <- fetch_ssurgo_mukey_raster(aoi_vect)
  if (is.null(mukey_raster)) return(NULL)

  draws_by_window <- NULL
  if (need_sim) {
    req_props <- unique(unlist(lapply(property_configs, cfg_id_prop), use.names = FALSE))
    span <- range(unlist(depth_windows))
    draws_by_window <- simulate_ssurgo_mapunit_draws(
      aoi_vect, top_depth = span[1], bottom_depth = span[2],
      n_mc = n_mc, parallel = parallel, n_cores = n_cores,
      mukey_raster = mukey_raster, depth_windows = depth_windows,
      requested_properties = req_props, seed = seed
    )
    if (is.null(draws_by_window)) return(NULL)
  }

  ## --- SOLUS fetch memo (per solus_variable x window) -----------------------
  # S1: one batched fetch_solus_percentiles_multiproperty() call per window covers every variable needed
  # by any standalone config or group member, instead of one fetch_solus_percentiles() call per
  # (variable, window). Falls back to the scalar per-variable path only if the batched request
  # itself errors outright (get_solus()'s memo-miss path below already does that transparently).
  vars_needed <- unique(c(
    vapply(property_configs[standalone_ids], cfg_id_prop, character(1)),
    unlist(lapply(group_names, function(g) {
      vapply(property_configs[group_members(g, composition_groups)], cfg_id_prop, character(1))
    }), use.names = FALSE)
  ))

  solus_memo <- list()
  get_solus <- function(solus_variable, w) {
    key <- paste0(solus_variable, "@", w[[1]], "-", w[[2]])
    if (!is.null(solus_memo[[key]])) {
      return(if (identical(solus_memo[[key]], "MISS")) NULL else solus_memo[[key]])
    }
    val <- fetch_solus_percentiles(aoi_vect, solus_variable, w[[1]], w[[2]])
    solus_memo[[key]] <<- if (is.null(val)) "MISS" else val
    val
  }
  prefetch_solus_batch <- function(w) {
    batch <- tryCatch(
      fetch_solus_percentiles_multiproperty(aoi_vect, vars_needed, w[[1]], w[[2]]),
      error = function(e) NULL
    )
    if (is.null(batch)) return(invisible(NULL))  # defensive: get_solus() falls back per-variable
    for (v in vars_needed) {
      key <- paste0(v, "@", w[[1]], "-", w[[2]])
      solus_memo[[key]] <<- if (is.null(batch[[v]])) "MISS" else batch[[v]]
    }
    invisible(NULL)
  }
  seed_solus_cache <- function(cfg_id, w, solus) {
    k <- build_cache_key(aoi_vect, cfg_id, w[[1]], w[[2]], "solus")
    if (is.null(cache_get_valid_percentiles(k))) cache_set(k, "solus", wrap_percentile_list(solus))
  }

  result <- stats::setNames(lapply(names(property_configs), function(id) {
    stats::setNames(vector("list", length(window_names)), window_names)
  }), names(property_configs))

  ## --- per window: standalone leaves, then compositional groups -------------
  for (wi in seq_along(depth_windows)) {
    w <- depth_windows[[wi]]
    wname <- window_names[wi]
    prefetch_solus_batch(w)
    dw <- if (is.null(draws_by_window)) NULL else draws_by_window[[wname]]

    for (id in standalone_ids) {
      cfg <- property_configs[[id]]
      ssid <- cfg_id_prop(cfg)
      cfg_raw <- resolve_want_raw_draws(cfg$prior_fusion_method, cfg$dist)
      skey <- build_cache_key(aoi_vect, cfg$id, w[[1]], w[[2]], "ssurgo")

      # raw-draws leaf: skip the percentile-cache READ so prior + draws come from one sim.
      prior <- if (cfg_raw) NULL else cache_get_valid_percentiles(skey)
      if (is.null(prior)) {
        if (is.null(dw)) next
        prior <- percentiles_from_draws(mukey_raster, dw, ssid)
        if (is.null(prior)) next
        cache_set(skey, "ssurgo", wrap_percentile_list(prior))
      }
      solus <- get_solus(cfg$solus_variable, w)
      if (is.null(solus)) next
      seed_solus_cache(cfg$id, w, solus)

      leaf <- stage1_fuse_from_prior_solus(
        cfg, prior, solus,
        draws = if (cfg_raw) dw else NULL,
        mukey_raster_native = mukey_raster, verbose = verbose, resampling = resampling
      )
      result[[id]][[wname]] <- if (isTRUE(simplify)) list(percentiles = leaf$posterior$percentiles) else leaf
    }

    for (g in group_names) {
      member_ids <- group_members(g, composition_groups)
      members <- property_configs[member_ids]
      want_raw <- any(vapply(members, function(m) identical(m$prior_fusion_method, "raw_draws"), logical(1)))
      group_kind  <- if (want_raw) "texture_group_raw_draws" else "texture_group"
      member_kind <- if (want_raw) "posterior_raw_draws" else "posterior"

      fetched <- lapply(members, function(m) {
        ssid <- cfg_id_prop(m)
        skey <- build_cache_key(aoi_vect, m$id, w[[1]], w[[2]], "ssurgo")
        prior <- if (want_raw) NULL else cache_get_valid_percentiles(skey)
        if (is.null(prior)) {
          if (is.null(dw)) return(NULL)
          prior <- percentiles_from_draws(mukey_raster, dw, ssid)
          if (!is.null(prior)) cache_set(skey, "ssurgo", wrap_percentile_list(prior))
        }
        solus <- get_solus(m$solus_variable, w)
        if (is.null(prior) || is.null(solus)) return(NULL)
        seed_solus_cache(m$id, w, solus)
        aligned <- align_prior_likelihood(prior$values, solus$values, resampling = resampling)
        list(id = m$id, prior = aligned$prior_values, prior_probs = prior$probs,
             lik = aligned$lik_values, lik_probs = solus$probs, reference_grid = aligned$reference_grid)
      })
      if (any(vapply(fetched, is.null, logical(1)))) next

      group_result <- stage1_fuse_texture_group_from_fetched(
        fetched, want_raw_draws = want_raw,
        shared_draws = if (want_raw) dw else NULL,
        shared_mukey_raster = if (want_raw) mukey_raster else NULL,
        reference_grid = fetched[[1]]$reference_grid
      )
      cache_set(build_cache_key(aoi_vect, g, w[[1]], w[[2]], group_kind), group_kind,
                wrap_nested_rasters(group_result))
      for (mid in member_ids) {
        cache_set(build_cache_key(aoi_vect, mid, w[[1]], w[[2]], member_kind), member_kind,
                  wrap_nested_rasters(group_result[[mid]]))
        if (mid %in% names(result)) {
          result[[mid]][[wname]] <- if (isTRUE(simplify)) {
            list(percentiles = group_result[[mid]]$posterior$percentiles)
          } else {
            group_result[[mid]]
          }
        }
      }
    }

    if (!is.null(draws_by_window)) draws_by_window[[wname]] <- NULL  # release as we go
  }

  rm(draws_by_window)
  gc(verbose = FALSE)
  result
}


# ============================================================================
# merged from core-fusion.R (P1 file reorg)
# ============================================================================

#' @title Per-Pixel Ensemble Re-Marginalization (raster-fusion <-> tabular MC bridge)
#'
#' @description Transforms the per-mukey joint Monte Carlo ensemble from
#'   `extract_mukey_joint_ensemble()` (`R/adapter-ssurgo-simulate.R`) onto the SOLUS grid so that each
#'   pixel's marginal distributions match `run_fusion()`'s SOLUS-fused posterior, while
#'   preserving the ensemble's empirical copula (cross-property and cross-depth-window rank
#'   structure). This is the "per-pixel" answer to wiring raster fusion into the modelling
#'   framework.
#'
#'   `zonal_distribution_from_posterior()` is the cheap mukey-collapsed comparison arm (feeds
#'   `observed_data_by_mukey` in `simulate_monte_carlo()`); `remarginalize_ensemble_to_posterior()`
#'   is the per-pixel transform; `remarginalize_awc()` is the first derived-quantity consumer
#'   (available water capacity, via `saxton_rawls_raster()` since SOLUS100 has no water-retention
#'   variable). Nothing here touches either existing pipeline - purely additive.
#' @name raster_fusion_bridge
NULL

#' Probabilities behind a `run_fusion()` posterior's percentile layers
#' @param posterior A `run_fusion()`-style `list(percentiles = <named P.. rasters>, ...)`.
#' @return Numeric probabilities in (0, 1), ascending, parsed from the `"P<pct>"` layer names.
#' @keywords internal
posterior_percentile_probs <- function(posterior) {
  nm <- names(posterior$percentiles)
  if (is.null(nm) || !all(grepl("^P[0-9]+$", nm))) {
    stop("posterior_percentile_probs(): posterior$percentiles must be a list named 'P1'/'P50'/... .")
  }
  p <- as.numeric(sub("^P", "", nm)) / 100
  if (is.unsorted(p)) stop("posterior_percentile_probs(): percentile layers are not in ascending order.")
  p
}

#' Zonal (per-mukey) Reduction of a Fused Posterior Raster - Benchmark Baseline
#'
#' Collapses `run_fusion()`'s per-pixel posterior to one distribution per mukey by taking
#' `terra::zonal()` means of its percentile layers, in the `low`/`rep`/`high` shape
#' `fuse_observed_data_into_priors()`'s `observed_data_by_mukey` accepts. **Percentile-based only -
#' there is deliberately no point-estimate ("mean") path**, which would narrow the fed-back prior
#' and understate uncertainty (see the design doc's Problem A redraft).
#'
#' Its sole purpose is the A.6 benchmark: the mukey-collapsed comparison arm for
#' `remarginalize_ensemble_to_posterior()`, run through the *same* tabular Monte Carlo machinery
#' via `observed_data_by_mukey`. Kept `@keywords internal` until that benchmark says whether the
#' cheap path is adequate for production (then promote) or not (then delete).
#'
#' @param posterior A `run_fusion()`-style `list(percentiles = , ...)`.
#' @param mukey_raster A categorical mukey `terra::SpatRaster` (any grid - resampled to the
#'   posterior grid via nearest-neighbour internally).
#' @param probs Length-3 ascending probabilities for the `low`/`rep`/`high` triplet; each must be
#'   one of the posterior's own percentile-layer probabilities (no interpolation).
#' @return Named list keyed by mukey code, each element `c(low =, rep =, high =)`.
#' @keywords internal
zonal_distribution_from_posterior <- function(posterior, mukey_raster,
                                              probs = c(0.05, 0.5, 0.95)) {
  if (length(probs) != 3 || is.unsorted(probs)) {
    stop("zonal_distribution_from_posterior(): `probs` must be 3 ascending probabilities.")
  }
  post_probs <- posterior_percentile_probs(posterior)
  if (!all(probs %in% post_probs)) {
    stop(sprintf("zonal_distribution_from_posterior(): every `probs` value must be an exact posterior percentile (have %s).",
                 paste(post_probs, collapse = ", ")))
  }
  layer_for <- function(p) posterior$percentiles[[paste0("P", round(p * 100))]]
  ref <- layer_for(probs[1])

  mk <- terra::resample(mukey_raster, ref, method = "near")
  mk_labels <- terra::cats(mk)[[1]]
  label_col <- names(mk_labels)[2]

  reduce_one <- function(p) {
    z <- terra::zonal(layer_for(p), mk, fun = "mean", na.rm = TRUE)
    stats::setNames(z[[2]], as.character(z[[1]]))
  }
  lo <- reduce_one(probs[1]); rp <- reduce_one(probs[2]); hi <- reduce_one(probs[3])

  mukeys <- Reduce(intersect, list(names(lo), names(rp), names(hi)))
  # zonal() keys by the factor label when the raster is categorical; map back to the mukey code.
  code_of <- function(key) {
    if (!is.null(mk_labels) && key %in% as.character(mk_labels[[1]])) {
      as.character(mk_labels[[label_col]][match(key, as.character(mk_labels[[1]]))])
    } else {
      key
    }
  }
  stats::setNames(
    lapply(mukeys, function(k) c(low = unname(lo[k]), rep = unname(rp[k]), high = unname(hi[k]))),
    vapply(mukeys, code_of, character(1))
  )
}

#' Invert a per-pixel posterior CDF at a per-mukey probability, vectorized over the grid
#'
#' For a scalar-per-mukey probability raster `u_ras`, returns the raster whose every cell is that
#' cell's posterior value at its `u`, by piecewise-linear interpolation on the posterior's
#' percentile knots with `rule = 2` constant extrapolation - the inverse-direction analogue of
#' `sim_linear_cdf_batch()`.
#' @param pct Named list of percentile-value `SpatRaster`s (ascending).
#' @param post_probs Matching ascending probabilities.
#' @param u_ras A `SpatRaster` of probabilities in (0, 1). **May have many layers** - each of the
#'   single-layer `pct` rasters is recycled across them, so one call inverts a whole realization
#'   stack in a fixed `~4 * (length(post_probs) - 1)` terra ops regardless of `nlyr(u_ras)`.
#' @return A `SpatRaster` of interpolated posterior values, `nlyr(u_ras)` layers.
#' @keywords internal
invert_posterior_cdf_raster <- function(pct, post_probs, u_ras) {
  k <- length(post_probs)
  # Work on a placeholder where u is NA (arithmetic on NA in terra doesn't reliably propagate to
  # the final ifel), then restore NA at the end via terra::mask().
  u <- terra::ifel(is.na(u_ras), 0.5, u_ras)
  out <- pct[[1]] + u * 0                           # rule = 2 low clamp, broadcast to nlyr(u_ras)
  for (i in seq_len(k - 1L)) {
    t <- (u - post_probs[i]) / (post_probs[i + 1L] - post_probs[i])
    seg_val <- pct[[i]] * (1 - t) + pct[[i + 1L]] * t
    in_seg <- (u >= post_probs[i]) & (u < post_probs[i + 1L])
    out <- terra::ifel(in_seg, seg_val, out)
  }
  out <- terra::ifel(u >= post_probs[k], pct[[k]], out)
  terra::mask(out, u_ras)                           # NA wherever the input u was NA
}

#' Saxton-Rawls water retention (field capacity, wilting point), raster / stack-native
#'
#' The vectorized, `terra`-native counterpart of `compute_saxton_rawls()` - identical
#' equations and clamps, but every operation is `terra` `Arith`/`Math`/`clamp`/`ifel` so it runs
#' on multi-layer realization stacks in one pass. Used by `remarginalize_awc()`'s default
#' (`method = "saxton_rawls"`) path, because SOLUS100 publishes **no** water-retention variable -
#' `wr_3b`/`wr_15b` can never be fused directly, only derived from the fusable
#' sand/silt/clay/`dbovendry`/`soc`/`fragvol`.
#'
#' Because `remarginalize_awc()` re-marginalizes `sand`/`silt`/`clay` to their own fused
#' posteriors independently, a realization's texture triple may not sum to 100. This
#' reproduces `compute_saxton_rawls()`'s renormalization: where the texture sum is off
#' by more than 5 points, the three fractions are rescaled to sum to 100 before the equations run
#' (the equations themselves use only the sand and clay fractions).
#'
#' @param sand,clay,silt,db,rfv,om Percent (or g/cm^3 for `db`) `SpatRaster`s. Any may be single
#'   or multi-layer; single-layer inputs are recycled.
#' @return `list(fc =, wp =)` - volumetric % `SpatRaster`s, `wp < fc` enforced.
#' @keywords internal
saxton_rawls_raster <- function(sand, clay, silt, db, rfv, om) {
  sc <- terra::clamp(sand, 0, 100)
  cc <- terra::clamp(clay, 0, 100)
  sic <- terra::clamp(silt, 0, 100)
  tot <- sc + cc + sic
  renorm <- abs(tot - 100) > 5                       # matches compute_saxton_rawls()
  s  <- terra::ifel(renorm, sc / tot, sc / 100)
  cl <- terra::ifel(renorm, cc / tot, cc / 100)
  dbc <- terra::clamp(db, 0.6, 2.5)
  rfvf <- terra::clamp(rfv, 0, 95) / 100
  omf <- terra::clamp(om, 0.1, 50) / 100

  # Shared coefficient core - identical arithmetic to compute_saxton_rawls(), one copy.
  core <- .saxton_rawls_gravimetric(s, cl, omf)

  fc_v <- terra::clamp(core$fc_g * dbc * (1 - rfvf) * 100, 3, 65)
  wp_v <- terra::clamp(core$wp_g * dbc * (1 - rfvf) * 100, 1, 45)
  wp_v <- terra::ifel(wp_v >= fc_v, fc_v * 0.6, wp_v)
  list(fc = fc_v, wp = wp_v)
}

#' Per-Pixel Ensemble Re-Marginalization to a Fused Posterior
#'
#' Rank-preserving marginal transform: for each property and depth window, maps every retained
#' source realization through `F_prior` (the mukey ensemble's empirical CDF) then `F_posterior^-1`
#' (the pixel's fused posterior, piecewise-linear on its percentile knots). Each realization keeps
#' its rank position in every property and window, so the ensemble's empirical copula
#' (cross-property, cross-depth rank correlations) carries over; only the marginals become the
#' pixel's posterior.
#'
#' @section Documented approximations:
#' The copula / correlation structure stays mukey-level (no sub-mukey spatially-varying
#' dependence); within-depth-window vertical shape comes from the source realization (window
#' aggregates only - no horizon push-down here); the copula is assumed invariant under
#' re-marginalization and is the SSURGO tabular one (KSSL + joint-copula vertical correlation) -
#' SOLUS carries no joint information, so nothing is re-estimated; the preserved correlation is
#' **rank (Spearman)**, not Pearson; composition (comppct) weighting stays mukey-level (already in
#' the source ensemble's differential replicate counts). Still strictly better than zonal
#' aggregation. See `docs/09_multi_source_raster_fusion_pipeline.md` -> "Statistical structure of
#' the per-pixel bridge" for how this shapes a derived-quantity uncertainty band.
#'
#' @param mukey_ensemble An `extract_mukey_joint_ensemble()` result.
#' @param posterior_by_property_window Named list `[[property]][[window]]` of `run_fusion()`
#'   posteriors (each `list(percentiles = , ...)`), on a common grid. Only properties present in
#'   BOTH this and the ensemble, and windows present in both, are produced.
#' @param n_out Max source realizations to carry per pixel (default 250). If a mukey has more
#'   ensemble replicates, the first `n_out` are used (a fixed subsample - same rows across every
#'   property/window, so the copula is preserved).
#' @param summarize `TRUE` (default) returns per-pixel percentile rasters (`probs`); `FALSE`
#'   returns the raw `n_kept`-layer realization stacks (for downstream use, e.g. `remarginalize_awc()`).
#' @param probs Percentile probabilities for the `summarize = TRUE` output.
#' @return `list(percentiles = [[property]][[window]] = <named list of P.. SpatRasters>)` when
#'   `summarize`, else `list(ensemble = [[property]][[window]] = <n_kept-layer SpatRaster>)`. Plus
#'   `properties`, `window_names`, `n_out`, `n_kept`.
#' @seealso `extract_mukey_joint_ensemble()`, `zonal_distribution_from_posterior()`, `remarginalize_awc()`
#' @family raster-fusion
#' @export
remarginalize_ensemble_to_posterior <- function(mukey_ensemble, posterior_by_property_window,
                                                 n_out = 250, summarize = TRUE,
                                                 probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  ens <- mukey_ensemble$by_mukey
  props <- intersect(mukey_ensemble$properties, names(posterior_by_property_window))
  if (length(props) == 0) {
    stop("remarginalize_ensemble_to_posterior(): no property is present in both the ensemble and posterior_by_property_window.")
  }

  # ONE retained realization count across every mukey AND window, so realization r is the same
  # source draw everywhere - preserves the copula across properties and across depth windows.
  n_kept <- min(
    n_out,
    min(vapply(ens, function(e) min(vapply(e$windows, nrow, integer(1))), integer(1)))
  )
  if (n_kept < 1L) stop("remarginalize_ensemble_to_posterior(): the ensemble has no retained realizations.")

  mukey_codes <- suppressWarnings(as.numeric(names(ens)))
  ok_code <- !is.na(mukey_codes)

  out_pct <- list()
  out_ens <- list()

  for (p in props) {
    windows <- intersect(mukey_ensemble$window_names, names(posterior_by_property_window[[p]]))
    for (w in windows) {
      posterior <- posterior_by_property_window[[p]][[w]]
      post_probs <- posterior_percentile_probs(posterior)
      pct <- posterior$percentiles

      # levels(mk) <- NULL makes the cell values the numeric mukey codes themselves, matching the
      # terra::subst(from = <numeric codes>) broadcast pattern used in R/core-fusion.R.
      mk <- terra::resample(mukey_ensemble$mukey_raster, pct[[1]], method = "near")
      levels(mk) <- NULL

      # [n_mukey x n_kept] rank-probabilities (empirical F_prior), one row per mukey.
      # na.last = "keep": a mukey with a missing value for this property/window (e.g. some cokey's
      # texture sim failed) keeps that realization NA rather than letting base rank()'s default
      # na.last = TRUE rank the NAs 1..n and fabricate a probability for them.
      u_mat <- do.call(rbind, lapply(names(ens), function(mkc) {
        v <- ens[[mkc]]$windows[[w]][seq_len(n_kept), p]
        nf <- sum(is.finite(v))
        if (nf == 0L) return(rep(NA_real_, n_kept))
        (rank(v, ties.method = "average", na.last = "keep") - 0.5) / nf
      }))

      # ONE terra::subst -> an n_kept-layer per-cell u stack; ONE CDF inversion over the 9 knots.
      u_stack <- terra::subst(mk, from = mukey_codes[ok_code],
                              to = u_mat[ok_code, , drop = FALSE], others = NA_real_)
      # terra::subst() with an NA in `to` LEAVES the cell at its raw mukey code rather than NA-ing
      # it; every real rank-probability is strictly in (0, 1), so clamp anything else to NA (mukey
      # codes are >> 1) - otherwise invert_posterior_cdf_raster() rule-2-clamps it to the P99 layer.
      u_stack <- terra::ifel(u_stack > 0 & u_stack < 1, u_stack, NA_real_)
      inv <- invert_posterior_cdf_raster(pct, post_probs, u_stack)
      names(inv) <- paste0("r", seq_len(n_kept))

      if (summarize) {
        qs <- terra::quantile(inv, probs = probs, na.rm = TRUE)
        out_pct[[p]][[w]] <- stats::setNames(
          lapply(seq_along(probs), function(i) qs[[i]]),
          paste0("P", round(probs * 100))
        )
      } else {
        out_ens[[p]][[w]] <- inv
      }
    }
  }

  res <- if (summarize) list(percentiles = out_pct) else list(ensemble = out_ens)
  c(res, list(properties = props, window_names = mukey_ensemble$window_names,
              n_out = n_out, n_kept = n_kept))
}

#' Per-Pixel Available Water Capacity from a Re-Marginalized Ensemble
#'
#' Plant-available water capacity (cm, summed over the ensemble's depth windows), per realization,
#' then summarized per pixel. Deliberately does NOT use `compute_aws()` (ROSETTA + van
#' Genuchten): that needs a live network POST and its own Monte Carlo per call, so it cannot run
#' per (pixel, realization).
#'
#' @section Method:
#' - **`"saxton_rawls"`** (default): `SOLUS100 publishes no water-retention variable`, so `wr_3b`/
#'   `wr_15b` can never be fused directly. Instead the fusable
#'   `sand_total`/`silt_total`/`clay_total`/`db` (plus optional `soc`, `rfv`) are re-marginalized
#'   to their per-pixel fused posteriors, then `saxton_rawls_raster()` derives field capacity /
#'   wilting point per (pixel, realization). `AWC = sum_windows (fc - wp)/100 * thickness` (the
#'   Saxton-Rawls conversion already applies the rock-fragment correction internally).
#' - **`"direct"`**: uses `wr_3b`/`wr_15b` posteriors supplied in `posterior_by_property_window`
#'   directly (`AWC = sum_windows (wr_3b - wr_15b)/100 * thickness * (1 - rfv/100)`). Only usable
#'   if the caller has a water-retention likelihood from some non-SOLUS source.
#'
#' @section Interpreting the per-pixel AWC distribution:
#' `AWC` is computed **per source realization** and combined by `terra::quantile()` per pixel, so
#' `awc_cm$P95 - awc_cm$P5` is a genuine 90% credible interval. Because
#' `remarginalize_ensemble_to_posterior()` threads one realization index through every property and
#' every window (rank-preserving), each realization's inputs are cross-property *and* vertically
#' consistent - a "wet" profile is wet at all depths. What the band does and does not capture:
#' \itemize{
#'   \item **Marginal** uncertainty of each input property: the full SSURGO x SOLUS fused posterior
#'     spread per pixel, propagated through Saxton-Rawls. Captured.
#'   \item **Dependence** (cross-property, cross-depth): the SSURGO tabular copula
#'     (KSSL + joint-copula vertical correlation), *held fixed*. SOLUS carries no joint information,
#'     so nothing is re-estimated. The preserved correlation is **rank (Spearman)**, not Pearson.
#'   \item **Pedotransfer-function error**: NOT propagated. Saxton-Rawls is applied deterministically
#'     per realization; the +/-15\% band `compute_saxton_rawls()` returns is unused.
#'   \item **Sub-mukey variation in dependence**: none. Every pixel in a mukey shares the source
#'     ensemble; only the fused marginals vary pixel-to-pixel.
#' }
#' See `docs/09_multi_source_raster_fusion_pipeline.md` -> "Statistical structure of the per-pixel
#' bridge" for the full treatment.
#'
#' @param mukey_ensemble An `extract_mukey_joint_ensemble()` result.
#' @param posterior_by_property_window As in `remarginalize_ensemble_to_posterior()`, keyed by the
#'   ensemble's property (sim-column) names. `"saxton_rawls"` needs
#'   `sand_total`/`silt_total`/`clay_total`/`db` (+ optional `soc`, `rfv`); `"direct"` needs
#'   `wr_3b`/`wr_15b` (+ optional `rfv`).
#' @param method `"saxton_rawls"` (default) or `"direct"` - see Method.
#' @param n_out,probs As in `remarginalize_ensemble_to_posterior()`.
#' @param rock_fragment If `TRUE` (default) apply the rock-fragment correction when an `rfv`
#'   posterior is available.
#' @param soc_to_om Multiplier converting the `soc` posterior (organic carbon %) to organic matter
#'   % for Saxton-Rawls (default 1.724, the Van Bemmelen factor). Only used by `"saxton_rawls"`.
#' @param tile_rows Optional. Process the AOI in horizontal row-strips, releasing each strip's
#'   realization stack before the next. Per-pixel quantiles are cell-independent, so a tiled run is
#'   numerically identical to `tile_rows = NULL`. `NULL` (default) / `>= nrow` disables tiling.
#' @param restriction_depth Optional single-layer `terra::SpatRaster` of restriction depth in cm -
#'   typically `fetch_solus_restriction_depth(aoi_vect)`,
#'   already right-censoring-guarded (cells with no detected restriction are `Inf`, never
#'   `SOLUS_RESTRICTION_CENSOR_CM` itself). `NULL` (default) preserves today's behavior exactly -
#'   every window contributes its full nominal thickness, regardless of bedrock. When supplied,
#'   resampled once (bilinear - a continuous depth value, not a category) onto the working grid,
#'   and each window's contribution is scaled by its **effective** thickness
#'   `pmax(0, pmin(bottom, restriction_depth) - top)` instead of the nominal `bottom - top`: a
#'   window entirely above the restriction keeps its full thickness, one straddling it gets partial
#'   credit, one entirely below it contributes exactly **0** (not `NA` - `NA` would incorrectly
#'   blank the whole pixel's AWC, including valid shallower windows, rather than just this window's
#'   contribution). A pixel with `restriction_depth` missing/`NA` (e.g. no coverage at the AOI edge)
#'   is treated the same as `Inf` - no truncation - the same "no evidence -> don't truncate"
#'   convention this pipeline already uses for other SSURGO/SOLUS coverage gaps, rather than
#'   propagating `NA` into the AWC total.
#' @return `list(awc_cm = <named list of P.. SpatRasters>, n_kept =, windows_used =, n_tiles =)`.
#'   AWC clamped at 0.
#' @seealso `remarginalize_ensemble_to_posterior()`, `saxton_rawls_raster()`, `compute_aws()`,
#'   `fetch_solus_restriction_depth()`
#' @family raster-fusion
#' @export
remarginalize_awc <- function(mukey_ensemble, posterior_by_property_window,
                               method = c("saxton_rawls", "direct"),
                               n_out = 250, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                               rock_fragment = TRUE, soc_to_om = 1.724, tile_rows = NULL,
                               restriction_depth = NULL) {
  method <- match.arg(method)
  req <- if (method == "saxton_rawls") c("sand_total", "silt_total", "clay_total", "db")
         else c("wr_3b", "wr_15b")
  opt <- if (method == "saxton_rawls") c("soc", "rfv") else c("rfv")

  if (!all(req %in% mukey_ensemble$properties)) {
    stop(sprintf("remarginalize_awc(method = '%s'): the ensemble must carry %s.",
                 method, paste(req, collapse = "/")))
  }
  if (!all(req %in% names(posterior_by_property_window))) {
    stop(sprintf("remarginalize_awc(method = '%s'): posterior_by_property_window needs %s posteriors.",
                 method, paste(req, collapse = "/")))
  }
  keys <- intersect(c(req, opt), names(posterior_by_property_window))
  if (!rock_fragment) keys <- setdiff(keys, "rfv")
  pbpw <- posterior_by_property_window[keys]

  block <- if (method == "saxton_rawls") {
    function(me, p, rd) .awc_block(me, p, n_out, probs, method = "saxton_rawls", soc_to_om = soc_to_om,
                                   restriction_depth = rd)
  } else {
    function(me, p, rd) .awc_block(me, p, n_out, probs, method = "direct", restriction_depth = rd)
  }

  w1 <- intersect(mukey_ensemble$window_names, names(pbpw[[req[1]]]))
  if (length(w1) == 0) stop("remarginalize_awc(): no window has a posterior for the required properties.")
  ref <- pbpw[[req[1]]][[w1[1]]]$percentiles[[1]]

  # S2: resample restriction_depth onto the working grid ONCE here (not per window inside
  # .awc_block(), not per tile below) - a continuous depth value, bilinear, matching how every
  # other continuous SOLUS/SSURGO raster in this pipeline is resampled onto a working grid. NA
  # cells (e.g. no coverage at the AOI edge) are treated as Inf - "no evidence, don't truncate" -
  # rather than propagating NA into the AWC total.
  restriction_depth_aligned <- if (is.null(restriction_depth)) {
    NULL
  } else {
    rd <- terra::resample(restriction_depth, ref, method = "bilinear")
    terra::ifel(is.na(rd), Inf, rd)
  }

  if (is.null(tile_rows) || tile_rows >= terra::nrow(ref)) {
    return(c(block(mukey_ensemble, pbpw, restriction_depth_aligned), list(n_tiles = 1L)))
  }

  nr <- terra::nrow(ref)
  ymax0 <- terra::ymax(ref); yres <- terra::yres(ref)
  x0 <- terra::xmin(ref); x1 <- terra::xmax(ref)
  starts <- seq(1L, nr, by = as.integer(tile_rows))

  blocks <- lapply(starts, function(s) {
    e <- min(s + as.integer(tile_rows) - 1L, nr)
    ext_b <- terra::ext(x0, x1, ymax0 - e * yres, ymax0 - (s - 1L) * yres)
    me_b <- mukey_ensemble
    me_b$mukey_raster <- terra::crop(mukey_ensemble$mukey_raster, ext_b, snap = "out")
    pbpw_b <- lapply(pbpw, function(prop) lapply(prop, function(post) {
      post$percentiles <- lapply(post$percentiles, terra::crop, y = ext_b, snap = "near")
      post
    }))
    rd_b <- if (is.null(restriction_depth_aligned)) NULL else terra::crop(restriction_depth_aligned, ext_b, snap = "near")
    block(me_b, pbpw_b, rd_b)
  })

  merged <- stats::setNames(lapply(seq_along(probs), function(j) {
    do.call(terra::merge, lapply(blocks, function(b) b$awc_cm[[j]]))
  }), paste0("P", round(probs * 100)))

  list(awc_cm = merged, n_kept = blocks[[1]]$n_kept,
       windows_used = blocks[[1]]$windows_used, n_tiles = length(starts))
}

# One-block (whole-grid) AWC core; called directly, or once per row-strip when tile_rows is set.
# restriction_depth (S2): NULL (default) = today's unchanged behavior, every window uses its full
# nominal thickness `th`. Otherwise a single-layer SpatRaster (already resampled onto the working
# grid, already NA->Inf'd, by the caller) - each window's thickness becomes the per-pixel effective
# thickness pmax(0, pmin(bottom, restriction_depth) - top), via two terra::clamp() calls (clamp
# with only `upper` set is exactly pmin(x, upper); with only `lower` set, exactly pmax(x, lower)).
.awc_block <- function(mukey_ensemble, pbpw, n_out, probs, method, soc_to_om = 1.724,
                       restriction_depth = NULL) {
  rm_res <- remarginalize_ensemble_to_posterior(mukey_ensemble, pbpw, n_out = n_out, summarize = FALSE)
  e <- rm_res$ensemble

  awc <- NULL
  used <- character(0)
  for (i in seq_along(mukey_ensemble$window_names)) {
    w <- mukey_ensemble$window_names[i]
    win <- mukey_ensemble$depth_windows[[i]]
    top <- win[1]; bottom <- win[2]
    th <- if (is.null(restriction_depth)) {
      bottom - top
    } else {
      terra::clamp(terra::clamp(restriction_depth, upper = bottom, values = TRUE) - top,
                   lower = 0, values = TRUE)
    }

    if (method == "saxton_rawls") {
      if (any(vapply(c("sand_total", "silt_total", "clay_total", "db"),
                     function(p) is.null(e[[p]]) || is.null(e[[p]][[w]]), logical(1)))) next
      tmpl <- e[["sand_total"]][[w]][[1]]
      om  <- if (!is.null(e[["soc"]]) && !is.null(e[["soc"]][[w]])) e[["soc"]][[w]] * soc_to_om
             else terra::setValues(tmpl, 2)
      rfv <- if (!is.null(e[["rfv"]]) && !is.null(e[["rfv"]][[w]])) e[["rfv"]][[w]]
             else terra::setValues(tmpl, 0)
      sr <- saxton_rawls_raster(e[["sand_total"]][[w]], e[["clay_total"]][[w]],
                                e[["silt_total"]][[w]], e[["db"]][[w]], rfv, om)
      layer <- (sr$fc - sr$wp) / 100 * th
    } else {
      if (is.null(e[["wr_3b"]][[w]]) || is.null(e[["wr_15b"]][[w]])) next
      layer <- (e[["wr_3b"]][[w]] - e[["wr_15b"]][[w]]) / 100 * th
      if (!is.null(e[["rfv"]]) && !is.null(e[["rfv"]][[w]])) layer <- layer * (1 - e[["rfv"]][[w]] / 100)
    }

    awc <- if (is.null(awc)) layer else awc + layer
    used <- c(used, w)
  }
  if (is.null(awc)) stop("remarginalize_awc(): no usable window.")

  awc <- terra::clamp(awc, lower = 0, values = TRUE)
  qs <- terra::quantile(awc, probs = probs, na.rm = TRUE)
  list(
    awc_cm = stats::setNames(lapply(seq_along(probs), function(j) qs[[j]]), paste0("P", round(probs * 100))),
    n_kept = rm_res$n_kept,
    windows_used = used
  )
}
