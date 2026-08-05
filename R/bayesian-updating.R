#' @title Bayesian Updating: Conjugate Fusion and General Grid-KDE Fusion
#'
#' @description Standalone toolkit for fusing a prior belief distribution with
#'   a likelihood (e.g. an SSURGO-derived prior against field-measured data),
#'   ported from `code_ref/brdf/bayesian_updating.R` and the validated
#'   `code_ref/reanalysis-platform/{bayes_fuse.R, reanalysis-platform-fusion.R,
#'   texture_ilr_fusion.R}` bundle (there raster-native; here plain-vector/
#'   scalar, since the underlying math is ordinary base-R arithmetic once
#'   `terra::` is stripped out).
#'
#'   **This module is intentionally NOT called anywhere in `monte-carlo.R`.**
#'   It is a standalone, independently-tested set of building blocks for a
#'   future "update an SSURGO-derived prior against observed/field data"
#'   feature - not wired into the main Monte Carlo simulation pipeline in
#'   this pass.
#'
#'   Three tiers, in increasing order of generality/cost:
#'   1. `bayes_update_normal_normal()` - exact, closed-form, Normal-only.
#'   2. `fuse_beta()`/`fuse_gamma()` - exact, closed-form, same-family-only
#'      (each family's kernel parameters simply ADD - the same principle
#'      `bayes_update_normal_normal()` uses for Normal's natural parameters,
#'      applied to Beta/Gamma), with a `feasible` flag and a documented
#'      moment-based fallback route when infeasible.
#'   3. `bayesian_update()` - fully general grid-KDE Bayes' rule (any shape,
#'      any family, even mismatched sides), at the cost of needing raw
#'      samples (not distributional parameters) on both sides.
#'
#'   `fuse_bivariate_normal()`/`fuse_texture_group_from_triplets()` extend
#'   tier 1 to the joint 3-part-composition ILR case (clay/sand/silt texture,
#'   in whichever role order `monte-carlo.R`'s `composition_groups$texture$members`
#'   configures - see `distributions.R`'s ILR section header), since
#'   independent per-fraction fusion (e.g. three separate `fuse_beta()`
#'   calls) measurably breaks sum-to-100 (documented upstream: up to 10.5
#'   percentage points on realistic synthetic data).
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
#' @export
bayes_update_normal_normal <- function(prior_mu, prior_sigma, lik_mu, lik_sigma) {
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
#' impossible. `bayes_update_normal_normal()` is reused in log-space instead.
#'
#' @param mu,sigma Raw-space mean/sd.
#' @return `list(mu = log-space mu, sigma = log-space sigma)`.
#' @export
normal_to_lognormal_params <- function(mu, sigma) {
  sigma_log <- sqrt(log(1 + (sigma / mu)^2))
  mu_log <- log(mu) - sigma_log^2 / 2
  list(mu = mu_log, sigma = sigma_log)
}

#' Inverse of `normal_to_lognormal_params()`
#'
#' Converts log-space Normal(mu, sigma) parameters back to the raw-space
#' Lognormal's mean/sd.
#' @param mu_log,sigma_log Log-space mean/sd.
#' @return `list(mu = raw-space mean, sigma = raw-space sd)`.
#' @export
lognormal_to_normal_params <- function(mu_log, sigma_log) {
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
#' @export
moments_to_gamma <- function(mean, var) {
  list(shape = mean^2 / var, rate = mean / var)
}

#' Method-of-moments Beta fit
#'
#' Used by the closed-form same-family fusion route's infeasible-cell
#' fallback (re-express as Normal moments, fuse as Normal, convert back).
#' @param mean,var Mean/variance to convert.
#' @return `list(alpha=, beta=)`.
#' @export
moments_to_beta <- function(mean, var) {
  common <- mean * (1 - mean) / var - 1
  list(alpha = mean * common, beta = (1 - mean) * common)
}

#' Beta's mean/variance as a function of its own (alpha, beta)
#'
#' The inverse direction of `moments_to_beta()` - used to re-express an
#' infeasible same-family fusion's inputs as Normal moments before falling
#' back to `bayes_update_normal_normal()`.
#' @param alpha,beta Beta shape parameters.
#' @return `list(mean=, var=)`.
#' @export
beta_to_moments <- function(alpha, beta) {
  list(mean = alpha / (alpha + beta), var = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
}

#' Gamma's mean/variance as a function of its own (shape, rate)
#'
#' The inverse direction of `moments_to_gamma()` - used to re-express an
#' infeasible same-family fusion's inputs as Normal moments before falling
#' back to `bayes_update_normal_normal()`.
#' @param shape,rate Gamma shape/rate parameters.
#' @return `list(mean=, var=)`.
#' @export
gamma_to_moments <- function(shape, rate) {
  list(mean = shape / rate, var = shape / rate^2)
}

#' Fuse a prior and likelihood belief distribution of the same family
#'
#' @param prior_params,lik_params Named lists of that family's parameters:
#'   `list(mu=,sigma=)` for `"normal"`, `list(alpha=,beta=)` for `"beta"`,
#'   `list(shape=,rate=)` for `"gamma"`.
#' @param family One of `"normal"`, `"beta"`, `"gamma"`. Both sides must
#'   already be fit in this same family - there is no closed form for fusing
#'   mismatched families; use `bayesian_update()` for that.
#' @return The matching `fuse_*()` function's return value.
#' @export
bayes_fuse <- function(prior_params, lik_params, family = c("normal", "beta", "gamma")) {
  family <- match.arg(family)
  switch(family,
    normal = bayes_update_normal_normal(prior_params$mu, prior_params$sigma, lik_params$mu, lik_params$sigma),
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
#' @return A numeric vector of `n` samples drawn from the posterior distribution.
#' @export
bayesian_update <- function(prior_distribution, likelihood_distribution, grid_range = NULL,
                             grid_resolution = 0.01, n = 1000) {
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

  sample(value_grid, size = n, prob = posterior_prob, replace = TRUE)
}

# ==============================================================================
# Joint (compositional) fusion for texture (clay/sand/silt)
# ==============================================================================

#' Fuse two independent bivariate Normal beliefs about the same 2D quantity
#'
#' Precision MATRICES add - the direct multivariate generalization of
#' `bayes_update_normal_normal()`'s scalar-precision addition.
#'
#' @param mu1,mu2 Length-2 mean vectors.
#' @param Sigma1,Sigma2 2x2 covariance matrices.
#' @return `list(mu = fused mean vector, Sigma = fused covariance matrix)`.
#' @export
fuse_bivariate_normal <- function(mu1, Sigma1, mu2, Sigma2) {
  P1 <- solve(Sigma1)
  P2 <- solve(Sigma2)
  Sigma_post <- solve(P1 + P2)
  list(mu = as.numeric(Sigma_post %*% (P1 %*% mu1 + P2 %*% mu2)), Sigma = Sigma_post)
}

#' Fuse two full clay/sand/silt low-rep-high triplets jointly via ILR fusion
#'
#' Wraps `distributions.R`'s `estimate_ilr_moments_mc()` (once per side) +
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
#'   ILR sequential binary partition - see `distributions.R`'s ILR section
#'   header), not an identity requirement; callers (e.g.
#'   `monte-carlo.R`'s `fuse_observed_data_into_priors()`) build them from
#'   `composition_groups$texture$members` in configured order.
#' @param z_prior,z_lik Standard-normal quantile matching each side's
#'   low/high interval (e.g. `qnorm(0.95)` for a 5th/95th-percentile
#'   low/high).
#' @param n_mc Monte Carlo sample size passed to `estimate_ilr_moments_mc()`.
#' @param n_samples Number of posterior composition samples to draw.
#' @param total Composition target sum (default 100).
#' @return `list(posterior_samples = <n_samples x 3 matrix, columns clay/sand/silt,
#'   guaranteed sum-to-`total` by construction>, ilr_mu=, ilr_Sigma=)`.
#' @export
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
#' general `bayesian_update()`; named lists of family-native parameters
#' route to the closed-form `bayes_fuse()`. This removes the "which
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
#' @param n_samples,grid_resolution Passed to `bayesian_update()` for the
#'   general route.
#' @return For the general route: a numeric vector of posterior samples
#'   (`bayesian_update()`'s native output). For the closed-form route: the
#'   family-native posterior parameter list (`bayes_fuse()`'s native output).
#'   These are NOT the same shape - documented deliberately rather than
#'   forced into a fake-uniform contract, since the caller already knows
#'   which shape it's passing in and therefore which shape it gets back.
#' @export
fuse_property <- function(prior, likelihood, family = NULL, bounds = NULL, method = NULL,
                           n_samples = 1000, grid_resolution = 0.01) {
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
    return(bayesian_update(prior, likelihood, grid_resolution = grid_resolution, n = n_samples))
  }

  if (is.null(family)) {
    stop("fuse_property(): family is required when prior/likelihood are parameter lists (closed_form route).")
  }
  bayes_fuse(prior, likelihood, family = family)
}
