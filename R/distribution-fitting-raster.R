#' @title Raster-Native Percentile-Based Distribution Fitting
#'
#' @description Fits a distribution to a property's percentile-value rasters
#'   (e.g. SSURGO low/rep/high, SOLUS low/pred/high) as pure `terra::SpatRaster`
#'   arithmetic and/or `terra::app()`/`terra::lapp()` calls to *vectorized*
#'   base-R math functions (`qnorm`, `qbeta`, `digamma`/`trigamma`) - never a
#'   per-cell optimizer, which is what makes these fast across whole rasters.
#'   Ported from `code_ref/reanalysis-platform/distribution_fitting_raster.R`.
#'
#'   Reuses `R/distributions.R`'s existing `metalog_basis_matrix()`/
#'   `metalog_to_z()`/`metalog_from_z()` directly rather than duplicating
#'   them: those functions are pure elementwise arithmetic (`log`/`exp`/`/`
#'   plus building a small, non-spatial coefficient matrix from the fixed
#'   probability grid), so they already work unchanged whether their `x`/`z`
#'   arguments are plain numerics or `terra::SpatRaster`s - confirmed by a
#'   dedicated smoke test before this file was written, not assumed.
#'
#'   Status, from the original source:
#'   - `fit_normal_raster()` / `quantile_normal_raster()` - EXACT.
#'   - `quantile_linear_cdf_raster()` - EXACT (nonparametric fallback target).
#'   - `fit_beta_mom_raster()` - method-of-moments Beta fit, APPROXIMATE
#'     relative to MLE, used only as the Newton-Raphson seed below.
#'   - `fit_beta_mle_newton_raster()` / `quantile_beta_mle_newton_raster()` -
#'     vectorized Newton-Raphson Beta MLE, validated (upstream) directly
#'     against `fitdistrplus::fitdist()` (~1e-3 to 1e-5).
#'   - `fit_metalog_linear_raster()` / `quantile_metalog_linear_raster()` /
#'     `check_metalog_feasibility_raster()` / `quantile_metalog_linear_with_fallback()` -
#'     exact linear-solve metalog reformulation (matches `rmetalog::metalog()`
#'     to ~1e-12 when its own fit is feasible, per upstream validation), with
#'     a validated fallback to `linear_cdf` for the cells where it isn't.
#'
#'   `spline` is intentionally excluded (per the original project's own
#'   decision - its accuracy relative to R's `splinefun(method="monoH.FC")`
#'   was never confirmed bit-exact).
#' @name distribution_fitting_raster
NULL

# ---------------------------------------------------------------------------
# EXACT: normal (raster arithmetic)
# ---------------------------------------------------------------------------

#' Fit a Normal distribution's mu/sigma rasters from three percentile-value
#' rasters, then evaluate at a fixed target quantile `q` - pure raster
#' arithmetic throughout.
#' @param p_lo_r,p50_r,p_hi_r SpatRasters for the low/median/high percentile values.
#' @param p_lo,p_hi The probabilities `p_lo_r`/`p_hi_r` represent (e.g. 0.025/0.975).
#' @return list(mu = SpatRaster, sigma = SpatRaster)
#' @export
fit_normal_raster <- function(p_lo_r, p50_r, p_hi_r, p_lo, p_hi) {
  list(mu = p50_r, sigma = (p_hi_r - p_lo_r) / (stats::qnorm(p_hi) - stats::qnorm(p_lo)))
}
#' @rdname fit_normal_raster
#' @param fit Output of `fit_normal_raster()`.
#' @param q Target quantile probability.
#' @export
quantile_normal_raster <- function(fit, q) {
  fit$mu + fit$sigma * stats::qnorm(q)
}

# ---------------------------------------------------------------------------
# EXACT: linear_cdf (fixed-breakpoint interpolation weight, pure raster arithmetic)
# ---------------------------------------------------------------------------

#' Evaluate the piecewise-linear inverse-CDF at a fixed quantile `q` across a
#' raster. Because the percentile breakpoints (`probs`) are the same at
#' every cell, only the two value-rasters bracketing `q` and a single scalar
#' interpolation weight are needed - no per-cell branching required since
#' the bracket is resolved once, outside the raster arithmetic.
#' @param value_rasters List of SpatRasters, sorted the same way as `probs`.
#' @param probs Numeric probabilities matching `value_rasters`' order.
#' @param q Target quantile probability.
#' @export
quantile_linear_cdf_raster <- function(value_rasters, probs, q) {
  stopifnot(q >= min(probs), q <= max(probs))
  i <- max(which(probs <= q))
  if (i == length(probs)) i <- length(probs) - 1
  p_lo <- probs[i]; p_hi <- probs[i + 1]
  weight <- if (p_hi == p_lo) 0 else (q - p_lo) / (p_hi - p_lo)
  value_rasters[[i]] + weight * (value_rasters[[i + 1]] - value_rasters[[i]])
}

# ---------------------------------------------------------------------------
# APPROXIMATE (method-of-moments): beta - used as the Newton-Raphson seed
# ---------------------------------------------------------------------------

#' Method-of-moments Beta fit: alpha/beta rasters are pure arithmetic on
#' mean/variance rasters.
#' @param mean_r,var_r SpatRasters of the rescaled-to-unit-interval mean/variance.
#' @export
fit_beta_mom_raster <- function(mean_r, var_r) {
  common <- mean_r * (1 - mean_r) / var_r - 1
  list(alpha = mean_r * common, beta = (1 - mean_r) * common)
}
#' @rdname fit_beta_mom_raster
#' @param fit Output of `fit_beta_mom_raster()`.
#' @param q Target quantile probability.
#' @export
quantile_beta_mom_raster <- function(fit, q) {
  terra::lapp(c(fit$alpha, fit$beta), fun = function(a, b) stats::qbeta(q, a, b))
}

# ---------------------------------------------------------------------------
# Vectorized Newton-Raphson Beta MLE
# ---------------------------------------------------------------------------

#' `fitdistrplus::fitdist()`'s generic `optim()`-based MLE has no vectorized
#' form. The Beta log-likelihood's score/Hessian equations have a known
#' closed form in terms of `digamma()`/`trigamma()` (both vectorized base-R
#' functions), so a fixed-iteration-count Newton-Raphson solver updates ALL
#' cells' `(alpha, beta)` estimates simultaneously per iteration.
#' @keywords internal
raster_digamma <- function(r) terra::app(r, fun = digamma)
#' @keywords internal
raster_trigamma <- function(r) terra::app(r, fun = trigamma)

#' Fit a Beta distribution's (alpha, beta) rasters via vectorized
#' Newton-Raphson MLE.
#'
#' If a percentile value lands exactly on `bounds` (rescaled x = 0 or 1
#' exactly), `log(x)`/`log(1-x)` in the score equations would be `-Inf`,
#' degenerating the fit to the clamp floor for every cell uniformly -
#' unlikely with real empirical percentiles but a real risk for synthetic or
#' explicitly bound-recorded data. Clamped away from the exact boundary by a
#' tiny epsilon - a no-op for realistic input.
#' @param value_rasters List of percentile-value SpatRasters (any count >= 3).
#' @param bounds c(lower, upper) physical bounds.
#' @param n_iter Fixed Newton-Raphson iteration count (validated sufficient at 15 upstream).
#' @param eps Clamp epsilon away from the exact 0/1 boundary.
#' @return list(alpha = SpatRaster, beta = SpatRaster)
#' @export
fit_beta_mle_newton_raster <- function(value_rasters, bounds, n_iter = 15, eps = 1e-6) {
  k <- length(value_rasters)
  scaled <- lapply(value_rasters, function(r) terra::clamp((r - bounds[1]) / (bounds[2] - bounds[1]), lower = eps, upper = 1 - eps, values = TRUE))
  g1 <- Reduce(`+`, lapply(scaled, function(r) log(r))) / k
  g2 <- Reduce(`+`, lapply(scaled, function(r) log(1 - r))) / k

  mean_r <- Reduce(`+`, scaled) / k
  var_r <- Reduce(`+`, lapply(scaled, function(r) (r - mean_r)^2)) / k
  mom <- fit_beta_mom_raster(mean_r, var_r)
  alpha_r <- terra::ifel(mom$alpha < 0.1, 0.1, mom$alpha)
  beta_r <- terra::ifel(mom$beta < 0.1, 0.1, mom$beta)

  for (i in seq_len(n_iter)) {
    a_plus_b <- alpha_r + beta_r
    dg_apb <- raster_digamma(a_plus_b); dg_a <- raster_digamma(alpha_r); dg_b <- raster_digamma(beta_r)
    tg_apb <- raster_trigamma(a_plus_b); tg_a <- raster_trigamma(alpha_r); tg_b <- raster_trigamma(beta_r)

    S1 <- dg_apb - dg_a + g1
    S2 <- dg_apb - dg_b + g2
    H11 <- tg_apb - tg_a; H12 <- tg_apb; H21 <- tg_apb; H22 <- tg_apb - tg_b
    det <- H11 * H22 - H12 * H21

    d_alpha <- (H22 * S1 - H12 * S2) / det
    d_beta <- (-H21 * S1 + H11 * S2) / det

    alpha_r <- terra::ifel(alpha_r - d_alpha < 0.01, 0.01, alpha_r - d_alpha)
    beta_r <- terra::ifel(beta_r - d_beta < 0.01, 0.01, beta_r - d_beta)
  }
  list(alpha = alpha_r, beta = beta_r)
}
#' @rdname fit_beta_mle_newton_raster
#' @param fit Output of `fit_beta_mle_newton_raster()`.
#' @param q Target quantile probability.
#' @export
quantile_beta_mle_newton_raster <- function(fit, q) {
  terra::lapp(c(fit$alpha, fit$beta), fun = function(a, b) stats::qbeta(q, a, b))
}

# ---------------------------------------------------------------------------
# EXACT (when feasible): metalog, linear-solve reformulation
# ---------------------------------------------------------------------------

#' Fit a metalog distribution via linear solve, vectorized as raster arithmetic.
#'
#' For a standard percentile setup where the number of interior percentiles
#' (excluding p=0/p=1, where logit is undefined) equals the number of
#' metalog terms, the fit is an EXACTLY-DETERMINED linear system
#' (`solve(Y, z)`), not an optimization - same closed-form math as
#' `R/distributions.R`'s `fit_metalog_linear()`, vectorized here as raster
#' arithmetic: `Y`/`Y_inv` are computed once (fixed, non-spatial probability
#' grid, via the existing `metalog_basis_matrix()`), then applied to
#' different right-hand-side rasters per coefficient.
#'
#' CAVEAT (from upstream validation): this is the same fast path `rmetalog`
#' itself takes when its solution is already feasible (implied density
#' non-negative); `rmetalog` also has an LP-based feasibility-correction
#' fallback for when it isn't - not reproduced here. Use
#' `check_metalog_feasibility_raster()` + `quantile_metalog_linear_with_fallback()`
#' below rather than trusting this fit unconditionally.
#' @param value_rasters List of INTERIOR percentile-value rasters (excluding
#'   p=0/p=1 if present).
#' @param probs Matching interior probabilities.
#' @param bounds `c(lower, upper)`; required unless `boundedness = "u"`.
#' @param boundedness One of `"u"`/`"sl"`/`"su"`/`"b"` - see `metalog_to_z()`.
#' @export
fit_metalog_linear_raster <- function(value_rasters, probs, bounds, boundedness) {
  term <- length(probs)
  Y <- metalog_basis_matrix(probs, term)
  Y_inv <- solve(Y)
  z_rasters <- lapply(value_rasters, metalog_to_z, bounds = bounds, boundedness = boundedness)
  a_rasters <- lapply(seq_len(term), function(j) {
    Reduce(`+`, lapply(seq_len(term), function(i) Y_inv[j, i] * z_rasters[[i]]))
  })
  list(a = a_rasters, term = term)
}
#' @rdname fit_metalog_linear_raster
#' @param fit Output of `fit_metalog_linear_raster()`.
#' @param q Target quantile probability.
#' @export
quantile_metalog_linear_raster <- function(fit, q, bounds, boundedness) {
  Y_q <- metalog_basis_matrix(q, fit$term)[1, ]
  z_q <- Reduce(`+`, lapply(seq_len(fit$term), function(j) Y_q[j] * fit$a[[j]]))
  metalog_from_z(z_q, bounds, boundedness)
}

#' Vectorized feasibility check: a valid metalog quantile function must be
#' monotonically increasing in y (equivalent to its density staying
#' non-negative everywhere). Probes `quantile_metalog_linear_raster()` at a
#' fixed grid of y-values and flags cells where consecutive probe values
#' decrease - a probe, not a proof, but fully vectorized raster arithmetic.
#' Streams one probe raster at a time rather than materializing all of them
#' (validated upstream to avoid an allocation failure at large cell counts).
#' @param fit Output of `fit_metalog_linear_raster()`.
#' @param bounds,boundedness Same as `fit_metalog_linear_raster()`'s arguments for this `fit`.
#' @param y_grid Probability grid to probe.
#' @export
check_metalog_feasibility_raster <- function(fit, bounds, boundedness, y_grid = seq(0.02, 0.98, by = 0.02)) {
  prev_probe <- quantile_metalog_linear_raster(fit, y_grid[1], bounds, boundedness)
  is_infeasible <- prev_probe < prev_probe # all-FALSE raster, same extent/dims/type
  for (y in y_grid[-1]) {
    curr_probe <- quantile_metalog_linear_raster(fit, y, bounds, boundedness)
    is_infeasible <- is_infeasible | (curr_probe < prev_probe)
    prev_probe <- curr_probe
  }
  is_infeasible
}

#' Metalog quantile with automatic fallback to `linear_cdf` for infeasible
#' cells - validated upstream: zero effect on feasible cells, exact
#' `linear_cdf` match on infeasible ones.
#'
#' Feasibility is computed ONCE per fit (via `check_metalog_feasibility_raster()`)
#' and passed in rather than recomputed per call, since it depends only on
#' the fitted coefficients, not on which quantile `q` is being evaluated.
#'
#' @param fit Output of `fit_metalog_linear_raster()` (fit on INTERIOR percentiles).
#' @param infeasible_r Output of `check_metalog_feasibility_raster()` for this same `fit`.
#' @param full_value_rasters,full_probs The FULL percentile set (including
#'   p=0/p=1 if available) - `linear_cdf` has no restriction there and
#'   benefits from every available knot.
#' @param q Target quantile probability.
#' @param bounds,boundedness Same as `fit_metalog_linear_raster()`'s arguments.
#' @return list(value = blended quantile raster, used_fallback = infeasible_r verbatim)
#' @export
quantile_metalog_linear_with_fallback <- function(fit, infeasible_r, full_value_rasters, full_probs,
                                                   q, bounds, boundedness) {
  metalog_q <- quantile_metalog_linear_raster(fit, q, bounds, boundedness)
  fallback_q <- quantile_linear_cdf_raster(full_value_rasters, full_probs, q)
  list(value = terra::ifel(infeasible_r, fallback_q, metalog_q), used_fallback = infeasible_r)
}
