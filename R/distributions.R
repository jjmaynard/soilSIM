#' @title Percentile-Triplet Distribution Fitting, ILR Composition Utilities, and
#'   Matrix Utilities
#'
#' @description Shared foundation used by both `statistics.R` and
#'   `monte-carlo.R` (and, for the ILR pieces, `bayesian-updating.R`).
#'
#'   The percentile-triplet fitting/quantile functions below are adapted from
#'   validated closed-form math originally written as `terra::SpatRaster`
#'   arithmetic (`code_ref/reanalysis-platform/distribution_fitting_raster.R`):
#'   stripped of all `terra::` calls, the underlying math is ordinary base-R
#'   vector arithmetic, so it ports directly. The Beta fit (vectorized
#'   Newton-Raphson MLE) was validated there against `fitdistrplus::fitdist()`
#'   to ~1e-3 to 1e-5; the metalog fit (exact linear solve when the number of
#'   interior percentiles equals the number of terms) was validated to ~1e-12
#'   against `rmetalog::metalog()` - this project therefore never needs the
#'   `rmetalog` package (which `code_ref/brdf/distribution_fitting.R`
#'   documents as prone to hangs and a reproduced segfault).
#'
#'   The ILR (isometric log-ratio) transform is ported verbatim from
#'   `code_ref/reanalysis-platform/texture_ilr_fusion.R`, itself validated to
#'   match `compositions::ilr()` exactly - so this package needs no
#'   `compositions` dependency either.
#' @name distributions
NULL

# ==============================================================================
# 1. PERCENTILE-TRIPLET FITTING / QUANTILE DISPATCH
# ==============================================================================

#' Fit a Normal distribution from three percentiles (exact closed form)
#'
#' @param p_lo_val,p50_val,p_hi_val Values at probabilities `p_lo`, 0.5, `p_hi`.
#' @param p_lo,p_hi The low/high probabilities `p_lo_val`/`p_hi_val` represent.
#' @return `list(mean=, sd=)`.
#' @export
fit_normal_triplet <- function(p_lo_val, p50_val, p_hi_val, p_lo, p_hi) {
  list(
    mean = p50_val,
    sd = (p_hi_val - p_lo_val) / (qnorm(p_hi) - qnorm(p_lo))
  )
}

#' @rdname fit_normal_triplet
#' @param fit A fit from `fit_normal_triplet()`.
#' @param q Vector of probabilities to evaluate the quantile function at.
#' @export
quantile_normal <- function(fit, q) {
  fit$mean + fit$sd * qnorm(q)
}

#' Method-of-moments Beta fit (closed form; used as the Newton-Raphson seed)
#'
#' @param mean,var Mean/variance on the rescaled-to-\eqn{[0,1]} support.
#' @return `list(shape1=, shape2=)`.
#' @export
fit_beta_mom <- function(mean, var) {
  common <- mean * (1 - mean) / var - 1
  list(shape1 = mean * common, shape2 = (1 - mean) * common)
}

#' Vectorized Newton-Raphson Beta MLE from raw-scale percentile values
#'
#' Fits a single `(shape1, shape2)` pair from `k` percentile values of one
#' property at one horizon. The Beta log-likelihood's score/Hessian have a
#' closed form in `digamma()`/`trigamma()` (standard Beta MLE theory),
#' iterated for a fixed number of steps rather than via a general-purpose
#' optimizer - validated (upstream) against `fitdistrplus::fitdist()` to
#' ~1e-3-1e-5 after 15 iterations.
#'
#' @param values Numeric vector of `k` percentile values (raw units).
#' @param bounds `c(lower, upper)` physical bounds.
#' @param n_iter Fixed Newton-Raphson iteration count.
#' @param eps Clamp distance from the exact 0/1 rescaled boundary (avoids
#'   `-Inf` in the score equations when a percentile lands exactly on a bound).
#' @return `list(shape1=, shape2=)`.
#' @export
fit_beta_mle_newton <- function(values, bounds, n_iter = 15, eps = 1e-6) {
  scaled <- pmin(pmax((values - bounds[1]) / (bounds[2] - bounds[1]), eps), 1 - eps)

  g1 <- mean(log(scaled))
  g2 <- mean(log(1 - scaled))

  mean_s <- mean(scaled)
  var_s <- mean((scaled - mean_s)^2)
  mom <- fit_beta_mom(mean_s, var_s)
  shape1 <- max(mom$shape1, 0.1)
  shape2 <- max(mom$shape2, 0.1)

  for (i in seq_len(n_iter)) {
    a_plus_b <- shape1 + shape2
    dg_apb <- digamma(a_plus_b); dg_a <- digamma(shape1); dg_b <- digamma(shape2)
    tg_apb <- trigamma(a_plus_b); tg_a <- trigamma(shape1); tg_b <- trigamma(shape2)

    S1 <- dg_apb - dg_a + g1
    S2 <- dg_apb - dg_b + g2
    H11 <- tg_apb - tg_a; H12 <- tg_apb; H22 <- tg_apb - tg_b
    det <- H11 * H22 - H12 * H12

    d_shape1 <- (H22 * S1 - H12 * S2) / det
    d_shape2 <- (-H12 * S1 + H11 * S2) / det

    shape1 <- max(shape1 - d_shape1, 0.01)
    shape2 <- max(shape2 - d_shape2, 0.01)
  }

  list(shape1 = shape1, shape2 = shape2)
}

#' @rdname fit_beta_mle_newton
#'
#' @description `fit_beta_mle_newton_vec()` runs the identical update across
#'   every row of `value_matrix` simultaneously (`digamma()`/`trigamma()` are
#'   already vectorized base-R functions), so it is the version
#'   `prepare_simulation_parameters()` should call once per property across
#'   all horizons, rather than looping `fit_beta_mle_newton()` per row.
#'
#' @param value_matrix Numeric matrix, one row per horizon, one column per
#'   percentile value (e.g. 3 columns for an l/r/h triplet).
#' @return A data frame with columns `shape1`, `shape2` (one row per horizon).
#' @export
fit_beta_mle_newton_vec <- function(value_matrix, bounds, n_iter = 15, eps = 1e-6) {
  scaled <- pmin(pmax((value_matrix - bounds[1]) / (bounds[2] - bounds[1]), eps), 1 - eps)

  g1 <- rowMeans(log(scaled))
  g2 <- rowMeans(log(1 - scaled))

  mean_s <- rowMeans(scaled)
  var_s <- rowMeans((scaled - mean_s)^2)
  common <- mean_s * (1 - mean_s) / var_s - 1
  shape1 <- pmax(mean_s * common, 0.1)
  shape2 <- pmax((1 - mean_s) * common, 0.1)

  for (i in seq_len(n_iter)) {
    a_plus_b <- shape1 + shape2
    dg_apb <- digamma(a_plus_b); dg_a <- digamma(shape1); dg_b <- digamma(shape2)
    tg_apb <- trigamma(a_plus_b); tg_a <- trigamma(shape1); tg_b <- trigamma(shape2)

    S1 <- dg_apb - dg_a + g1
    S2 <- dg_apb - dg_b + g2
    H11 <- tg_apb - tg_a; H12 <- tg_apb; H22 <- tg_apb - tg_b
    det <- H11 * H22 - H12 * H12

    d_shape1 <- (H22 * S1 - H12 * S2) / det
    d_shape2 <- (-H12 * S1 + H11 * S2) / det

    shape1 <- pmax(shape1 - d_shape1, 0.01)
    shape2 <- pmax(shape2 - d_shape2, 0.01)
  }

  data.frame(shape1 = shape1, shape2 = shape2)
}

#' Evaluate a fitted Beta quantile function on \eqn{[0,1]}
#' @param fit A fit with `shape1`/`shape2` (e.g. from `fit_beta_mle_newton()`).
#' @param q Vector of probabilities.
#' @export
quantile_beta <- function(fit, q) {
  qbeta(q, fit$shape1, fit$shape2)
}

#' Piecewise-linear inverse-CDF quantile function, exact at the given knots
#'
#' Tails are clamped to the min/max value (`rule = 2`) rather than
#' extrapolated, matching `code_ref/brdf/distribution_fitting.R`'s
#' `sim_linear_cdf()`.
#'
#' @param probs,values Matching, sorted percentile probabilities/values.
#' @param q Vector of probabilities to evaluate at.
#' @export
quantile_linear_cdf <- function(probs, values, q) {
  o <- order(probs)
  approxfun(probs[o], values[o], method = "linear", rule = 2)(q)
}

#' Degenerate quantile function for a min/mode/max triangular fit
#' @param fit `list(min=, mode=, max=)`.
#' @param q Vector of probabilities.
#' @export
quantile_triangular <- function(fit, q) {
  min_v <- fit$min; mode_v <- fit$mode; max_v <- fit$max
  if (!is.finite(max_v - min_v) || max_v <= min_v) {
    return(rep(mode_v %||% min_v, length(q)))
  }
  fc <- (mode_v - min_v) / (max_v - min_v)
  ifelse(q < fc,
         min_v + sqrt(q * (max_v - min_v) * (mode_v - min_v)),
         max_v - sqrt((1 - q) * (max_v - min_v) * (max_v - mode_v)))
}

#' Draw random samples from a triangular distribution
#'
#' Ported verbatim from `code_ref/brdf/property_simulation.R`'s `tri_dist()`
#' (originally adapted from the `triangle` package), used by
#' `R/depth-simulation.R`'s profile-depth simulators. Kept as its own
#' random-draw implementation rather than layered on `quantile_triangular()`
#' (a deterministic inverse-CDF evaluator with different degenerate-input
#' handling) so its exact edge-case behavior is preserved: it errors on
#' invalid `n`, but returns `NaN` (not an error) when the mode falls outside
#' `[a, b]` or any parameter is infinite/`NA`.
#'
#' @param n Number of samples to draw. If a vector is passed, its length is
#'   used (matching `stats::rnorm()`-style conventions).
#' @param a,b,c Minimum, maximum, and mode of the triangular distribution.
#' @return A numeric vector of length `n`; `NaN` for every element if `a`,
#'   `b`, `c` describe a degenerate/invalid triangle.
#' @export
tri_dist <- function(n = 1, a = 0, b = 1, c = (a + b) / 2) {
  if (length(n) > 1) n <- length(n)
  if (n < 1 | is.na(n)) stop(paste("invalid argument: n =", n))
  n <- floor(n)

  if (any(is.na(c(a, b, c)))) return(rep(NaN, times = n))
  if (a > c | b < c) return(rep(NaN, times = n))
  if (any(is.infinite(c(a, b, c)))) return(rep(NaN, times = n))

  p <- stats::runif(n)

  if (a != c) {
    i <- which((a + sqrt(p * (b - a) * (c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) > c)
  } else {
    i <- which((a + sqrt(p * (b - a) * (c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >= c)
  }

  if (length(i) != 0) p[i] <- a + sqrt(p[i] * (b - a) * (c - a))
  if (length(j) != 0) p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))

  p
}

# ---------------------------------------------------------------------------
# Metalog: exact linear-solve reformulation (no rmetalog dependency)
# ---------------------------------------------------------------------------

#' Metalog basis matrix
#'
#' Column `j` is the j-th metalog basis function evaluated at each
#' probability in `y`.
#' @param y Numeric vector of probabilities (0,1).
#' @param term Number of metalog terms.
#' @keywords internal
metalog_basis_matrix <- function(y, term) {
  Y <- matrix(0, nrow = length(y), ncol = term)
  Y[, 1] <- 1
  if (term >= 2) Y[, 2] <- log(y / (1 - y))
  if (term >= 3) Y[, 3] <- (y - 0.5) * Y[, 2]
  if (term >= 4) Y[, 4] <- (y - 0.5)
  if (term >= 5) {
    o <- 2; e <- 2
    for (i in 5:term) {
      if (i %% 2 != 0) { o <- o + 1; Y[, i] <- (y - 0.5)^o } else { e <- e + 1; Y[, i] <- Y[, 2] * (y - 0.5)^e }
    }
  }
  Y
}

#' @keywords internal
metalog_to_z <- function(x, bounds, boundedness) {
  switch(boundedness,
    u = x,
    sl = log(x - bounds[1]),
    su = -log(bounds[2] - x),
    b = log((x - bounds[1]) / (bounds[2] - x))
  )
}

#' @keywords internal
metalog_from_z <- function(z, bounds, boundedness) {
  switch(boundedness,
    u = z,
    sl = bounds[1] + exp(z),
    su = bounds[2] - exp(-z),
    b = (bounds[1] + bounds[2] * exp(z)) / (1 + exp(z))
  )
}

#' Fit a metalog distribution via exact linear solve
#'
#' When the number of interior percentiles (`p` strictly in `(0,1)`) equals
#' the number of metalog terms, the fit is an EXACTLY-DETERMINED linear
#' system (`solve(Y, z)`), not an optimization - validated (upstream) to
#' ~1e-12 against `rmetalog::metalog()`. This is the fast path `rmetalog`
#' itself takes when its solution is already feasible; see
#' `check_metalog_feasible()`/`quantile_metalog_with_fallback()` for the
#' guard `rmetalog`'s own LP feasibility-correction would otherwise provide.
#'
#' @param interior_values,interior_probs Percentile values/probabilities,
#'   both strictly interior (excluding p=0/p=1 if present).
#' @param bounds `c(lower, upper)`; required unless `boundedness="u"`.
#' @param boundedness One of `"u"` (unbounded), `"sl"` (semi-bounded below),
#'   `"su"` (semi-bounded above), `"b"` (bounded).
#' @return `list(a=, term=, bounds=, boundedness=)`.
#' @export
fit_metalog_linear <- function(interior_values, interior_probs, bounds = NULL, boundedness = "u") {
  term <- length(interior_probs)
  Y <- metalog_basis_matrix(interior_probs, term)
  z <- metalog_to_z(interior_values, bounds, boundedness)
  a <- as.vector(solve(Y, z))
  list(a = a, term = term, bounds = bounds, boundedness = boundedness)
}

#' Evaluate a metalog quantile function
#'
#' Unlike the raster-native original (which evaluates one shared `q` across
#' many cells' coefficient rasters), here `fit` is a single set of
#' coefficients and `q` may be a vector - so this multiplies the basis
#' matrix evaluated at `q` by the coefficient vector directly, rather than
#' the raster version's "select row 1" pattern.
#'
#' @param fit Output of `fit_metalog_linear()`.
#' @param q Vector of probabilities in (0,1).
#' @export
quantile_metalog_linear <- function(fit, q) {
  Y_q <- metalog_basis_matrix(q, fit$term)
  z_q <- as.vector(Y_q %*% fit$a)
  metalog_from_z(z_q, fit$bounds, fit$boundedness)
}

#' Check whether a metalog fit's quantile function is monotonic (feasible)
#'
#' Probes `quantile_metalog_linear()` at a grid of probabilities and flags
#' non-monotonicity (equivalent to the implied density going negative
#' somewhere) - a probe, not a proof, but cheap and effective.
#'
#' @param fit Output of `fit_metalog_linear()`.
#' @param y_grid Probability grid to probe.
#' @return `TRUE` if infeasible (non-monotonic), else `FALSE`.
#' @export
check_metalog_feasible <- function(fit, y_grid = seq(0.02, 0.98, by = 0.02)) {
  probe <- quantile_metalog_linear(fit, y_grid)
  any(diff(probe) < 0)
}

#' Metalog quantile with automatic fallback to `linear_cdf`
#'
#' @param fit Output of `fit_metalog_linear()`.
#' @param infeasible Output of `check_metalog_feasible()` for this same `fit`.
#' @param full_probs,full_values The FULL percentile set (including p=0/p=1
#'   if available) used for the `linear_cdf` fallback.
#' @param q Vector of probabilities.
#' @export
quantile_metalog_with_fallback <- function(fit, infeasible, full_probs, full_values, q) {
  if (isTRUE(infeasible)) {
    quantile_linear_cdf(full_probs, full_values, q)
  } else {
    quantile_metalog_linear(fit, q)
  }
}

# ---------------------------------------------------------------------------
# "auto" family resolution and the single fit_percentile_triplet() dispatcher
# ---------------------------------------------------------------------------

#' Resolve a property's distribution family from its own percentile skew
#'
#' Adapted from `code_ref/reanalysis-platform/property_fusion_dispatch.R`'s
#' `resolve_property_dist()`: computed directly from this row's own l/r/h
#' (no spatial/AOI aggregation needed in a tabular context). Deliberately
#' narrow - never resolves to `"metalog"` (that stays an explicit config
#' choice), matching the reference's own documented restriction.
#'
#' @param l,r,h Low/representative/high values.
#' @param bounds Optional `c(lower, upper)`; when supplied, always resolves
#'   to `"beta"`.
#' @param skew_threshold Absolute skew-proxy threshold below which `"normal"`
#'   is chosen over `"lognormal"`.
#' @return `list(family=, skew_proxy=)`.
#' @export
resolve_property_family <- function(l, r, h, bounds = NULL, skew_threshold = 0.15) {
  skew_proxy <- ((h - r) - (r - l)) / (h - l)
  family <- if (!is.null(bounds)) {
    "beta"
  } else if (is.finite(skew_proxy) && abs(skew_proxy) < skew_threshold) {
    "normal"
  } else {
    "lognormal"
  }
  list(family = family, skew_proxy = skew_proxy)
}

#' Fit a distribution to a low/representative/high percentile triplet
#'
#' The single dispatcher `extract_property_parameters()` (`monte-carlo.R`)
#' calls to turn a SSURGO `_l/_r/_h` triplet into a family-appropriate fit.
#'
#' @param l,r,h Low/representative/high values.
#' @param family One of `"triangular"`, `"uniform"`, `"normal"`,
#'   `"lognormal"`, `"beta"`, `"metalog"`, `"linear_cdf"`, `"auto"`.
#' @param lh_probs Length-2 vector giving the probabilities `l`/`h` are
#'   assumed to represent (SSURGO's low/high are professional-judgment
#'   bounds, not confirmed exact percentiles - this makes that assumption
#'   explicit and overridable rather than silently hardcoded).
#' @param bounds Optional `c(lower, upper)` physical bounds; required for
#'   `"beta"` and for `"metalog"` when `boundedness != "u"`.
#' @param boundedness Metalog boundedness; defaults to `"u"` when `bounds`
#'   is `NULL`, else `"b"`.
#' @return `list(family=, fit=, valid=)`. `valid=FALSE` when any of
#'   `l`/`r`/`h` is non-finite.
#' @export
fit_percentile_triplet <- function(l, r, h, family,
                                    lh_probs = c(0.05, 0.95),
                                    bounds = NULL,
                                    boundedness = if (is.null(bounds)) "u" else "b") {
  if (!all(is.finite(c(l, r, h)))) {
    return(list(family = family, fit = NULL, valid = FALSE))
  }

  # A zero-width interval is degenerate for every parametric family (division
  # by zero in the SD/shape-parameter formulas) - the only honest
  # representation is "always draw exactly r", regardless of the requested
  # family.
  if (isTRUE(all.equal(h, l))) {
    return(list(family = "triangular", fit = list(min = r, mode = r, max = r), valid = TRUE))
  }

  if (identical(family, "auto")) {
    resolved <- resolve_property_family(l, r, h, bounds = bounds)
    family <- resolved$family
    # log() is undefined for non-positive values - fall back to normal
    # exactly as mod04's own get_appropriate_distributions() already skips
    # lognormal/gamma for non-positive data.
    if (identical(family, "lognormal") && any(c(l, r, h) <= 0)) family <- "normal"
  }

  probs <- c(lh_probs[1], 0.5, lh_probs[2])
  values <- c(l, r, h)

  fit <- switch(family,
    triangular = list(min = l, mode = r, max = h),
    uniform = list(min = l, max = h),
    normal = fit_normal_triplet(l, r, h, lh_probs[1], lh_probs[2]),
    lognormal = {
      if (any(values <= 0)) {
        return(list(family = "normal", fit = fit_normal_triplet(l, r, h, lh_probs[1], lh_probs[2]), valid = TRUE))
      }
      fit_normal_triplet(log(l), log(r), log(h), lh_probs[1], lh_probs[2])
    },
    beta = {
      if (is.null(bounds)) stop("fit_percentile_triplet(): bounds = c(lower, upper) is required for family = 'beta'.")
      c(fit_beta_mle_newton(values, bounds), list(lower = bounds[1], upper = bounds[2]))
    },
    metalog = {
      metalog_fit <- fit_metalog_linear(values, probs, bounds = bounds %||% c(0, 1), boundedness = boundedness)
      c(metalog_fit, list(
        infeasible = check_metalog_feasible(metalog_fit),
        fallback_probs = probs,
        fallback_values = values
      ))
    },
    linear_cdf = list(probs = probs, values = values),
    stop(sprintf("fit_percentile_triplet(): unknown family '%s'.", family))
  )

  list(family = family, fit = fit, valid = TRUE)
}

#' Evaluate the quantile function of a fitted percentile-triplet distribution
#'
#' Replaces the old `transform_to_distribution()` design: a thin dispatcher
#' from `family` to the matching `quantile_*()` function.
#'
#' @param u Vector of probabilities (e.g. correlated uniform draws).
#' @param family One of `fit_percentile_triplet()`'s resolved families.
#' @param fit The `fit` element of `fit_percentile_triplet()`'s return value.
#' @return Numeric vector, same length as `u`. `NA_real_` (not a silent wrong
#'   value or error) when `fit` is `NULL`/invalid.
#' @export
quantile_from_fit <- function(u, family, fit) {
  if (is.null(fit)) return(rep(NA_real_, length(u)))

  switch(family,
    triangular = quantile_triangular(fit, u),
    uniform = qunif(u, fit$min, fit$max),
    normal = quantile_normal(fit, u),
    lognormal = exp(quantile_normal(fit, u)),
    beta = fit$lower + quantile_beta(fit, u) * (fit$upper - fit$lower),
    metalog = quantile_metalog_with_fallback(fit, fit$infeasible, fit$fallback_probs, fit$fallback_values, u),
    linear_cdf = quantile_linear_cdf(fit$probs, fit$values, u),
    rep(NA_real_, length(u))
  )
}

#' Validate a fitted distribution's parameters
#'
#' Real replacement for the always-`valid=TRUE` `validate_distribution_parameters()`
#' stub previously in `mod05_monte_carlo.R`.
#'
#' @param family One of `fit_percentile_triplet()`'s resolved families.
#' @param fit The `fit` element of `fit_percentile_triplet()`'s return value.
#' @return `list(valid=, message=)`.
#' @export
validate_fit_parameters <- function(family, fit) {
  if (is.null(fit)) return(list(valid = FALSE, message = "fit is NULL"))

  switch(family,
    triangular = ,
    uniform = {
      mn <- fit$min; mx <- fit$max
      if (!is.finite(mn) || !is.finite(mx) || mn > mx) {
        return(list(valid = FALSE, message = "min must be finite and <= max"))
      }
      if (!is.null(fit$mode) && (fit$mode < mn || fit$mode > mx)) {
        return(list(valid = FALSE, message = "mode outside [min, max]"))
      }
      list(valid = TRUE, message = "")
    },
    normal = ,
    lognormal = {
      if (!is.finite(fit$mean) || !is.finite(fit$sd) || fit$sd < 0) {
        return(list(valid = FALSE, message = "mean/sd must be finite and sd >= 0"))
      }
      list(valid = TRUE, message = "")
    },
    beta = {
      if (!is.finite(fit$shape1) || !is.finite(fit$shape2) || fit$shape1 <= 0 || fit$shape2 <= 0) {
        return(list(valid = FALSE, message = "shape1/shape2 must be positive and finite"))
      }
      if (!is.finite(fit$lower) || !is.finite(fit$upper) || fit$lower >= fit$upper) {
        return(list(valid = FALSE, message = "lower must be < upper"))
      }
      list(valid = TRUE, message = "")
    },
    metalog = {
      # Deliberately NOT gated on `!fit$infeasible`: an infeasible metalog fit
      # still produces correct output via the automatic linear_cdf fallback
      # in quantile_metalog_with_fallback()/quantile_from_fit() - flagging it
      # invalid here would cause callers (e.g. mod05's setup_distributions())
      # to discard a perfectly usable (if degraded-to-nonparametric) fit in
      # favor of a cruder triangular guess. Only check structural completeness.
      required <- c("a", "term", "bounds", "boundedness", "infeasible", "fallback_probs", "fallback_values")
      if (!all(required %in% names(fit))) {
        return(list(valid = FALSE, message = "metalog fit is missing required fields"))
      }
      list(valid = TRUE, message = if (isTRUE(fit$infeasible)) "infeasible fit; using linear_cdf fallback" else "")
    },
    linear_cdf = {
      if (length(fit$probs) < 2 || length(fit$values) != length(fit$probs)) {
        return(list(valid = FALSE, message = "need >= 2 matching probs/values"))
      }
      list(valid = TRUE, message = "")
    },
    list(valid = FALSE, message = paste("unknown family:", family))
  )
}

# ==============================================================================
# 2. ILR (ISOMETRIC LOG-RATIO) TRANSFORM + MONTE CARLO MOMENTS
# ==============================================================================
#
# Dependency-free port of code_ref/reanalysis-platform/texture_ilr_fusion.R,
# validated there to match compositions::ilr() exactly and round-trip to
# ~1e-14. NOTE: the balance hierarchy (position 1 vs {2,3}, then 2 vs 3) is
# FIXED by these formulas. The `clay`/`sand`/`silt` parameter and output
# names below are POSITIONAL-ROLE placeholders (position 1/2/3 in the
# sequential binary partition), not an identity requirement - callers
# determine which real property occupies which position, not this file.
# `monte-carlo.R`'s `composition_groups$texture$members` is the single
# source of truth for that mapping (default: sandtotal, silttotal,
# claytotal - chosen to match the sequential binary partition
# `compositions::ilr()` used for the optional KSSL reference correlation
# matrices in `kssl-reference-correlations.R`: sand vs {silt,clay}, then
# silt vs clay).

#' Forward ILR transform: 2 unconstrained coordinates from a 3-part composition
#'
#' `clay`/`sand`/`silt` name positions 1/2/3 of the sequential binary
#' partition (position 1 vs positions 2 and 3, then 2 vs 3) - which real
#' property occupies which position is entirely up to the caller (see the
#' file header comment above); the names are historical, not an identity
#' requirement.
#' @param clay,sand,silt Numeric vectors, same units (e.g. percent), same length.
#' @return A 2-column matrix (`z1`, `z2`).
#' @export
ilr_forward <- function(clay, sand, silt) {
  cbind(z1 = sqrt(2 / 3) * log(clay / sqrt(sand * silt)), z2 = sqrt(1 / 2) * log(sand / silt))
}

#' Inverse ILR transform: 2 unconstrained coordinates -> a valid composition
#'
#' Output columns `clay`/`sand`/`silt` are positions 1/2/3 of the same
#' sequential binary partition `ilr_forward()` uses - see that function's doc
#' and the file header comment above for the positional-role vs identity
#' distinction.
#' @param z1,z2 Numeric vectors of ILR coordinates.
#' @param total The composition's target sum (default 100).
#' @return A 3-column matrix (`clay`, `sand`, `silt`) summing to `total`.
#' @export
ilr_inverse <- function(z1, z2, total = 100) {
  clr_clay <- z1 * sqrt(2 / 3)
  clr_sand <- z1 * (-sqrt(1 / 6)) + z2 * sqrt(1 / 2)
  clr_silt <- z1 * (-sqrt(1 / 6)) + z2 * (-sqrt(1 / 2))
  e_clay <- exp(clr_clay); e_sand <- exp(clr_sand); e_silt <- exp(clr_silt)
  s <- e_clay + e_sand + e_silt
  cbind(clay = total * e_clay / s, sand = total * e_sand / s, silt = total * e_silt / s)
}

#' Estimate one horizon's ILR-space mean/covariance from marginal l/r/h triplets
#'
#' Draws each fraction's own marginal Normal (clamped `> 0.01`), CLOSES each
#' draw (renormalizes the 3 draws to sum to 100 - the standard
#' compositional-data "closure" operation, which is what induces the correct
#' negative cross-fraction correlation, not an assumption), ILR-transforms,
#' and takes the empirical mean/covariance. A cheaper closed-form delta-method
#' shortcut was tried upstream and REJECTED: it underestimates the true
#' covariance by 5-10x once per-fraction CV exceeds ~30% (common near
#' boundary fractions) - do not substitute it for this Monte Carlo version.
#'
#' @param low_clay,rep_clay,high_clay,low_sand,rep_sand,high_sand,low_silt,rep_silt,high_silt
#'   Low/representative/high values for each fraction.
#' @param z Standard-normal quantile matching the low/high interval (e.g.
#'   `qnorm(0.95)` for a 5th/95th-percentile low/high).
#' @param n_mc Monte Carlo sample size.
#' @return `list(mu = c(z1, z2), Sigma = 2x2 matrix)`.
#' @export
estimate_ilr_moments_mc <- function(low_clay, rep_clay, high_clay,
                                     low_sand, rep_sand, high_sand,
                                     low_silt, rep_silt, high_silt,
                                     z, n_mc = 2000) {
  sd_clay <- (high_clay - low_clay) / (2 * z)
  sd_sand <- (high_sand - low_sand) / (2 * z)
  sd_silt <- (high_silt - low_silt) / (2 * z)

  clay_s <- pmax(rnorm(n_mc, rep_clay, sd_clay), 0.01)
  sand_s <- pmax(rnorm(n_mc, rep_sand, sd_sand), 0.01)
  silt_s <- pmax(rnorm(n_mc, rep_silt, sd_silt), 0.01)
  total <- clay_s + sand_s + silt_s
  z_samples <- ilr_forward(clay_s / total * 100, sand_s / total * 100, silt_s / total * 100)

  list(mu = colMeans(z_samples), Sigma = stats::cov(z_samples))
}

#' Draw posterior composition samples from a fused/fitted ILR-space (mu, Sigma)
#'
#' Guaranteed (by construction, via `ilr_inverse()`) to sum to `total` and
#' stay within `[0, total]` for every draw.
#'
#' @param mu Length-2 mean vector.
#' @param Sigma 2x2 covariance matrix.
#' @param n Number of samples to draw.
#' @param total Composition target sum.
#' @return An `n x 3` matrix (`clay`, `sand`, `silt`).
#' @export
sample_ilr_posterior <- function(mu, Sigma, n = 1000, total = 100) {
  L <- chol(Sigma)
  z_draws <- matrix(rnorm(2 * n), ncol = 2) %*% L
  z_draws[, 1] <- z_draws[, 1] + mu[1]
  z_draws[, 2] <- z_draws[, 2] + mu[2]
  ilr_inverse(z_draws[, 1], z_draws[, 2], total = total)
}

# ==============================================================================
# 3. MATRIX UTILITIES (shared correlation-matrix repair/validation/estimation)
# ==============================================================================

#' Repair a near-correlation matrix to be positive definite
#'
#' Floors small/negative eigenvalues, reconstructs, and rescales back to a
#' unit-diagonal correlation matrix. Relocated verbatim from
#' `mod05_monte_carlo.R` (already correct there - not one of the confirmed
#' bugs, just moved so `statistics.R` can share it too).
#'
#' @param matrix A square numeric matrix.
#' @param min_eigenvalue Eigenvalue floor.
#' @return A positive-definite correlation matrix, same dimensions.
#' @export
ensure_positive_definite_matrix <- function(matrix, min_eigenvalue = 1e-6) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
    return(matrix)
  }

  # eigen()$vectors carries no dimnames of its own, so the reconstructed
  # matrix below would otherwise silently lose whatever row/column names the
  # input had (a real, previously-unnoticed gap: callers indexing the result
  # by property name, e.g. result["propA","propB"], would fail even though
  # the underlying values were correct) - preserve them explicitly.
  original_dimnames <- dimnames(matrix)

  eigen_decomp <- eigen(matrix)
  eigenvals <- eigen_decomp$values
  eigenvecs <- eigen_decomp$vectors

  eigenvals[eigenvals < min_eigenvalue] <- min_eigenvalue

  adjusted_matrix <- eigenvecs %*% diag(eigenvals) %*% t(eigenvecs)

  diag_sqrt <- sqrt(diag(adjusted_matrix))
  result <- adjusted_matrix / outer(diag_sqrt, diag_sqrt)
  dimnames(result) <- original_dimnames
  result
}

#' Validate a correlation matrix's shape and positive-definiteness
#'
#' Relocated verbatim from `mod05_monte_carlo.R`.
#'
#' @param corr_matrix A candidate correlation matrix.
#' @param properties Character vector the matrix's dimensions should match.
#' @return `list(valid=, message=)`.
#' @export
validate_correlation_matrix <- function(corr_matrix, properties) {
  validation <- list(valid = TRUE, message = "")

  if (!is.matrix(corr_matrix)) {
    return(list(valid = FALSE, message = "Not a matrix"))
  }
  if (nrow(corr_matrix) != ncol(corr_matrix)) {
    return(list(valid = FALSE, message = "Not square matrix"))
  }
  if (nrow(corr_matrix) != length(properties)) {
    return(list(valid = FALSE, message = "Matrix dimensions don't match number of properties"))
  }
  if (!all(abs(diag(corr_matrix) - 1) < 1e-10)) {
    return(list(valid = FALSE, message = "Diagonal elements are not 1"))
  }

  eigenvals <- eigen(corr_matrix, only.values = TRUE)$values
  if (min(eigenvals) <= 0) {
    return(list(valid = FALSE, message = "Matrix is not positive definite"))
  }

  validation
}

#' Estimate a correlation matrix robustly from data, with grouped and global fallbacks
#'
#' Ports `code_ref/brdf/property_simulation.R`'s `simulate_soil_properties()`
#' correlation pattern: `Hmisc::rcorr()` per group when the group has enough
#' complete observations, else falling back; forces symmetry; repairs
#' non-positive-definite matrices via `Matrix::nearPD()`. When `group_var` is
#' supplied, each qualifying group's empirical matrix is combined into a
#' single size-weighted average (mod05's simulation architecture uses one
#' correlation matrix for the whole run, not a per-stratum matrix, so this
#' folds stratification information in without requiring a larger
#' architecture change).
#'
#' @param data A data frame of numeric property columns (one row per
#'   horizon), optionally with a grouping column named by `group_var`.
#' @param group_var Optional column name to stratify by (e.g. `"genhz"`).
#' @param min_group_n Minimum complete-observation count required to trust a
#'   group's (or the pooled data's) empirical correlation over the fallback.
#' @param global_fallback Optional fallback correlation matrix; defaults to
#'   the identity matrix over the numeric columns found in `data`.
#' @param group_fallback_matrices Optional named list of correlation
#'   matrices, keyed by `group_var` value, used in place of dropping a group
#'   that fails to meet `min_group_n` (today's default behavior when this is
#'   `NULL`). A substituted group contributes to the same sample-size-weighted
#'   combination as successfully-estimated groups, with nominal weight
#'   `n_obs = min_group_n` (these matrices carry no real sample-size
#'   metadata, so `min_group_n` - the same bar local data must already clear
#'   to stand alone - is used as a legible, already-plumbed-through stand-in;
#'   any single successful empirical group already outweighs one substituted
#'   group, and multiple dominate further). A group with no empirical result
#'   AND no matching entry here is still dropped, exactly as before.
#' @param global_fallback_method Label for the final-tier `method` when
#'   `global_fallback` is used because no group (empirical or substituted)
#'   and no pooled fit succeeded. Defaults to `"global_fallback"` (today's
#'   exact string); callers supplying a non-identity `global_fallback` (e.g.
#'   a KSSL-derived prior) should override this so the returned `method`
#'   honestly distinguishes "genuinely uninformative" from "a real borrowed
#'   prior."
#' @return `list(matrix=, method=, n_obs=)`. `method` is one of
#'   `"empirical_grouped"`, `"kssl_fallback_grouped"`,
#'   `"empirical_grouped_kssl_blended"`, `"empirical_pooled"`, or
#'   `global_fallback_method`'s value (`"global_fallback"` by default). The
#'   two `kssl_*`-named values are generic outcomes of any
#'   `group_fallback_matrices` use, not KSSL-specific logic in this function
#'   - named for the fallback source this function was built to support
#'   (`kssl-reference-correlations.R`).
#' @export
estimate_correlation_matrix_robust <- function(data, group_var = NULL, min_group_n = 5,
                                                global_fallback = NULL,
                                                group_fallback_matrices = NULL,
                                                global_fallback_method = "global_fallback") {
  numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, group_var)
  n_vars <- length(numeric_cols)

  if (is.null(global_fallback)) {
    global_fallback <- diag(n_vars)
  }
  rownames(global_fallback) <- colnames(global_fallback) <- numeric_cols

  repair_and_symmetrize <- function(m) {
    m <- (m + t(m)) / 2
    is_pd <- tryCatch({ chol(m); TRUE }, error = function(e) FALSE)
    if (!is_pd) {
      m <- tryCatch(as.matrix(Matrix::nearPD(m, corr = TRUE)$mat), error = function(e) NULL)
    }
    m
  }

  compute_one <- function(sub_data) {
    if (n_vars < 2) return(NULL)
    sub_numeric <- as.matrix(sub_data[, numeric_cols, drop = FALSE])
    complete_rows <- rowSums(!is.na(sub_numeric)) == n_vars
    sub_numeric <- sub_numeric[complete_rows, , drop = FALSE]
    if (nrow(sub_numeric) <= min_group_n) return(NULL)

    local_cor <- tryCatch(Hmisc::rcorr(sub_numeric, type = "pearson")$r, error = function(e) NULL)
    if (is.null(local_cor) || anyNA(local_cor) || any(!is.finite(local_cor))) return(NULL)

    local_cor <- repair_and_symmetrize(local_cor)
    if (is.null(local_cor)) return(NULL)
    list(matrix = local_cor, n_obs = nrow(sub_numeric))
  }

  if (!is.null(group_var) && group_var %in% names(data)) {
    groups <- unique(data[[group_var]][!is.na(data[[group_var]])])
    per_group <- lapply(groups, function(g) {
      result <- compute_one(data[!is.na(data[[group_var]]) & data[[group_var]] == g, , drop = FALSE])
      if (!is.null(result)) {
        result$source <- "empirical"
        return(result)
      }
      fb <- group_fallback_matrices[[as.character(g)]]
      if (!is.null(fb) && all(numeric_cols %in% rownames(fb))) {
        fb <- repair_and_symmetrize(fb[numeric_cols, numeric_cols, drop = FALSE])
        if (!is.null(fb)) {
          return(list(matrix = fb, n_obs = min_group_n, source = "kssl_fallback"))
        }
      }
      NULL
    })
    names(per_group) <- as.character(groups)
    valid <- Filter(Negate(is.null), per_group)

    if (length(valid) > 0) {
      weights <- vapply(valid, function(x) x$n_obs, numeric(1))
      combined <- Reduce(`+`, Map(function(x, w) x$matrix * w, valid, weights)) / sum(weights)
      combined <- repair_and_symmetrize(combined)
      if (!is.null(combined)) {
        sources <- vapply(valid, function(x) x$source, character(1))
        method <- if (all(sources == "empirical")) {
          "empirical_grouped"
        } else if (all(sources == "kssl_fallback")) {
          "kssl_fallback_grouped"
        } else {
          "empirical_grouped_kssl_blended"
        }
        return(list(matrix = combined, method = method, n_obs = weights))
      }
    }
  }

  pooled <- compute_one(data)
  if (!is.null(pooled)) {
    return(list(matrix = pooled$matrix, method = "empirical_pooled", n_obs = pooled$n_obs))
  }

  list(matrix = global_fallback, method = global_fallback_method, n_obs = 0)
}

# ==============================================================================
# 4. COMPOSITION-GROUP ORCHESTRATION (used by monte-carlo.R)
# ==============================================================================

#' Expand a requested properties vector around active composition groups
#'
#' If ALL members of a configured composition group (e.g. `"texture"`:
#' clay/sand/silt) are present in `properties`, they are replaced (position
#' of the first member preserved) with the group's pseudo-property names
#' (e.g. `"ilr1"`/`"ilr2"`) so they can flow through the existing generic
#' Cholesky-copula simulation machinery like any other property. If only
#' some members are present, the group stays inactive (a logged WARN) and
#' those properties simulate independently via the legacy path.
#'
#' `members` declares which real property occupies `ilr_forward()`/
#' `ilr_inverse()`'s position 1/2/3 (their `clay`/`sand`/`silt` naming is a
#' positional-role placeholder, not an identity requirement - see the header
#' comment in the ILR section of this file). `restore_composition_properties()`
#' maps `ilr_inverse()`'s fixed `clay`/`sand`/`silt` output columns onto
#' `members[1]`/`members[2]`/`members[3]` by POSITION, not by name matching.
#'
#' @param properties Character vector of requested property names.
#' @param config A Monte Carlo config list; reads
#'   `config$monte_carlo$composition_groups` (a named list, each entry
#'   `list(members=, pseudo=)`).
#' @return `list(sim_properties=, groups=)` where `groups` is a named list of
#'   `list(members=, pseudo=, active=)`.
#' @export
resolve_composition_groups <- function(properties, config) {
  group_defs <- config$monte_carlo$composition_groups %||% list()
  sim_properties <- properties
  groups <- list()

  for (group_name in names(group_defs)) {
    members <- group_defs[[group_name]]$members
    pseudo <- group_defs[[group_name]]$pseudo
    present <- intersect(members, sim_properties)

    if (length(present) == length(members)) {
      inserted <- FALSE
      new_props <- character(0)
      for (p in sim_properties) {
        if (p %in% members) {
          if (!inserted) {
            new_props <- c(new_props, pseudo)
            inserted <- TRUE
          }
        } else {
          new_props <- c(new_props, p)
        }
      }
      sim_properties <- new_props
      groups[[group_name]] <- list(members = members, pseudo = pseudo, active = TRUE)
    } else {
      if (length(present) > 0) {
        log_message(
          "WARN",
          paste0(group_name, " composition_group requires all of ", paste(members, collapse = "/"),
                 " (in that exact configured role order); found ", length(present), "/", length(members),
                 " - falling back to independent per-property simulation for this group"),
          category = "MonteCarlo"
        )
      }
      groups[[group_name]] <- list(members = members, pseudo = pseudo, active = FALSE)
    }
  }

  list(sim_properties = sim_properties, groups = groups)
}

#' Collapse simulated ILR pseudo-properties back to raw composition members
#'
#' For each active composition group, applies `ilr_inverse()` per horizon
#' (vectorized over realizations) and reassembles a
#' `[horizon, property, realization]` array dimnamed over the ORIGINAL
#' `properties` vector - non-group properties are copied through unchanged,
#' group members are filled from the inverse-ILR output (guaranteed
#' sum-to-`total` and in-bounds for every realization by construction), and
#' the pseudo-properties are dropped.
#'
#' @param simulation_results A `[horizon, property, realization]` array
#'   dimnamed over `sim_properties`.
#' @param sim_properties The (possibly ILR-substituted) property vector
#'   `simulation_results` is dimnamed over.
#' @param properties The caller-facing, original property vector.
#' @param groups The `groups` element of `resolve_composition_groups()`'s return value.
#' @return An array dimnamed over `properties`.
#' @export
restore_composition_properties <- function(simulation_results, sim_properties, properties, groups) {
  active_groups <- Filter(function(g) isTRUE(g$active), groups)
  if (length(active_groups) == 0) {
    return(simulation_results)
  }

  dims <- dim(simulation_results)
  n_horizons <- dims[1]
  n_realizations <- dims[3]

  out <- array(
    NA_real_,
    dim = c(n_horizons, length(properties), n_realizations),
    dimnames = list(
      horizon = dimnames(simulation_results)$horizon %||% seq_len(n_horizons),
      property = properties,
      realization = dimnames(simulation_results)$realization %||% seq_len(n_realizations)
    )
  )

  pseudo_all <- unlist(lapply(active_groups, `[[`, "pseudo"))
  passthrough <- setdiff(sim_properties, pseudo_all)
  for (p in passthrough) {
    out[, p, ] <- simulation_results[, match(p, sim_properties), ]
  }

  for (group in active_groups) {
    ilr1_idx <- match(group$pseudo[1], sim_properties)
    ilr2_idx <- match(group$pseudo[2], sim_properties)

    for (h in seq_len(n_horizons)) {
      composition <- ilr_inverse(simulation_results[h, ilr1_idx, ], simulation_results[h, ilr2_idx, ], total = 100)
      out[h, group$members[1], ] <- composition[, "clay"]
      out[h, group$members[2], ] <- composition[, "sand"]
      out[h, group$members[3], ] <- composition[, "silt"]
    }
  }

  out
}
