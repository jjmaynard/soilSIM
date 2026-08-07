#' @title Raster-Native Bayesian Fusion (Prior x Likelihood, per Cell)
#'
#' @description Per-cell Bayesian fusion of a prior belief distribution
#'   (e.g. an SSURGO-derived percentile raster set) with a likelihood (e.g.
#'   an independent percentile raster set from another source), across a
#'   whole `terra::SpatRaster` AOI at once. This is the raster-native
#'   counterpart of `R/bayesian-updating.R` (itself "intentionally NOT
#'   called anywhere in `monte-carlo.R`" - a standalone toolkit), ported
#'   from `code_ref/reanalysis-platform/{bayes_fuse.R,
#'   property_fusion_dispatch.R, reanalysis-platform-fusion.R}`.
#'
#'   Reuses `bayesian-updating.R`'s existing scalar/vector fusion functions
#'   directly rather than duplicating them - `bayes_update_normal_normal()`,
#'   `fuse_beta()`, `fuse_gamma()`, `moments_to_gamma()`/`moments_to_beta()`/
#'   `beta_to_moments()`/`gamma_to_moments()`, `normal_to_lognormal_params()`/
#'   `lognormal_to_normal_params()`, `bayesian_update()`, and
#'   `R/distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`, plus
#'   `fuse_bivariate_normal()` (also in `bayesian-updating.R`) - are all pure
#'   elementwise arithmetic, so they already work unchanged on `SpatRaster`
#'   inputs (confirmed by a dedicated smoke test before this file was
#'   written; see `bayes_update_normal_normal()`'s own doc comment upstream,
#'   which claims exactly this property).
#'
#' @section Deliberately out of scope:
#' `run_stage1_fusion()`/`run_stage1_fusion_group()` (the original top-level
#' orchestrators) were initially **not** ported, since at that point their
#' dependencies - `build_cache_key()`, `cache_get()`/`cache_set()`,
#' `fetch_ssurgo_percentiles()`, `fetch_solus_percentiles()` (itself calling
#' the live `soilDB::fetchSOLUS()` network API) - didn't exist anywhere in
#' this repo. The source bundle's own `HANDOFF_NOTES.md` confirms it was
#' built as a handoff package for a *different, sibling project*
#' ("reanalysis-platform"), not soilSIM, with these exact fetch/cache
#' functions listed as unresolved integration gaps for that project. Those
#' gaps have since been closed (`R/raster-cache.R`, `R/ssurgo-simulation.R`,
#' `R/solus-simulation.R`), and `run_stage1_fusion()`/
#' `run_stage1_fusion_group()` are now fully implemented below, wiring
#' everything together. `fuse_property_adaptive()` remains the lower-level
#' entry point for callers who already have their own pre-fetched
#' `prior_value_rasters`/`lik_value_rasters` and want to skip the
#' fetch-and-cache wrapper.
#'
#' Also not ported: the global `PROPERTIES` config-list registry from the
#' source bundle's `config.R`. `group_members()`/`fuse_property_adaptive()`
#' below take a per-call config list/member vector directly instead,
#' reusing soilSIM's own existing `config$monte_carlo$composition_groups`
#' convention (see `R/distributions.R`'s `resolve_composition_groups()`)
#' rather than introducing a second, parallel config schema.
#' @name raster_fusion
NULL

# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

#' Method-of-moments Gamma fit from a list of percentile-value rasters
#'
#' Thin raster wrapper around the existing scalar `moments_to_gamma()`
#' (`bayesian-updating.R`) - the mean/variance computation is the only
#' genuinely raster-specific part (combining a *list* of rasters via `Reduce()`).
#' @param value_rasters List of percentile-value SpatRasters.
#' @return `list(shape = SpatRaster, rate = SpatRaster)`.
#' @export
fit_gamma_mom_raster <- function(value_rasters) {
  k <- length(value_rasters)
  mean_r <- Reduce(`+`, value_rasters) / k
  var_r <- Reduce(`+`, lapply(value_rasters, function(r) (r - mean_r)^2)) / k
  moments_to_gamma(mean_r, var_r)
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

#' Closed-form same-family fusion path for `fuse_adaptive()`, with a
#' per-cell fallback to Normal-moment fusion (re-expressed back into the
#' requested family) where the naive same-family fusion is infeasible.
#' @param prior_value_rasters,lik_value_rasters Lists of percentile-value SpatRasters.
#' @param percentile_probs Numeric probabilities matching the value-raster list order.
#' @param family One of "normal", "beta", "gamma".
#' @param bounds Required for `family = "beta"`.
#' @return `list(posterior=, route_detail=, n_fallback_cells=)`.
#' @keywords internal
fuse_closed_form <- function(prior_value_rasters, lik_value_rasters, percentile_probs, family, bounds) {
  if (family == "normal") {
    idx <- percentile_index(percentile_probs)
    fit_prior <- fit_normal_raster(prior_value_rasters[[idx$lo_idx]], prior_value_rasters[[idx$p50_idx]], prior_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
    fit_lik <- fit_normal_raster(lik_value_rasters[[idx$lo_idx]], lik_value_rasters[[idx$p50_idx]], lik_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
    posterior <- bayes_update_normal_normal(fit_prior$mu, fit_prior$sigma, fit_lik$mu, fit_lik$sigma)
    return(list(posterior = posterior, route_detail = NULL, n_fallback_cells = 0))
  }

  if (family == "beta") {
    fit_prior <- fit_beta_mle_newton_raster(prior_value_rasters, bounds)
    fit_lik <- fit_beta_mle_newton_raster(lik_value_rasters, bounds)
    fused <- fuse_beta(fit_prior$alpha, fit_prior$beta, fit_lik$alpha, fit_lik$beta)

    prior_m <- beta_to_moments(fit_prior$alpha, fit_prior$beta)
    lik_m <- beta_to_moments(fit_lik$alpha, fit_lik$beta)
    span <- bounds[2] - bounds[1]
    fallback_normal <- bayes_update_normal_normal(
      bounds[1] + prior_m$mean * span, sqrt(prior_m$var) * span,
      bounds[1] + lik_m$mean * span, sqrt(lik_m$var) * span
    )
    fallback_mean <- (fallback_normal$mu - bounds[1]) / span
    fallback_var <- (fallback_normal$sigma / span)^2
    fallback_beta <- moments_to_beta(fallback_mean, fallback_var)

    alpha_final <- terra::ifel(fused$feasible, fused$alpha, fallback_beta$alpha)
    beta_final <- terra::ifel(fused$feasible, fused$beta, fallback_beta$beta)
    route_detail <- terra::ifel(fused$feasible, 1, 0)
    n_fallback <- sum(!terra::values(fused$feasible), na.rm = TRUE)

    return(list(posterior = list(alpha = alpha_final, beta = beta_final), route_detail = route_detail, n_fallback_cells = n_fallback))
  }

  if (family == "gamma") {
    fit_prior <- fit_gamma_mom_raster(prior_value_rasters)
    fit_lik <- fit_gamma_mom_raster(lik_value_rasters)
    fused <- fuse_gamma(fit_prior$shape, fit_prior$rate, fit_lik$shape, fit_lik$rate)

    prior_m <- gamma_to_moments(fit_prior$shape, fit_prior$rate)
    lik_m <- gamma_to_moments(fit_lik$shape, fit_lik$rate)
    fallback_normal <- bayes_update_normal_normal(prior_m$mean, sqrt(prior_m$var), lik_m$mean, sqrt(lik_m$var))
    fallback_gamma <- moments_to_gamma(fallback_normal$mu, fallback_normal$sigma^2)

    shape_final <- terra::ifel(fused$feasible, fused$shape, fallback_gamma$shape)
    rate_final <- terra::ifel(fused$feasible, fused$rate, fallback_gamma$rate)
    route_detail <- terra::ifel(fused$feasible, 1, 0)
    n_fallback <- sum(!terra::values(fused$feasible), na.rm = TRUE)

    return(list(posterior = list(shape = shape_final, rate = rate_final), route_detail = route_detail, n_fallback_cells = n_fallback))
  }
}

#' Fully general fusion path for `fuse_adaptive()`: per cell, draw samples
#' from both sides' percentiles (via `simulate_from_percentiles(method=
#' "linear_cdf")`), fuse via `bayesian_update()`, and moment-match the
#' posterior samples back into the requested family so the output contract
#' matches the closed-form route.
#' @inheritParams fuse_closed_form
#' @param n_samples Samples drawn per side, per cell.
#' @param grid_resolution Passed to `bayesian_update()`.
#' @keywords internal
fuse_general_kde <- function(prior_value_rasters, lik_value_rasters, percentile_probs,
                              family, bounds, n_samples, grid_resolution) {
  percentile_cols <- paste0("P", round(percentile_probs * 100))
  k <- length(percentile_probs)
  prior_stack <- terra::rast(prior_value_rasters); names(prior_stack) <- paste0("prior_", percentile_cols)
  lik_stack <- terra::rast(lik_value_rasters); names(lik_stack) <- paste0("lik_", percentile_cols)
  combined <- c(prior_stack, lik_stack)

  # terra::app() calls fun with a plain matrix (cells x layers) for a
  # chunk, and makes an initial probe call with a bare vector (not a
  # matrix) to infer the output layer count - fun must normalize both shapes.
  #
  # Real-world discovery (not in the original ported source): cells with too many NA percentile
  # values (e.g. at raster edges after terra::resample()/align_percentile_probs() alignment
  # between two differently-extented grids) make extract_percentile_pairs() hard-stop() with
  # "Need at least 2 valid quantile columns." terra::app() is not per-cell-error-tolerant by
  # default, so one such cell would abort the entire raster operation - wrapped per-cell so a bad
  # cell degrades to NA output instead of crashing the whole fusion.
  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1)
    t(apply(row_mat, 1, function(row) {
      tryCatch({
        prior_row <- as.data.frame(stats::setNames(as.list(row[1:k]), percentile_cols))
        lik_row <- as.data.frame(stats::setNames(as.list(row[(k + 1):(2 * k)]), percentile_cols))
        prior_samples <- simulate_from_percentiles(prior_row, method = "linear_cdf", percentile_cols = percentile_cols, n = n_samples)
        lik_samples <- simulate_from_percentiles(lik_row, method = "linear_cdf", percentile_cols = percentile_cols, n = n_samples)
        post <- if (is.null(grid_resolution)) bayesian_update(prior_samples, lik_samples) else bayesian_update(prior_samples, lik_samples, grid_resolution = grid_resolution)
        c(mean = mean(post), var = stats::var(post))
      }, error = function(e) c(mean = NA_real_, var = NA_real_))
    }))
  }
  moments_r <- terra::app(combined, fun)
  mean_r <- moments_r[[1]]; var_r <- moments_r[[2]]

  posterior <- switch(family,
    normal = list(mu = mean_r, sigma = sqrt(var_r)),
    beta = { span <- bounds[2] - bounds[1]; moments_to_beta((mean_r - bounds[1]) / span, var_r / span^2) },
    gamma = moments_to_gamma(mean_r, var_r)
  )
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
#'   `bayesian_update()` route is used; above it, the closed-form route.
#' @param n_samples Samples drawn per side for the general route's KDE fusion.
#' @param grid_resolution Passed to `bayesian_update()` for the general route.
#' @param verbose If TRUE (default), print the chosen route and why.
#' @return A list:
#'   - `posterior`: family-native parameter list (`mu`/`sigma` for "normal",
#'     `alpha`/`beta` for "beta", `shape`/`rate` for "gamma") - always this
#'     shape regardless of which route ran.
#'   - `route`: one of `"bayesian_update_general"`, `"closed_form_normal"`,
#'     `"closed_form_beta"`, `"closed_form_gamma"`.
#'   - `route_detail`: for `family %in% c("beta","gamma")` on the
#'     closed-form route, a `SpatRaster` of 1 (native same-family fusion) /
#'     0 (per-cell fallback). `NULL` for "normal" and for the general route.
#'   - `n_fallback_cells`: count of cells that used the fallback.
#'   - `diagnostics`: `list(ncell=, threshold_cells=, family=, elapsed_sec=)`.
#' @export
fuse_adaptive <- function(prior_value_rasters, lik_value_rasters, percentile_probs,
                           family = c("normal", "beta", "gamma"), bounds = NULL,
                           threshold_cells = 80000, n_samples = 500,
                           grid_resolution = NULL, verbose = TRUE) {
  family <- match.arg(family)
  if (family == "beta" && is.null(bounds)) stop("bounds = c(lower, upper) is required for family = 'beta'.")

  ncell <- terra::ncell(prior_value_rasters[[1]])
  use_general <- ncell <= threshold_cells
  route <- if (use_general) "bayesian_update_general" else paste0("closed_form_", family)

  if (verbose) {
    cat(sprintf(
      "[fuse_adaptive] AOI: %d cells (threshold: %d) -> route: %s\n",
      ncell, threshold_cells, route
    ))
  }

  t0 <- proc.time()[["elapsed"]]
  result <- if (use_general) {
    fuse_general_kde(prior_value_rasters, lik_value_rasters, percentile_probs,
                      family, bounds, n_samples, grid_resolution)
  } else {
    fuse_closed_form(prior_value_rasters, lik_value_rasters, percentile_probs, family, bounds)
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

#' @keywords internal
fuse_lognormal_adaptive <- function(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs,
                                     ncell, threshold_cells, n_samples = 500, grid_resolution = NULL, verbose = TRUE) {
  if (ncell <= threshold_cells) {
    # bayesian_update()'s general route makes no distributional assumption
    # at all, so it needs no log-space detour - force fuse_adaptive() into
    # its general route (threshold_cells=Inf) directly on raw values;
    # family="normal" only selects the (mu,sigma) OUTPUT shape.
    r <- fuse_adaptive(prior_value_rasters, lik_value_rasters, prior_probs,
                        family = "normal", threshold_cells = Inf,
                        n_samples = n_samples, grid_resolution = grid_resolution, verbose = verbose)
    return(list(posterior = r$posterior, route = r$route, route_detail = r$route_detail, n_fallback_cells = r$n_fallback_cells))
  }

  if (verbose) {
    cat(sprintf("[fuse_property_adaptive] AOI: %d cells (threshold: %d) -> route: closed_form_lognormal\n", ncell, threshold_cells))
  }
  idx <- percentile_index(prior_probs)
  fit_prior <- fit_normal_raster(prior_value_rasters[[idx$lo_idx]], prior_value_rasters[[idx$p50_idx]], prior_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
  fit_lik <- fit_normal_raster(lik_value_rasters[[idx$lo_idx]], lik_value_rasters[[idx$p50_idx]], lik_value_rasters[[idx$hi_idx]], idx$p_lo, idx$p_hi)
  prior_log <- normal_to_lognormal_params(fit_prior$mu, fit_prior$sigma)
  lik_log <- normal_to_lognormal_params(fit_lik$mu, fit_lik$sigma)
  posterior_log <- bayes_update_normal_normal(prior_log$mu, prior_log$sigma, lik_log$mu, lik_log$sigma)
  posterior <- lognormal_to_normal_params(posterior_log$mu, posterior_log$sigma)
  list(posterior = posterior, route = "closed_form_lognormal", route_detail = NULL, n_fallback_cells = 0)
}

#' Compute a metalog fit's mean/sd via numerical quadrature over its own raw
#' quantile function (`mean = integral of Q(p) dp`, `var = integral of
#' Q(p)^2 dp - mean^2`).
#'
#' @section Known limitation:
#' This is new glue code with no prior validated version (unlike the rest of
#' this file's math, which was validated upstream against `fitdistrplus`/
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
#' @export
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

#' @keywords internal
fuse_metalog_adapter <- function(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs, property_config, verbose = TRUE) {
  interior <- percentile_interior(prior_probs)
  bounds <- property_config$bounds
  boundedness <- property_config$boundedness %||% "b"

  if (verbose) {
    cat("[fuse_property_adaptive] -> route: closed_form_metalog_adapter (always used regardless of AOI size)\n")
  }

  fit_prior <- fit_metalog_linear_raster(prior_value_rasters[interior$idx], interior$probs, bounds, boundedness)
  fit_lik <- fit_metalog_linear_raster(lik_value_rasters[interior$idx], interior$probs, bounds, boundedness)

  infeasible_prior <- check_metalog_feasibility_raster(fit_prior, bounds, boundedness)
  infeasible_lik <- check_metalog_feasibility_raster(fit_lik, bounds, boundedness)

  prior_m <- metalog_moments_raster(fit_prior, infeasible_prior, prior_value_rasters, prior_probs, bounds, boundedness)
  lik_m <- metalog_moments_raster(fit_lik, infeasible_lik, lik_value_rasters, lik_probs, bounds, boundedness)

  posterior <- bayes_update_normal_normal(prior_m$mu, prior_m$sigma, lik_m$mu, lik_m$sigma)
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
#' The raster-native counterpart of `bayesian-updating.R`'s `fuse_property()`,
#' and this toolkit's top-level entry point for non-compositional properties
#' (see `fuse_texture_group()` for the compositional/texture case) - the
#' raster analogue of `run_stage1_fusion()` minus the SSURGO/SOLUS
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
#' @param ... Additional arguments passed to `fuse_adaptive()`/
#'   `fuse_lognormal_adaptive()` (e.g. `n_samples`, `grid_resolution`, `verbose`).
#' @return `list(posterior=, route=, route_detail=, n_fallback_cells=, dist=,
#'   dist_source=, skew_proxy=)`.
#' @export
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
    fuse_metalog_adapter(prior_value_rasters, prior_probs, lik_value_rasters, lik_probs, property_config)
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
#' (see `R/distributions.R`'s `resolve_composition_groups()`) rather than
#' the original source bundle's separate `PROPERTIES`/`config.R` global
#' registry - order is trusted from the config, not hardcoded to literal
#' "clay"/"sand"/"silt" names, exactly like the already-ported
#' `fuse_texture_group_from_triplets()` (`bayesian-updating.R`) already
#' documents ("positional-role placeholders... not an identity requirement").
#'
#' @param group Composition group name (e.g. `"texture"`).
#' @param composition_groups A `config$monte_carlo$composition_groups`-shaped
#'   list, e.g. `list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))`.
#' @return Character vector of member property ids, in configured order.
#' @export
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

#' Fuse a compositional group's members JOINTLY via ILR fusion
#' (`R/distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`,
#' `R/bayesian-updating.R`'s `fuse_bivariate_normal()`), rather than
#' independently via `fuse_beta()` per member - independent fusion
#' measurably breaks sum-to-100 (up to 10.5 percentage points on realistic
#' synthetic data, per upstream validation). The raster counterpart of the
#' already-ported scalar `fuse_texture_group_from_triplets()`.
#'
#' @param fetched A list (one entry per group member, in `group_members()`'s
#'   order) of `list(id=, prior=<named list of ALIGNED percentile-value
#'   rasters>, prior_probs=, lik=<value rasters>, lik_probs=)`.
#' @return A named list keyed by each member's `id`:
#'   `list(posterior = list(value = <inverse-ILR point-estimate raster for
#'   that fraction>, ilr_mu = <2-layer raster, shared across the group>,
#'   ilr_Sigma = <3-layer raster (S11,S12,S22), shared>), dist =
#'   "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL,
#'   n_fallback_cells = 0)`. To draw posterior samples for a specific
#'   fraction/cell, extract that cell's `ilr_mu`/`ilr_Sigma` and pass to
#'   `sample_ilr_posterior()`.
#' @export
fuse_texture_group <- function(fetched) {
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

  fun <- function(row_mat) {
    row_mat <- if (is.matrix(row_mat)) row_mat else matrix(row_mat, nrow = 1, dimnames = list(NULL, col_names))
    t(apply(row_mat, 1, function(row) {
      row <- stats::setNames(as.list(row), col_names)
      prior_m <- estimate_ilr_moments_mc(
        row[[paste0(clay_id, "_prior_lo")]], row[[paste0(clay_id, "_prior_p50")]], row[[paste0(clay_id, "_prior_hi")]],
        row[[paste0(sand_id, "_prior_lo")]], row[[paste0(sand_id, "_prior_p50")]], row[[paste0(sand_id, "_prior_hi")]],
        row[[paste0(silt_id, "_prior_lo")]], row[[paste0(silt_id, "_prior_p50")]], row[[paste0(silt_id, "_prior_hi")]],
        z = prior_z
      )
      lik_m <- estimate_ilr_moments_mc(
        row[[paste0(clay_id, "_lik_lo")]], row[[paste0(clay_id, "_lik_p50")]], row[[paste0(clay_id, "_lik_hi")]],
        row[[paste0(sand_id, "_lik_lo")]], row[[paste0(sand_id, "_lik_p50")]], row[[paste0(sand_id, "_lik_hi")]],
        row[[paste0(silt_id, "_lik_lo")]], row[[paste0(silt_id, "_lik_p50")]], row[[paste0(silt_id, "_lik_hi")]],
        z = lik_z
      )
      fused <- fuse_bivariate_normal(prior_m$mu, prior_m$Sigma, lik_m$mu, lik_m$Sigma)
      c(mu1 = fused$mu[1], mu2 = fused$mu[2],
        S11 = fused$Sigma[1, 1], S12 = fused$Sigma[1, 2], S22 = fused$Sigma[2, 2])
    }))
  }
  result_r <- terra::app(combined, fun)
  names(result_r) <- c("mu1", "mu2", "S11", "S12", "S22")

  mu1_vals <- terra::values(result_r[["mu1"]])[, 1]
  mu2_vals <- terra::values(result_r[["mu2"]])[, 1]
  point_est <- ilr_inverse(mu1_vals, mu2_vals)

  template <- fetched[[1]]$prior[[1]]
  make_layer <- function(v) { r <- template; terra::values(r) <- v; r }

  stats::setNames(lapply(seq_along(ids), function(i) {
    list(
      posterior = list(
        value = make_layer(point_est[, i]),
        ilr_mu = result_r[[c("mu1", "mu2")]],
        ilr_Sigma = result_r[[c("S11", "S12", "S22")]]
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

#' Run Stage 1 Fusion for One Property/Depth over an AOI
#'
#' Fetches the (cached) SSURGO prior (`R/ssurgo-simulation.R`'s `fetch_ssurgo_percentiles()`) and
#' SOLUS likelihood (`R/solus-simulation.R`'s `fetch_solus_percentiles()`) percentiles, aligns the
#' SSURGO grid onto the SOLUS grid via `terra::resample()`, then fuses via
#' `fuse_property_adaptive()`.
#'
#' Compositional properties (`property_config$composition_group` set) are detected up front and
#' delegated to `run_stage1_fusion_group()` instead, since they must be fetched/fit jointly.
#'
#' @param aoi_vect A `terra::SpatVector` AOI (projected, e.g. EPSG:5070).
#' @param property_config A list with `id` (used as `fetch_ssurgo_percentiles()`'s property id and
#'   as the cache key), `solus_variable` (a `soilDB::fetchSOLUS()`-recognized variable name), and
#'   `dist`/`bounds`/`boundedness`/`auto_skew_threshold`/`composition_group` as needed by
#'   `fuse_property_adaptive()`.
#' @param top_depth,bottom_depth Numeric depth bounds in cm.
#' @param composition_groups,property_configs Only needed when `property_config$composition_group`
#'   is set - passed through to `run_stage1_fusion_group()` (see its docs).
#' @param parallel,n_cores Passed through to `simulate_ssurgo_mapunit_draws()`'s
#'   `parallel`/`n_cores` (via `fetch_ssurgo_percentiles()` or `run_stage1_fusion_group()`,
#'   whichever path this call takes) - the SSURGO adapter's per-cokey depth-trend GP fitting is
#'   this pipeline's dominant cost for AOIs with many cokeys. Default `parallel = FALSE` matches
#'   prior behavior exactly.
#' @return `list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=,
#'   route_detail=, n_fallback_cells=)`, or `NULL` if the SSURGO or SOLUS side failed.
#'   `posterior`'s shape depends on `dist` - see `fuse_property_adaptive()`'s docs.
#' @export
run_stage1_fusion <- function(aoi_vect, property_config, top_depth, bottom_depth,
                               composition_groups = NULL, property_configs = NULL,
                               parallel = FALSE, n_cores = NULL) {
  if (!is.null(property_config$composition_group)) {
    group_result <- run_stage1_fusion_group(
      aoi_vect, property_config$composition_group, composition_groups, property_configs,
      top_depth, bottom_depth, parallel = parallel, n_cores = n_cores
    )
    return(if (is.null(group_result)) NULL else group_result[[property_config$id]])
  }

  ssurgo_key <- build_cache_key(aoi_vect, property_config$id, top_depth, bottom_depth, "ssurgo")
  prior <- cache_get(ssurgo_key)
  if (!is.null(prior)) prior <- unwrap_percentile_list(prior)
  if (is.null(prior)) {
    prior <- fetch_ssurgo_percentiles(aoi_vect, property_config$id, top_depth, bottom_depth,
                                       parallel = parallel, n_cores = n_cores)
    if (is.null(prior)) return(NULL)
    cache_set(ssurgo_key, "ssurgo", wrap_percentile_list(prior))
  }

  solus_key <- build_cache_key(aoi_vect, property_config$id, top_depth, bottom_depth, "solus")
  solus <- cache_get(solus_key)
  if (!is.null(solus)) solus <- unwrap_percentile_list(solus)
  if (is.null(solus)) {
    solus <- fetch_solus_percentiles(aoi_vect, property_config$solus_variable, top_depth, bottom_depth)
    if (is.null(solus)) return(NULL)
    cache_set(solus_key, "solus", wrap_percentile_list(solus))
  }

  # SSURGO prior (~30m, from mukey.wcs()) and SOLUS likelihood (100m) are on different grids -
  # resample the prior onto the SOLUS grid before combining, since every fusion route requires
  # cell-aligned rasters.
  reference_grid <- solus$values[[1]]
  prior_aligned <- lapply(prior$values, terra::resample, y = reference_grid, method = "bilinear")

  fused <- fuse_property_adaptive(
    prior_value_rasters = prior_aligned, prior_probs = prior$probs,
    lik_value_rasters = solus$values, lik_probs = solus$probs,
    property_config = property_config
  )

  list(
    prior = list(values = prior_aligned, probs = prior$probs),
    likelihood = list(values = solus$values, probs = solus$probs),
    posterior = fused$posterior,
    dist = fused$dist, dist_source = fused$dist_source, skew_proxy = fused$skew_proxy,
    route = fused$route, route_detail = fused$route_detail, n_fallback_cells = fused$n_fallback_cells
  )
}

#' Run Stage 1 Fusion for a Whole Compositional Group Jointly
#'
#' Fuses a compositional group's members (currently only `"texture"`: clay/sand/silt) jointly via
#' `fuse_texture_group()`, instead of one independent `run_stage1_fusion()` call per member -
#' independent fusion measurably breaks sum-to-100. Each member's SSURGO/SOLUS percentiles are
#' still fetched/cached per-property (so a later request for just one member's raw percentiles
#' still hits cache), but the fusion itself runs once for the whole group, under a `"texture_group"`
#' cache kind, with each member's resulting posterior also seeded into a per-property `"posterior"`
#' cache kind - so three sequential `run_stage1_fusion()` calls for clay, then sand, then silt
#' trigger the joint fetch+fusion exactly once, not three times.
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
#' @return The full group result: a named list keyed by member id, each element
#'   `list(posterior = list(value=, ilr_mu=, ilr_Sigma=), dist = "texture_ilr", route =
#'   "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)` (see
#'   `fuse_texture_group()`'s docs - NOT the uniform `(mu,sigma)`/`(alpha,beta)` contract
#'   non-compositional properties get) - or `NULL` if any member's SSURGO/SOLUS fetch failed.
#'   `run_stage1_fusion()`'s own dispatch slices this down to the single requested member
#'   (`group_result[[property_config$id]]`) to keep its own per-property return contract
#'   consistent regardless of `dist`.
#' @export
run_stage1_fusion_group <- function(aoi_vect, group, composition_groups, property_configs,
                                     top_depth, bottom_depth, parallel = FALSE, n_cores = NULL) {
  member_ids <- group_members(group, composition_groups)
  members <- property_configs[member_ids]

  group_key <- build_cache_key(aoi_vect, group, top_depth, bottom_depth, "texture_group")
  group_result <- cache_get(group_key)
  # cache_set()'s own docs flag that SpatRasters must be wrap()/unwrap()ed around a cache
  # round-trip - group_result nests them inside each member's $posterior$value/ilr_mu/ilr_Sigma,
  # not the simple list(values=, probs=) shape wrap_percentile_list() handles, so the generic
  # wrap_nested_rasters()/unwrap_nested_rasters() pair is needed here instead.
  if (!is.null(group_result)) group_result <- unwrap_nested_rasters(group_result)

  if (is.null(group_result)) {
    # Each member's own SSURGO percentile cache is still checked first (so a later single-
    # property request, or a partially-cached group, still hits cache per member) - but if any
    # member is missing, run the shared mukey-raster fetch + Monte Carlo simulation ONCE here
    # rather than once per member. simulate_cokey_generalized() already simulates every
    # recognized property (clay/sand/silt/db/ph/...) jointly per cokey in one pass, so calling
    # fetch_ssurgo_percentiles() independently per texture member used to re-run that same
    # (expensive - the dominant cost of the whole fusion pipeline) simulation 3 times for a
    # 3-member group, to extract 3 columns that a single simulation pass already produces
    # together. percentiles_from_draws() (R/ssurgo-simulation.R) is the shared quantile/
    # rasterize step, factored out of fetch_ssurgo_percentiles() for exactly this reuse.
    ssurgo_keys <- lapply(members, function(m) build_cache_key(aoi_vect, m$id, top_depth, bottom_depth, "ssurgo"))
    ssurgo_cached <- lapply(ssurgo_keys, function(k) {
      p <- cache_get(k)
      if (!is.null(p)) unwrap_percentile_list(p) else NULL
    })

    shared_mukey_raster <- NULL
    shared_draws <- NULL
    if (any(vapply(ssurgo_cached, is.null, logical(1)))) {
      shared_mukey_raster <- fetch_ssurgo_mukey_raster(aoi_vect)
      if (!is.null(shared_mukey_raster)) {
        shared_draws <- simulate_ssurgo_mapunit_draws(aoi_vect, top_depth, bottom_depth,
                                                        parallel = parallel, n_cores = n_cores)
      }
    }

    fetched <- lapply(seq_along(members), function(idx) {
      m <- members[[idx]]

      prior <- ssurgo_cached[[idx]]
      if (is.null(prior)) {
        prior <- if (is.null(shared_mukey_raster) || is.null(shared_draws)) {
          NULL
        } else {
          percentiles_from_draws(shared_mukey_raster, shared_draws, m$id)
        }
        if (!is.null(prior)) cache_set(ssurgo_keys[[idx]], "ssurgo", wrap_percentile_list(prior))
      }

      solus_key <- build_cache_key(aoi_vect, m$id, top_depth, bottom_depth, "solus")
      solus <- cache_get(solus_key)
      if (!is.null(solus)) solus <- unwrap_percentile_list(solus)
      if (is.null(solus)) {
        solus <- fetch_solus_percentiles(aoi_vect, m$solus_variable, top_depth, bottom_depth)
        if (!is.null(solus)) cache_set(solus_key, "solus", wrap_percentile_list(solus))
      }

      if (is.null(prior) || is.null(solus)) return(NULL)
      prior_aligned <- lapply(prior$values, terra::resample, y = solus$values[[1]], method = "bilinear")
      list(id = m$id, prior = prior_aligned, prior_probs = prior$probs, lik = solus$values, lik_probs = solus$probs)
    })
    if (any(vapply(fetched, is.null, logical(1)))) return(NULL)

    group_result <- fuse_texture_group(fetched)
    cache_set(group_key, "texture_group", wrap_nested_rasters(group_result))

    for (id in member_ids) {
      member_key <- build_cache_key(aoi_vect, id, top_depth, bottom_depth, "posterior")
      cache_set(member_key, "posterior", wrap_nested_rasters(group_result[[id]]))
    }
  }

  group_result
}
