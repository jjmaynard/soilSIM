#' @title Raster SOLUS Percentile Likelihood
#'
#' @description The SOLUS100 half of the raster fusion prior/likelihood pipeline (see
#'   `R/core-fusion.R`): fetches SOLUS100 low/prediction/high rasters via
#'   `soilDB::fetchSOLUS()` for a requested depth window, in the `list(values=, probs=)` shape
#'   `fuse_property_adaptive()` expects.
#'
#'   `fetch_solus_low_pred_high()`'s output-layer-naming assumption
#'   (`paste0(solus_variable, "_", depth_slice, "_cm_", suffix)`) is verified against a live
#'   `soilDB::fetchSOLUS()` call for `output_type` `"prediction"` / `"95% low prediction
#'   interval"` / `"95% high prediction interval"` (suffixes `p`/`l`/`h`).
#'
#'   `solus_variable` is a direct per-call argument (there is no global property registry).
#'   Note that `awc` is not a valid `fetchSOLUS()` variable.
#' @name solus_simulation
NULL

#' Snap a Depth Window to the Nearest Native SOLUS Depth Slice
#'
#' `soilDB::fetchSOLUS()`'s `depth_slices` are single depth points (default `c(0, 5, 15, 30, 60,
#' 100, 150)`, live-verified), not ranges, so a `[top_depth, bottom_depth]` window must be
#' approximated. This function snaps to the native slice closest to the window's midpoint.
#'
#' @section Relationship to the fusion pipeline:
#' `fetch_solus_low_pred_high()` uses `solus_depth_window_weights()` for a trapezoidal
#' depth-average across the native slices spanning the window, so the SOLUS likelihood and the
#' SSURGO prior both represent a thickness-weighted window mean. This function is public API for
#' callers that want a single nearest slice.
#'
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param available_slices Native SOLUS depth points to snap to.
#' @return A single numeric value from `available_slices`.
#' @seealso `solus_depth_window_weights()`
#' @keywords internal
closest_solus_depth_slice <- function(top_depth, bottom_depth, available_slices = c(0, 5, 15, 30, 60, 100, 150)) {
  midpoint <- (top_depth + bottom_depth) / 2
  available_slices[which.min(abs(available_slices - midpoint))]
}

#' Trapezoidal Depth-Average Weights over Native SOLUS Depth Slices
#'
#' Returns the weights that turn SOLUS100's point predictions at the native depth slices into an
#' estimate of the **depth-average** of the property over `[top_depth, bottom_depth]` -
#' `sum(weights * slice_values)` approximates `(1 / (b - a)) * integral_a^b f(d) dd`, where `f` is
#' taken piecewise-linear between native slices (with `rule = 2` constant extrapolation outside the
#' `[0, 150]` span SOLUS predicts over). This is the raster analogue of the thickness-weighted
#' window mean the SSURGO prior side already computes (`aggregate_depth_window_by_replicate()`), so
#' the two sides of `fuse_property_adaptive()` represent the same estimand.
#'
#' Exact for a linear depth profile (reduces to the midpoint value, matching
#' `closest_solus_depth_slice()`); unbiased to second order for curved profiles (OC decay, clay
#' bulge, pH trend). Cost is a single weighted raster-stack sum - `length(available_slices)` is
#' always small, so this scales with the number of slices, not the number of pixels.
#'
#' Windows lying entirely at/beyond the deepest slice (or above the shallowest) fall back to the
#' single nearest slice (constant extrapolation).
#'
#' @param top_depth,bottom_depth Numeric depth window bounds in cm (`bottom_depth > top_depth`).
#' @param available_slices Native SOLUS depth points (default `c(0, 5, 15, 30, 60, 100, 150)`).
#' @return A named numeric vector over `sort(unique(available_slices))` (names are the slice depths
#'   as character), summing to 1. Zero for slices that do not contribute.
#' @seealso `closest_solus_depth_slice()`, `fetch_solus_low_pred_high()`
#' @keywords internal
solus_depth_window_weights <- function(top_depth, bottom_depth,
                                       available_slices = c(0, 5, 15, 30, 60, 100, 150)) {
  if (!is.finite(top_depth) || !is.finite(bottom_depth) || bottom_depth <= top_depth) {
    stop("solus_depth_window_weights(): need finite bounds with bottom_depth > top_depth.")
  }
  slices <- sort(unique(available_slices))
  n <- length(slices)
  smin <- slices[1]
  smax <- slices[n]
  w <- stats::setNames(numeric(n), as.character(slices))

  a <- max(top_depth, smin)
  b <- min(bottom_depth, smax)
  if (b <= a) {
    # Window entirely outside the predicted span - nearest-slice constant extrapolation.
    w[which.min(abs(slices - (top_depth + bottom_depth) / 2))] <- 1
    return(w)
  }

  # Weights on `slices` reproducing f(d) by piecewise-linear interpolation, rule = 2 outside span.
  interp_w <- function(d) {
    v <- numeric(n)
    if (d <= smin) { v[1] <- 1; return(v) }
    if (d >= smax) { v[n] <- 1; return(v) }
    j <- min(max(findInterval(d, slices), 1L), n - 1L)
    t <- (d - slices[j]) / (slices[j + 1L] - slices[j])
    v[j] <- 1 - t
    v[j + 1L] <- t
    v
  }

  # Trapezoidal knots: window ends plus every native slice strictly inside the window.
  knots <- sort(unique(c(a, slices[slices > a & slices < b], b)))
  acc <- numeric(n)
  for (k in seq_len(length(knots) - 1L)) {
    lo <- knots[k]
    hi <- knots[k + 1L]
    acc <- acc + ((hi - lo) / 2) * (interp_w(lo) + interp_w(hi))
  }
  w[] <- acc / (b - a)
  w
}

#' Fetch SOLUS100 Low/Prediction/High Rasters for One Variable and Depth Window
#'
#' Fetches every native SOLUS depth slice spanning `[top_depth, bottom_depth]` and reduces them to
#' a single **depth-average** raster per output type via `solus_depth_window_weights()` - a
#' trapezoidal integral of SOLUS's point predictions over the window. See
#' `solus_depth_window_weights()` for why this replaces the previous single-nearest-slice snap
#' (`closest_solus_depth_slice()`): a point value only represents a window average for a linear
#' depth profile, and the SSURGO prior side is a genuine thickness-weighted window mean.
#'
#' The same window weights are applied to the `prediction` raster and both `95% ... prediction
#' interval` rasters. Averaging the interval bounds treats the prediction interval as perfectly
#' depth-correlated (comonotone) - a mildly conservative (wider) choice for a likelihood, and the
#' opposite assumption (depth-independence, narrow by ~sqrt(n)) is equally wrong the other way.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variable A `soilDB::fetchSOLUS()`-recognized variable name (e.g. `"claytotal"`,
#'   `"dbovendry"`, `"ph1to1h2o"`).
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return `list(pred=, low=, high=)`, each a single-layer `terra::SpatRaster` (the window
#'   depth-average) or `NULL` if that output type wasn't returned for every required slice.
#' @keywords internal
fetch_solus_low_pred_high <- function(aoi_vect, solus_variable, top_depth, bottom_depth) {
  weights <- solus_depth_window_weights(top_depth, bottom_depth)
  weights <- weights[weights != 0]
  slice_depths <- names(weights)
  layer_name <- function(depth, suffix) paste0(solus_variable, "_", depth, "_cm_", suffix)

  fetch_one <- function(output_type, suffix) {
    result <- tryCatch(
      soilDB::fetchSOLUS(aoi_vect, depth_slices = slice_depths,
                          variables = solus_variable, output_type = output_type, grid = TRUE),
      error = function(e) {
        warning(sprintf("fetch_solus_low_pred_high(): fetchSOLUS() failed for output_type '%s': %s",
                         output_type, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(result)) return(NULL)

    wanted <- vapply(slice_depths, layer_name, character(1), suffix = suffix)
    if (!all(wanted %in% names(result))) {
      warning(sprintf("fetch_solus_low_pred_high(): fetchSOLUS() did not return every required depth slice for output_type '%s' (missing %s).",
                       output_type, paste(setdiff(wanted, names(result)), collapse = ", ")))
      return(NULL)
    }

    stack <- result[[wanted]]
    # Weighted stack sum = the trapezoidal depth-average. na.rm = FALSE: an NA at any contributing
    # slice (edge of coverage) propagates rather than silently dropping depth mass.
    terra::app(stack * unname(weights), fun = "sum", na.rm = FALSE)
  }

  list(
    pred = fetch_one("prediction", "p"),
    low = fetch_one("95% low prediction interval", "l"),
    high = fetch_one("95% high prediction interval", "h")
  )
}

#' Fetch SOLUS100 Percentile-Value Rasters for an AOI
#'
#' The top-level SOLUS "likelihood" entry point for `R/core-fusion.R`'s
#' `fuse_property_adaptive()`.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variable A `soilDB::fetchSOLUS()`-recognized variable name.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return `list(values = list(P025 = <low raster>, P50 = <pred raster>, P975 = <high raster>),
#'   probs = c(0.025, 0.5, 0.975))`, or `NULL` if any of low/pred/high is unavailable.
#' @keywords internal
fetch_solus_percentiles <- function(aoi_vect, solus_variable, top_depth, bottom_depth) {
  lph <- fetch_solus_low_pred_high(aoi_vect, solus_variable, top_depth, bottom_depth)

  if (is.null(lph$pred) || is.null(lph$low) || is.null(lph$high)) {
    return(NULL)
  }

  list(
    values = list(P025 = lph$low, P50 = lph$pred, P975 = lph$high),
    probs = c(0.025, 0.5, 0.975)
  )
}

#' Fetch SOLUS100 Low/Prediction/High Rasters for Multiple Variables, One Depth Window (S1)
#'
#' Batched sibling of `fetch_solus_low_pred_high()`: `soilDB::fetchSOLUS()` accepts `variables`,
#' `depth_slices`, and `output_type` all as vectors simultaneously (live-verified 2026-09-04 - see
#'), so every variable's low/prediction/high
#' rasters for a window can be fetched in **one** network call instead of `3 * length(solus_variables)`
#' separate ones. The returned layer-naming convention (`paste0(variable, "_", depth, "_cm_",
#' suffix)`) was confirmed identical to the single-variable case for the fully-batched request
#' (`nlyr` matched `length(variables) * length(depth_slices) * length(output_type)` exactly, no
#' layers dropped, `names()` order not request order but that does not matter since layers are
#' looked up by exact name).
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variables Character vector of `soilDB::fetchSOLUS()`-recognized variable names.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return A named list, one entry per element of `solus_variables`, each `list(pred=, low=, high=)`
#'   in the same shape `fetch_solus_low_pred_high()` returns (a single-layer `terra::SpatRaster` or
#'   `NULL` per output type). If the single underlying `fetchSOLUS()` call itself errors, every
#'   variable's entry is `list(pred=NULL, low=NULL, high=NULL)` (a whole-request failure, not a
#'   crash) - callers should fall back to the scalar `fetch_solus_low_pred_high()` per variable in
#'   that case. A per-`(variable, output_type)` combo missing from an otherwise-successful response
#'   is `NULL` for just that piece, with a `warning()` naming it - a partial failure, not a whole-
#'   request one.
#' @seealso `fetch_solus_low_pred_high()`, `fetch_solus_percentiles_multiproperty()`
#' @keywords internal
fetch_solus_low_pred_high_multiproperty <- function(aoi_vect, solus_variables, top_depth, bottom_depth) {
  weights <- solus_depth_window_weights(top_depth, bottom_depth)
  weights <- weights[weights != 0]
  slice_depths <- names(weights)
  layer_name <- function(var, depth, suffix) paste0(var, "_", depth, "_cm_", suffix)

  empty_result <- stats::setNames(
    rep(list(list(pred = NULL, low = NULL, high = NULL)), length(solus_variables)),
    solus_variables
  )

  result <- tryCatch(
    soilDB::fetchSOLUS(aoi_vect, depth_slices = slice_depths, variables = solus_variables,
                        output_type = c("prediction", "95% low prediction interval",
                                        "95% high prediction interval"), grid = TRUE),
    error = function(e) {
      warning(sprintf("fetch_solus_low_pred_high_multiproperty(): fetchSOLUS() failed: %s",
                       conditionMessage(e)))
      NULL
    }
  )
  if (is.null(result)) return(empty_result)

  fetch_one <- function(var, suffix, output_type) {
    wanted <- vapply(slice_depths, layer_name, character(1), var = var, suffix = suffix)
    if (!all(wanted %in% names(result))) {
      warning(sprintf("fetch_solus_low_pred_high_multiproperty(): fetchSOLUS() did not return every required depth slice for variable '%s', output_type '%s' (missing %s).",
                       var, output_type, paste(setdiff(wanted, names(result)), collapse = ", ")))
      return(NULL)
    }
    stack <- result[[wanted]]
    terra::app(stack * unname(weights), fun = "sum", na.rm = FALSE)
  }

  stats::setNames(lapply(solus_variables, function(v) {
    list(
      pred = fetch_one(v, "p", "prediction"),
      low = fetch_one(v, "l", "95% low prediction interval"),
      high = fetch_one(v, "h", "95% high prediction interval")
    )
  }), solus_variables)
}

#' Fetch SOLUS100 Percentile-Value Rasters for Multiple Variables, One Depth Window (S1)
#'
#' Batched sibling of `fetch_solus_percentiles()`: thin per-variable `list(values=, probs=)`
#' packaging over `fetch_solus_low_pred_high_multiproperty()`'s one-call fetch.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variables Character vector of `soilDB::fetchSOLUS()`-recognized variable names.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return A named list, one entry per element of `solus_variables`, each either
#'   `list(values = list(P025=, P50=, P975=), probs = c(0.025, 0.5, 0.975))` (matching
#'   `fetch_solus_percentiles()`'s shape) or `NULL` if any of that variable's low/pred/high rasters
#'   is unavailable.
#' @seealso `fetch_solus_percentiles()`, `fetch_solus_low_pred_high_multiproperty()`
#' @keywords internal
fetch_solus_percentiles_multiproperty <- function(aoi_vect, solus_variables, top_depth, bottom_depth) {
  lph_by_var <- fetch_solus_low_pred_high_multiproperty(aoi_vect, solus_variables, top_depth, bottom_depth)

  stats::setNames(lapply(solus_variables, function(v) {
    lph <- lph_by_var[[v]]
    if (is.null(lph$pred) || is.null(lph$low) || is.null(lph$high)) return(NULL)
    list(
      values = list(P025 = lph$low, P50 = lph$pred, P975 = lph$high),
      probs = c(0.025, 0.5, 0.975)
    )
  }), solus_variables)
}

#' Fetch SOLUS100 Rasters for a Site-Level (Depth-Independent) Variable (S2)
#'
#' Some `soilDB::fetchSOLUS()` variables are **site-level** - a single prediction per pixel with no
#' depth dependence (e.g. `anylithicdpt`/`resdept`, depth to bedrock/restriction) - requested via
#' the special `depth_slices = "all"` value rather than a numeric depth. Confirmed live during the
#' S1 spike (2026-09-04): the returned layer is named `paste0(variable, "_all_cm_", suffix)`,
#' exactly one layer per output type, so unlike `fetch_solus_low_pred_high()` there is no
#' depth-window weighting to do - deliberately **not** built on `solus_depth_window_weights()`,
#' which only makes sense for numeric depth slices.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param variable A `soilDB::fetchSOLUS()`-recognized site-level variable name (e.g.
#'   `"anylithicdpt"`, `"resdept"`).
#' @param output_types Which output types to fetch - any of `"prediction"`, `"95% low prediction
#'   interval"`, `"95% high prediction interval"`. Default fetches all three.
#' @return `list(pred=, low=, high=)`, each a single-layer `terra::SpatRaster` (or `NULL` for an
#'   output type not requested, or not returned by `fetchSOLUS()`).
#' @seealso `fetch_solus_restriction_depth()`
#' @keywords internal
fetch_solus_site_level <- function(aoi_vect, variable,
                                   output_types = c("prediction", "95% low prediction interval",
                                                    "95% high prediction interval")) {
  suffix_for <- c("prediction" = "p", "95% low prediction interval" = "l",
                  "95% high prediction interval" = "h")

  result <- tryCatch(
    soilDB::fetchSOLUS(aoi_vect, depth_slices = "all", variables = variable,
                       output_type = output_types, grid = TRUE),
    error = function(e) {
      warning(sprintf("fetch_solus_site_level(): fetchSOLUS() failed for variable '%s': %s",
                       variable, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(result)) return(list(pred = NULL, low = NULL, high = NULL))

  fetch_one <- function(output_type) {
    if (!output_type %in% output_types) return(NULL)
    wanted <- paste0(variable, "_all_cm_", suffix_for[[output_type]])
    if (!wanted %in% names(result)) {
      warning(sprintf("fetch_solus_site_level(): fetchSOLUS() did not return the expected layer '%s' for variable '%s', output_type '%s'.",
                       wanted, variable, output_type))
      return(NULL)
    }
    result[[wanted]]
  }

  list(
    pred = fetch_one("prediction"),
    low = fetch_one("95% low prediction interval"),
    high = fetch_one("95% high prediction interval")
  )
}

#' SOLUS100's Right-Censoring Depth for `anylithicdpt`/`resdept` (S2)
#'
#' Both depth-to-bedrock (`anylithicdpt`) and depth-to-restriction (`resdept`) are right-censored
#' at 201 cm in SOLUS100's own training data, per gSSURGO data-definition convention (Nauman et al.
#' 2024, SOLUS100 documentation section 2.1.3, confirmed against the actual training-data-prep
#' source). A predicted value at/above this censor point means "no restriction detected within the
#' training depth range" - not a literal claim that bedrock/restriction sits at exactly 201 cm.
#' `fetch_solus_restriction_depth()` recodes values `>= SOLUS_RESTRICTION_CENSOR_CM` to `Inf`
#' (never truncates) rather than treating 201 as a real cutoff.
#' @keywords internal
SOLUS_RESTRICTION_CENSOR_CM <- 201

#' Fetch a SOLUS100 Restriction-Depth Raster, Right-Censoring-Guarded (S2)
#'
#' The truncation-depth signal for bedrock/restriction-aware AWC (`remarginalized_awc()`'s
#' `restriction_depth` argument). Default variable is
#' `"anylithicdpt"` - trained (per SOLUS100's own published methodology) **only** on
#' `reskind %in% c("Lithic bedrock", "Paralithic bedrock")` records, making it hard-bedrock-specific
#' by construction; no separate SSURGO-side hard/soft classifier is needed since SOLUS already
#' encodes that distinction across its two site-level variables. `"resdept"` is trained on *every*
#' restriction record regardless of kind (fragipans, duripans, petrocalcic, etc. included) - a
#' broader, more conservative (more truncation-prone) alternative, offered but not the default,
#' since it would truncate through partially-penetrable layers with no way to tell those apart from
#' true bedrock.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param variable `"anylithicdpt"` (default, hard-bedrock-specific) or `"resdept"` (broader - any
#'   restriction, hard or soft).
#' @return A single-layer `terra::SpatRaster` of restriction depth in cm, with every cell
#'   `>= SOLUS_RESTRICTION_CENSOR_CM` recoded to `Inf` (no restriction detected - never truncates),
#'   or `NULL` if the underlying `fetchSOLUS()` prediction fetch fails.
#' @seealso `SOLUS_RESTRICTION_CENSOR_CM`, `fetch_solus_site_level()`,
#'   [remarginalized_awc()][remarginalized_awc]'s `restriction_depth` argument
#' @keywords internal
fetch_solus_restriction_depth <- function(aoi_vect, variable = "anylithicdpt") {
  lph <- fetch_solus_site_level(aoi_vect, variable, output_types = "prediction")
  if (is.null(lph$pred)) return(NULL)
  terra::ifel(lph$pred >= SOLUS_RESTRICTION_CENSOR_CM, Inf, lph$pred)
}
