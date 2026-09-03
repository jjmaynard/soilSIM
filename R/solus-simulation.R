#' @title Raster SOLUS Percentile Likelihood
#'
#' @description The SOLUS100 half of the raster fusion prior/likelihood pipeline (see
#'   `R/raster-fusion.R`): fetches SOLUS100 low/prediction/high rasters via
#'   `soilDB::fetchSOLUS()` for a requested depth window, in the `list(values=, probs=)` shape
#'   `fuse_property_adaptive()` expects. Ported from `code_ref/reanalysis-platform/solus_raster.R`.
#'
#'   `fetch_solus_low_pred_high()`'s output-layer-naming assumption
#'   (`paste0(solus_variable, "_", depth_slice, "_cm_", suffix)`) was verified against a live,
#'   network-connected `soilDB::fetchSOLUS()` call before porting (not assumed from the source
#'   comments alone) - confirmed exact for `output_type` `"prediction"`/`"95% low prediction
#'   interval"`/`"95% high prediction interval"` (suffixes `p`/`l`/`h`).
#'
#'   Unlike the original bundle's `config.R`-driven `solus_variable` lookup (which had a
#'   confirmed-wrong `awc = "awc"` entry - not a valid `fetchSOLUS()` variable), `solus_variable`
#'   is a direct per-call argument here, consistent with `R/raster-fusion.R`'s decision not to
#'   port the `PROPERTIES`/`config.R` global registry.
#' @name solus_simulation
NULL

#' Snap a Depth Window to the Nearest Native SOLUS Depth Slice
#'
#' `soilDB::fetchSOLUS()`'s `depth_slices` are single depth points (default `c(0, 5, 15, 30, 60,
#' 100, 150)`, live-verified), not ranges, so a `[top_depth, bottom_depth]` window must be
#' approximated. This function did not exist anywhere in the source this was ported from
#' (referenced but never defined) - genuinely new design, not a port: snaps to the native slice
#' closest to the window's midpoint.
#'
#' @section Superseded for the fusion pipeline:
#' `fetch_solus_low_pred_high()` no longer uses this - a single midpoint point-value only equals
#' the window depth-average when the property is linear in depth, which most soil properties are
#' not, and it left the SOLUS likelihood representing a different quantity (a point) than the
#' SSURGO prior side (a genuine thickness-weighted window mean, see
#' `aggregate_depth_window_by_replicate()`). It now uses `solus_depth_window_weights()` for a
#' trapezoidal depth-average across the native slices spanning the window. This function is
#' retained as public API for back-compat and for callers that genuinely want a single nearest
#' slice.
#'
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param available_slices Native SOLUS depth points to snap to.
#' @return A single numeric value from `available_slices`.
#' @seealso `solus_depth_window_weights()`
#' @export
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
#' @export
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
#' @export
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
#' The top-level SOLUS "likelihood" entry point for `R/raster-fusion.R`'s
#' `fuse_property_adaptive()`.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variable A `soilDB::fetchSOLUS()`-recognized variable name.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return `list(values = list(P025 = <low raster>, P50 = <pred raster>, P975 = <high raster>),
#'   probs = c(0.025, 0.5, 0.975))`, or `NULL` if any of low/pred/high is unavailable.
#' @export
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
