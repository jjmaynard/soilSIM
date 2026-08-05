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
#' approximated by one native slice. This function did not exist anywhere in the source this was
#' ported from (referenced but never defined) - genuinely new design, not a port: snaps to the
#' native slice closest to the window's midpoint, the simplest representative-point choice for an
#' otherwise point-based data source.
#'
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param available_slices Native SOLUS depth points to snap to.
#' @return A single numeric value from `available_slices`.
#' @export
closest_solus_depth_slice <- function(top_depth, bottom_depth, available_slices = c(0, 5, 15, 30, 60, 100, 150)) {
  midpoint <- (top_depth + bottom_depth) / 2
  available_slices[which.min(abs(available_slices - midpoint))]
}

#' Fetch SOLUS100 Low/Prediction/High Rasters for One Variable and Depth Window
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param solus_variable A `soilDB::fetchSOLUS()`-recognized variable name (e.g. `"claytotal"`,
#'   `"dbovendry"`, `"ph1to1h2o"`).
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @return `list(pred=, low=, high=)`, each a single-layer `terra::SpatRaster` or `NULL` if that
#'   output type wasn't returned.
#' @export
fetch_solus_low_pred_high <- function(aoi_vect, solus_variable, top_depth, bottom_depth) {
  depth_slice <- closest_solus_depth_slice(top_depth, bottom_depth)
  layer_name <- function(suffix) paste0(solus_variable, "_", depth_slice, "_cm_", suffix)

  fetch_one <- function(output_type, suffix) {
    result <- tryCatch(
      soilDB::fetchSOLUS(aoi_vect, depth_slices = as.character(depth_slice),
                          variables = solus_variable, output_type = output_type, grid = TRUE),
      error = function(e) {
        warning(sprintf("fetch_solus_low_pred_high(): fetchSOLUS() failed for output_type '%s': %s",
                         output_type, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(result) || !(layer_name(suffix) %in% names(result))) {
      return(NULL)
    }
    result[[layer_name(suffix)]]
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
