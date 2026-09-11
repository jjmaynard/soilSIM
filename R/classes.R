# S3 result classes for soilSIM's user-facing engines.
#
# Each constructor is a thin wrapper: it attaches a class to the plain list the
# engine already returns, so existing `$` / `[[` / `names()` access is unchanged.
# `print`/`summary` give a readable overview instead of dumping a deeply nested
# list. See CONTRIBUTING.md ("Result objects").

# --- shared helpers --------------------------------------------------------

.new_soilSIM <- function(x, subclass) {
  if (!is.list(x)) {
    stop("soilSIM result objects wrap a list; got ", class(x)[1], call. = FALSE)
  }
  class(x) <- unique(c(subclass, "soilSIM_result", "list"))
  x
}

#' @export
print.soilSIM_result <- function(x, ...) {
  cat("<", class(x)[1], ">\n", sep = "")
  cat("Components:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.soilSIM_result <- function(object, ...) {
  print(object, ...)
}

.fmt_num <- function(x, digits = 3) {
  if (is.null(x) || !is.numeric(x) || length(x) != 1 || is.na(x)) return("NA")
  formatC(x, digits = digits, format = "g")
}

# --- soilSIM_simulation ---------------------------------------------------

#' Monte Carlo simulation result
#'
#' Wraps the list returned by [simulate_monte_carlo()]. Nothing
#' about the contents changes; the class only adds `print`/`summary`.
#'
#' @param x A list produced by the Monte Carlo engine.
#' @return An object of class `soilSIM_simulation`.
#' @export
new_soilSIM_simulation <- function(x) .new_soilSIM(x, "soilSIM_simulation")

#' @export
print.soilSIM_simulation <- function(x, ...) {
  md <- x$metadata %||% list()
  cat("<soilSIM_simulation>\n")
  cat("  properties     :", paste(md$properties %||% character(0), collapse = ", "), "\n")
  cat("  realizations   :", md$n_realizations %||% NA, "\n")
  cat("  horizons       :", md$n_horizons %||% NA, "\n")
  if (!is.null(md$success_rate)) {
    cat("  success rate   :", .fmt_num(md$success_rate * 100, 4), "%\n")
  }
  qs <- x$quality_assessment$overall_quality_score
  if (!is.null(qs)) cat("  quality score  :", .fmt_num(qs), "\n")
  if (!is.null(md$processing_time)) {
    cat("  processing time:", .fmt_num(as.numeric(md$processing_time), 3), "s\n")
  }
  invisible(x)
}

#' @export
summary.soilSIM_simulation <- function(object, ...) {
  print(object, ...)
  sd <- object$simulation_data
  if (!is.null(dim(sd))) {
    cat("  simulation_data dim:", paste(dim(sd), collapse = " x "), "\n")
  }
  invisible(object)
}

#' @export
as.data.frame.soilSIM_simulation <- function(x, ...) {
  sd <- x$simulation_data
  if (is.data.frame(sd)) return(sd)
  as.data.frame(sd, ...)
}

# --- soilSIM_statistics -------------------------------------------------

#' Descriptive-statistics result
#'
#' Wraps the list returned by [analyze_soil_statistics()].
#'
#' @param x A list produced by [analyze_soil_statistics()].
#' @return An object of class `soilSIM_statistics`.
#' @export
new_soilSIM_statistics <- function(x) .new_soilSIM(x, "soilSIM_statistics")

#' @export
print.soilSIM_statistics <- function(x, ...) {
  md <- x$analysis_metadata %||% list()
  cat("<soilSIM_statistics>\n")
  cat("  properties analyzed:", md$properties_analyzed %||% NA, "\n")
  qs <- x$quality_report$overall_quality_score
  if (!is.null(qs)) cat("  quality score      :", .fmt_num(qs), "\n")
  cm <- x$correlation_matrices$matrices
  if (!is.null(cm)) cat("  correlation groups :", length(cm), "\n")
  invisible(x)
}

# --- soilSIM_gp_models -------------------------------------------------

#' Stratified GP depth-model set
#'
#' Wraps the list returned by [fit_depth_gp_models()].
#'
#' @param x A list produced by [fit_depth_gp_models()].
#' @return An object of class `soilSIM_gp_models`.
#' @export
new_soilSIM_gp_models <- function(x) .new_soilSIM(x, "soilSIM_gp_models")

#' @export
print.soilSIM_gp_models <- function(x, ...) {
  props <- setdiff(names(x), "model_summary")
  cat("<soilSIM_gp_models>\n")
  cat("  properties:", paste(props, collapse = ", "), "\n")
  for (p in props) {
    ng <- x[[p]]$n_groups %||% length(x[[p]]$models %||% list())
    cat("   -", p, ":", ng, "group model(s)\n")
  }
  invisible(x)
}

# --- soilSIM_fusion --------------------------------------------------

#' Multi-source raster fusion result
#'
#' Wraps the list returned by [run_fusion()] / [run_fusion_group()].
#'
#' @param x A list produced by a fusion orchestrator.
#' @return An object of class `soilSIM_fusion`.
#' @export
new_soilSIM_fusion <- function(x) .new_soilSIM(x, "soilSIM_fusion")

# --- soilSIM_diagnostics --------------------------------------------

#' Diagnostic report
#'
#' Wraps the list a `diagnose_*()` / master `validate_*()` reporter returns.
#'
#' @param x A diagnostic-report list.
#' @return An object of class `soilSIM_diagnostics`.
#' @export
new_soilSIM_diagnostics <- function(x) .new_soilSIM(x, "soilSIM_diagnostics")

#' @export
print.soilSIM_diagnostics <- function(x, ...) {
  cat("<soilSIM_diagnostics>\n")
  oa <- x$overall_assessment %||% list()
  status <- x$overall_status %||% x$status %||% oa$workflow_status
  qs <- x$overall_quality_score %||% oa$quality_score
  if (!is.null(status)) cat("  status :", status, "\n")
  if (!is.null(qs)) cat("  quality score:", .fmt_num(qs),
                       if (!is.null(oa$quality_grade)) paste0(" (", oa$quality_grade, ")") else "", "\n")
  cat("  sections:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.soilSIM_diagnostics <- function(object, ...) {
  print(object, ...)
  oa <- object$overall_assessment %||% list()
  cs <- oa$component_scores
  if (!is.null(cs) && length(cs)) {
    cat("  component scores:\n")
    for (nm in names(cs)) cat("   -", nm, ":", .fmt_num(cs[[nm]]), "\n")
  }
  ci <- oa$critical_issues
  if (!is.null(ci) && length(ci)) {
    cat("  critical issues:", length(ci), "\n")
  }
  rec <- object$recommendations
  if (!is.null(rec) && length(rec)) {
    cat("  recommendations:\n")
    for (r in utils::head(rec, 10)) cat("   -", r, "\n")
  }
  invisible(object)
}

# --- soilSIM_property_config ---------------------------------------

#' Property configuration object
#'
#' Wraps a property-configuration list ([default_property_config()],
#' [create_custom_property_config()]).
#'
#' @param x A property-configuration list.
#' @return An object of class `soilSIM_property_config`.
#' @export
new_soilSIM_property_config <- function(x) .new_soilSIM(x, "soilSIM_property_config")

#' @export
print.soilSIM_property_config <- function(x, ...) {
  props <- x$properties %||% x
  cat("<soilSIM_property_config>\n")
  cat("  properties:", if (is.list(props)) length(props) else length(props), "\n")
  nm <- if (!is.null(names(props))) names(props) else NULL
  if (!is.null(nm)) cat("  ", paste(utils::head(nm, 20), collapse = ", "),
                        if (length(nm) > 20) " ..." else "", "\n", sep = "")
  invisible(x)
}

# --- predict/plot methods (D6.4) -------------------------------------------
#
# ggplot2 stays a Suggests dependency: every plot method below checks
# requireNamespace("ggplot2") and degrades to a message + invisible(NULL) if it
# isn't installed, rather than promoting ggplot2 to Imports or reimplementing
# in base graphics. plot.soilSIM_fusion is the one exception - its panels are
# terra::SpatRaster objects, which terra::plot() (already an Imports dependency)
# draws directly; ggplot2 has no native SpatRaster geom without the optional
# tidyterra package, so base/terra graphics is the correct tool there, not a
# deviation for its own sake.

#' Predict a fitted GP depth-trend model at new depths
#'
#' A convenience wrapper over [predict_gp_depth_trends()] (unchanged - this does not replace it,
#' since that function's real callers throughout the package operate on one leaf model at a time,
#' not the whole `soilSIM_gp_models` collection) that looks up the requested property/group's
#' fitted model inside a [fit_depth_gp_models()] result.
#'
#' @param object A `soilSIM_gp_models` object from [fit_depth_gp_models()].
#' @param new_depths Numeric vector of depths (cm) to predict at.
#' @param property Which property's model to use (a name in `object`).
#' @param group Which stratified group's model to use (default: the first group fitted for
#'   `property`, with a message when more than one exists - pass explicitly to silence it).
#' @param ... Passed to [predict_gp_depth_trends()].
#' @return Numeric vector of predicted means, one per `new_depths` (see
#'   [predict_gp_depth_trends()]'s own `@return` - it is the fitted mean only, not a prediction
#'   interval).
#' @family gp-modeling
#' @export
predict.soilSIM_gp_models <- function(object, new_depths, property, group = NULL, ...) {
  if (!property %in% setdiff(names(object), "model_summary")) {
    stop("predict.soilSIM_gp_models(): '", property, "' not found. Available: ",
         paste(setdiff(names(object), "model_summary"), collapse = ", "))
  }
  prop_entry <- object[[property]]
  groups <- names(prop_entry$models)
  if (is.null(group)) {
    group <- groups[1]
    if (length(groups) > 1) {
      message("predict.soilSIM_gp_models(): multiple groups fitted for '", property,
              "' (", paste(groups, collapse = ", "), ") - using '", group,
              "'; pass group = to choose another.")
    }
  } else if (!group %in% groups) {
    stop("predict.soilSIM_gp_models(): group '", group, "' not found for '", property,
         "'. Available: ", paste(groups, collapse = ", "))
  }
  predict_gp_depth_trends(prop_entry$models[[group]], new_depths, ...)
}

#' @export
plot.soilSIM_simulation <- function(x, property = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Install the 'ggplot2' package to use plot.soilSIM_simulation().")
    return(invisible(NULL))
  }
  sim <- x$simulation_data
  if (is.null(sim) || is.null(dim(sim))) {
    message("plot.soilSIM_simulation(): no simulation_data array found on this object.")
    return(invisible(NULL))
  }
  props <- dimnames(sim)[[2]]
  property <- if (is.null(property)) props[1] else match.arg(property, props)
  df <- data.frame(value = as.vector(sim[, property, ]))
  print(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(bins = 30, fill = "#4C72B0", colour = "white") +
      ggplot2::labs(title = paste("Simulated distribution:", property), x = property, y = "count") +
      ggplot2::theme_minimal()
  )
  invisible(x)
}

#' @export
plot.soilSIM_gp_models <- function(x, property = NULL, group = NULL, new_depths = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Install the 'ggplot2' package to use plot.soilSIM_gp_models().")
    return(invisible(NULL))
  }
  props <- setdiff(names(x), "model_summary")
  property <- if (is.null(property)) props[1] else match.arg(property, props)
  if (is.null(new_depths)) new_depths <- seq(0, 150, by = 5)
  preds <- predict.soilSIM_gp_models(x, new_depths, property = property, group = group)
  df <- data.frame(depth = new_depths, mean = preds)
  print(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$depth, y = .data$mean)) +
      ggplot2::geom_line(colour = "#4C72B0", linewidth = 1) +
      ggplot2::labs(title = paste("GP depth trend:", property),
                    x = "Depth (cm)", y = property,
                    caption = "Fitted mean only - predict_gp_depth_trends() does not return a prediction interval.") +
      ggplot2::theme_minimal()
  )
  invisible(x)
}

#' @export
plot.soilSIM_fusion <- function(x, ...) {
  pick_mid <- function(values, probs) {
    if (is.null(values) || !length(values)) return(NULL)
    if (is.null(probs)) return(values[[1]])
    values[[which.min(abs(probs - 0.5))]]
  }
  posterior_raster <- function(posterior) {
    if (is.null(posterior)) return(NULL)
    for (nm in c("mu", "value", "alpha", "shape")) {
      if (!is.null(posterior[[nm]]) && inherits(posterior[[nm]], "SpatRaster")) return(posterior[[nm]])
    }
    hit <- Filter(function(v) inherits(v, "SpatRaster"), posterior)
    if (length(hit)) hit[[1]] else NULL
  }
  panels <- list(
    Prior = pick_mid(x$prior$values, x$prior$probs),
    Likelihood = pick_mid(x$likelihood$values, x$likelihood$probs),
    Posterior = posterior_raster(x$posterior)
  )
  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) {
    message("plot.soilSIM_fusion(): no raster panels found on this object.")
    return(invisible(NULL))
  }
  op <- graphics::par(mfrow = c(1, length(panels)))
  on.exit(graphics::par(op))
  for (nm in names(panels)) terra::plot(panels[[nm]], main = nm)
  invisible(x)
}
