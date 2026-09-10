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
  status <- x$overall_status %||% x$status %||%
    (if (!is.null(x$overall_quality_score)) .fmt_num(x$overall_quality_score) else NULL)
  if (!is.null(status)) cat("  overall:", status, "\n")
  cat("  sections:", paste(names(x), collapse = ", "), "\n")
  invisible(x)
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
