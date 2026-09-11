# Generated from utils.R (P1b split): configuration & parameter management


#' Get Default Configuration
#'
#' Returns comprehensive default configuration for the soil workflow.
#'
#' @param config_type Configuration type ("full", "minimal", "validation")
#' @return Default configuration
#' @family utilities
#' @export
default_config <- function(config_type = "full") {

  base_config <- list(
    workflow = list(
      version = "1.0.0",
      description = "Soil simulation workflow configuration"
    ),

    data_processing = list(
      max_depth = DEFAULT_MAX_DEPTH_CM,
      exclude_unsuitable_horizons = TRUE,
      depth_units = "cm",
      missing_value_strategy = "interpolate"
    ),

    monte_carlo = list(
      n_simulations = 100,
      distribution_type = "triangular",
      confidence_level = 0.95,
      random_seed = 12345
    ),

    correlation = list(
      method = "pearson",
      minimum_samples = 30,
      handle_missing = "pairwise.complete.obs"
    ),

    gp_models = list(
      grouping_strategy = "auto",
      min_profiles_per_group = 3,
      min_observations_per_group = 15,
      optimize_hyperparameters = TRUE
    ),

    validation = list(
      run_validation = TRUE,
      quality_threshold = 0.7,
      generate_plots = TRUE,
      max_outlier_percentage = 5
    ),

    logging = list(
      log_level = "INFO",
      log_file = NULL,
      include_timestamp = TRUE
    )
  )

  if (config_type == "minimal") {
    # Return only essential configuration
    return(list(
      data_processing = base_config$data_processing,
      monte_carlo = list(n_simulations = base_config$monte_carlo$n_simulations),
      validation = list(run_validation = base_config$validation$run_validation)
    ))
  } else if (config_type == "validation") {
    # Return validation-focused configuration
    return(base_config$validation)
  } else {
    # Return full configuration
    return(base_config)
  }
}

#' Validate Parameters
#'
#' Validates input parameters for workflow functions.
#'
#' @param parameters Parameter list to validate
#' @param parameter_specs Parameter specifications
#' @param strict_mode Whether to use strict validation
#' @return Validation results
#' @keywords internal
validate_parameters <- function(parameters, parameter_specs, strict_mode = TRUE) {

  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    corrected_parameters = parameters
  )

  # Validate each parameter
  for (param_name in names(parameter_specs)) {
    param_spec <- parameter_specs[[param_name]]
    param_value <- parameters[[param_name]]

    # Check if required parameter is missing
    if (is.null(param_value) || length(param_value) == 0) {
      if (param_spec$required %||% FALSE) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Required parameter missing:", param_name))
        validation_result$valid <- FALSE
      }
      next
    }

    # Type validation
    if (!is.null(param_spec$type)) {
      expected_type <- param_spec$type

      type_valid <- switch(expected_type,
                           "character" = is.character(param_value),
                           "numeric" = is.numeric(param_value),
                           "integer" = is.integer(param_value) || (is.numeric(param_value) && all(param_value == as.integer(param_value))),
                           "logical" = is.logical(param_value),
                           "list" = is.list(param_value),
                           TRUE  # Default to valid for unknown types
      )

      if (!type_valid) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Type mismatch for", param_name, "- expected", expected_type, "got", class(param_value)[1]))
        validation_result$valid <- FALSE
        next
      }
    }

    # Range validation (for numeric parameters)
    if (!is.null(param_spec$range) && is.numeric(param_value)) {
      range_limits <- param_spec$range

      # Check each value in the vector
      out_of_range <- param_value < range_limits[1] | param_value > range_limits[2]

      if (any(out_of_range, na.rm = TRUE)) {
        if (strict_mode) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Value out of range for", param_name, "- expected",
                                              paste(range_limits, collapse = " to ")))
          validation_result$valid <- FALSE
        } else {
          validation_result$warnings <- c(validation_result$warnings,
                                          paste("Value out of range for", param_name))
        }
      }
    }

    # Choices validation (FIXED to handle vectors properly)
    if (!is.null(param_spec$choices)) {
      valid_choices <- param_spec$choices

      # Handle vector parameters properly
      if (length(param_value) == 1) {
        # Single value - check directly
        if (!param_value %in% valid_choices) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Invalid choice for", param_name, "- got", param_value,
                                              "expected one of:", paste(valid_choices, collapse = ", ")))
          validation_result$valid <- FALSE
        }
      } else {
        # Vector values - check each element
        invalid_choices <- param_value[!param_value %in% valid_choices]

        if (length(invalid_choices) > 0) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Invalid choices for", param_name, "- got",
                                              paste(invalid_choices, collapse = ", "),
                                              "expected from:", paste(valid_choices, collapse = ", ")))
          validation_result$valid <- FALSE
        }
      }
    }

    # Length validation
    if (!is.null(param_spec$min_length)) {
      if (length(param_value) < param_spec$min_length) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Parameter", param_name, "too short - got", length(param_value),
                                            "expected at least", param_spec$min_length))
        validation_result$valid <- FALSE
      }
    }

    if (!is.null(param_spec$max_length)) {
      if (length(param_value) > param_spec$max_length) {
        if (strict_mode) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Parameter", param_name, "too long - got", length(param_value),
                                              "expected at most", param_spec$max_length))
          validation_result$valid <- FALSE
        } else {
          validation_result$warnings <- c(validation_result$warnings,
                                          paste("Parameter", param_name, "longer than recommended"))
        }
      }
    }
  }

  return(validation_result)
}

#' Merge Configurations
#'
#' Merges user configuration with default configuration.
#'
#' @param default_config Default configuration
#' @param user_config User-provided configuration
#' @param deep_merge Whether to perform deep merging
#' @return Merged configuration
#' @keywords internal
merge_configurations <- function(default_config, user_config, deep_merge = TRUE) {

  if (is.null(user_config)) {
    return(default_config)
  }

  if (!deep_merge) {
    # Simple merge - user config overrides defaults
    merged_config <- default_config
    for (name in names(user_config)) {
      merged_config[[name]] <- user_config[[name]]
    }
    return(merged_config)
  }

  # Deep merge - recursively merge nested lists
  merge_recursive <- function(default, user) {
    if (!is.list(default) || !is.list(user)) {
      return(user)  # User value overrides default
    }

    result <- default
    for (name in names(user)) {
      if (name %in% names(default) && is.list(default[[name]]) && is.list(user[[name]])) {
        result[[name]] <- merge_recursive(default[[name]], user[[name]])
      } else {
        result[[name]] <- user[[name]]
      }
    }
    return(result)
  }

  return(merge_recursive(default_config, user_config))
}

#' Export Workflow Metadata
#'
#' Exports comprehensive metadata about the workflow execution.
#'
#' @param workflow_results Workflow results
#' @param output_path Output file path
#' @param include_data_summary Whether to include data summaries
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Success status
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need it, not dead/broken code.
#' @keywords internal
export_workflow_metadata <- function(workflow_results,
                                     output_path,
                                     include_data_summary = TRUE,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", paste("Exporting workflow metadata to:", output_path), category = "Metadata")

  metadata <- list(
    workflow_info = list(
      timestamp = Sys.time(),
      r_version = R.version.string,
      platform = R.version$platform,
      packages = get_package_versions()
    ),

    execution_summary = extract_execution_summary(workflow_results),

    configuration = extract_configuration_info(workflow_results),

    quality_assessment = extract_quality_summary(workflow_results)
  )

  if (include_data_summary) {
    metadata$data_summary <- extract_data_summary(workflow_results)
  }

  # Export as JSON
  tryCatch({
    jsonlite::write_json(metadata, output_path, pretty = TRUE, auto_unbox = TRUE)
    log_message("INFO", "Metadata exported successfully", category = "Metadata")
    return(TRUE)

  }, error = function(e) {
    log_message("ERROR", paste("Error exporting metadata:", e$message), category = "Metadata")
    return(FALSE)
  })
}

