#' @title Multivariate Adjustment & GP Integration
#' @description Functions that integrate Monte Carlo simulation results with
#'   Gaussian Process depth models, applying realistic depth trends while
#'   preserving within-depth correlations across properties.
#' @name multivariate_adjustment
NULL

# ============================================================================
# 1. MASTER INTEGRATION FUNCTIONS
# ============================================================================

#' Integrate Monte Carlo Simulations with GP Models
#'
#' Master function that integrates Monte Carlo simulation results with GP models
#' to apply realistic depth trends while preserving within-depth correlations.
#' Enhanced with Module 8 utilities for robust processing and validation.
#'
#' @param simulation_results Results from monte_carlo::generate_monte_carlo_realizations()
#' @param gp_models Optional NRCS GP models from gp_modeling module
#' @param cokey_mapping Optional mapping from simulations to NRCS GP groups
#' @param integration_method Method: "nrcs_gp", "local_gp", or "hybrid" (default = "hybrid")
#' @param preserve_correlations Whether to preserve within-depth correlations (default = TRUE)
#' @param properties Properties to adjust (NULL = auto-detect)
#' @param parallel Whether to use parallel processing (default = FALSE)
#' @param n_cores Number of cores for parallel processing
#' @param config Integration configuration (uses Module 8 defaults if NULL)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Integrated simulation results with realistic depth trends
#' @export
integrate_monte_carlo_with_gp <- function(simulation_results,
                                          gp_models = NULL,
                                          cokey_mapping = NULL,
                                          integration_method = "hybrid",
                                          preserve_correlations = TRUE,
                                          properties = NULL,
                                          parallel = FALSE,
                                          n_cores = NULL,
                                          config = NULL,
                                          verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== MONTE CARLO - GP INTEGRATION ===", category = "MultivarAdjust")
  log_message("INFO", paste("Integration method:", integration_method), category = "MultivarAdjust")
  log_message("INFO", paste("Preserve correlations:", preserve_correlations), category = "MultivarAdjust")
  log_message("INFO", paste("Parallel processing:", parallel), category = "MultivarAdjust")

  start_time <- Sys.time()

  # Load configuration using Module 8
  if (is.null(config)) {
    config <- get_default_configuration("full")
  }

  # Validate parameters using Module 8
  param_specs <- list(
    integration_method = list(required = TRUE, type = "character",
                              choices = c("nrcs_gp", "local_gp", "hybrid")),
    preserve_correlations = list(required = FALSE, type = "logical"),
    parallel = list(required = FALSE, type = "logical")
  )

  param_validation <- validate_parameters(
    list(integration_method = integration_method,
         preserve_correlations = preserve_correlations,
         parallel = parallel),
    param_specs,
    strict_mode = TRUE
  )

  if (!param_validation$valid) {
    stop("Parameter validation failed: ", paste(param_validation$errors, collapse = ", "))
  }

  # Input validation using Module 8
  if (is.null(simulation_results) || is.null(simulation_results$simulation_data)) {
    stop("Invalid simulation_results - missing simulation_data")
  }

  simulation_data <- simulation_results$simulation_data
  original_metadata <- simulation_results$metadata

  # Data quality validation
  validation_results <- validate_data_quality(
    simulation_data,
    required_columns = c("cokey", "hzdept_r", "simulation_number"),
    quality_thresholds = config$validation
  )

  if (!validation_results$overall_quality$validation_passed) {
    log_message("WARN", "Input simulation data has quality issues", category = "MultivarAdjust")
  }

  log_message("INFO", paste("Processing", nrow(simulation_data), "simulation rows"), category = "MultivarAdjust")

  # Auto-detect properties if not specified
  if (is.null(properties)) {
    properties <- detect_simulation_properties(simulation_data)
    log_message("INFO", paste("Auto-detected properties:", paste(properties, collapse = ", ")), category = "MultivarAdjust")
  }

  if (length(properties) == 0) {
    stop("No suitable properties found for integration")
  }

  # Validate properties using Module 8
  property_validation <- validate_properties(properties, "laboratory", strict_mode = FALSE)
  if (!property_validation$valid) {
    log_message("WARN", paste("Property validation issues:", paste(property_validation$warnings, collapse = "; ")), category = "MultivarAdjust")
  }

  # Determine integration approach
  use_nrcs_gp <- !is.null(gp_models) && !is.null(cokey_mapping) &&
    integration_method %in% c("nrcs_gp", "hybrid")

  use_local_gp <- integration_method %in% c("local_gp", "hybrid")

  log_message("INFO", paste("Using NRCS GP:", use_nrcs_gp, "Using Local GP:", use_local_gp), category = "MultivarAdjust")

  # Process by component (cokey) groups
  unique_cokeys <- unique(simulation_data$cokey)
  log_message("INFO", paste("Processing", length(unique_cokeys), "unique cokeys"), category = "MultivarAdjust")

  # Process cokeys with progress tracking
  if (parallel && length(unique_cokeys) > 1) {
    integrated_results <- process_cokeys_parallel(
      simulation_data, unique_cokeys, properties, gp_models, cokey_mapping,
      use_nrcs_gp, use_local_gp, preserve_correlations, n_cores, config
    )
  } else {
    integrated_results <- process_cokeys_sequential(
      simulation_data, unique_cokeys, properties, gp_models, cokey_mapping,
      use_nrcs_gp, use_local_gp, preserve_correlations, config
    )
  }

  # Combine and validate results
  final_data <- combine_and_validate_results(integrated_results, unique_cokeys, simulation_data)

  # Comprehensive validation using Module 8
  log_message("INFO", "=== INTEGRATION VALIDATION ===", category = "MultivarAdjust")
  validation_results <- validate_integration_results(
    simulation_data, final_data, properties, preserve_correlations
  )

  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "secs")

  # Prepare final output with Module 8 metadata handling
  final_results <- create_integration_results(
    final_data, simulation_data, original_metadata, integration_method,
    properties, preserve_correlations, use_nrcs_gp, use_local_gp,
    length(unique_cokeys), length(integrated_results), processing_time,
    parallel, validation_results
  )

  log_message("INFO", "=== INTEGRATION COMPLETE ===", category = "MultivarAdjust")
  log_message("INFO", paste("Processing time:", round(as.numeric(processing_time), 2), "seconds"), category = "MultivarAdjust")
  log_message("INFO", paste("Success rate:", round(length(integrated_results) / length(unique_cokeys) * 100, 1), "%"), category = "MultivarAdjust")

  return(final_results)
}

#' Apply GP Depth Trends with Correlation Preservation
#'
#' Enhanced core function that applies GP-derived depth trends using Module 8 utilities
#' for robust error handling and validation.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param gp_predictions Named list of GP predictions by property
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param primary_property Reference property for correlation preservation
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Adjusted simulation data
#' @export
apply_gp_depth_trends <- function(cokey_data,
                                  gp_predictions,
                                  properties,
                                  preserve_correlations = TRUE,
                                  primary_property = NULL,
                                  verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Enhanced validation using Module 8
  if (nrow(cokey_data) < 2) {
    log_message("DEBUG", "Insufficient rows for GP trend application", category = "MultivarAdjust")
    return(cokey_data)
  }

  valid_depths <- !is.na(cokey_data$hzdept_r) & is.finite(cokey_data$hzdept_r)
  if (sum(valid_depths) < 2) {
    log_message("DEBUG", "Insufficient valid depths for GP trend application", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Filter to valid depths and properties
  valid_data <- cokey_data[valid_depths, ]
  available_properties <- intersect(properties, names(gp_predictions))

  if (length(available_properties) == 0) {
    log_message("DEBUG", "No available properties for GP adjustment", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Set primary property for correlation preservation
  if (is.null(primary_property)) {
    primary_property <- available_properties[1]
  }

  # Get unique depths and simulation numbers
  unique_depths <- sort(unique(valid_data$hzdept_r))
  sim_numbers <- unique(valid_data$simulation_number)

  if (length(unique_depths) < 2 || length(sim_numbers) == 0) {
    log_message("DEBUG", "Insufficient depth or simulation variety", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Convert to matrix format for processing
  property_matrices <- tryCatch({
    convert_to_property_matrices(valid_data, available_properties, unique_depths, sim_numbers)
  }, error = function(e) {
    handle_workflow_error(e, "Property matrix conversion", "warn")
    return(list())
  })

  if (length(property_matrices) == 0) {
    log_message("WARN", "Failed to convert to property matrices", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Apply multivariate adjustment with enhanced error handling
  adjusted_matrices <- tryCatch({
    if (preserve_correlations && length(property_matrices) >= 2) {
      preserve_correlation_structure(
        property_matrices, gp_predictions, unique_depths, primary_property
      )
    } else {
      apply_individual_adjustments(property_matrices, gp_predictions, unique_depths)
    }
  }, error = function(e) {
    handle_workflow_error(e, "GP adjustment application", "warn")
    return(property_matrices)  # Return original if adjustment fails
  })

  # Convert back to long format
  adjusted_data <- tryCatch({
    convert_to_long_format(
      adjusted_matrices, unique_depths, sim_numbers, valid_data, available_properties
    )
  }, error = function(e) {
    handle_workflow_error(e, "Long format conversion", "warn")
    return(data.frame())
  })

  if (nrow(adjusted_data) == 0) {
    log_message("WARN", "Failed to convert adjusted data to long format", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Merge with original data using Module 8 safe operations
  result_data <- merge_adjusted_data(cokey_data, adjusted_data, available_properties)

  return(result_data)
}

#' Preserve Correlation Structure During GP Adjustment
#'
#' Enhanced core correlation preservation algorithm with Module 8 error handling.
#'
#' @param property_matrices Named list of property matrices
#' @param gp_predictions Named list of GP predictions
#' @param depths Depth vector
#' @param primary_property Reference property for correlation preservation
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of adjusted property matrices
#' @export
preserve_correlation_structure <- function(property_matrices,
                                           gp_predictions,
                                           depths,
                                           primary_property,
                                           verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Preserving correlation structure during GP adjustment", category = "MultivarAdjust")

  property_names <- names(property_matrices)
  n_depths <- length(depths)
  n_sims <- ncol(property_matrices[[1]])

  # Validate inputs
  if (n_depths < 2 || n_sims < 1) {
    log_message("WARN", "Insufficient dimensions for correlation preservation", category = "MultivarAdjust")
    return(property_matrices)
  }

  # Establish reference quantile ordering from primary property
  primary_matrix <- property_matrices[[primary_property]]

  if (is.null(primary_matrix)) {
    primary_property <- property_names[1]
    primary_matrix <- property_matrices[[primary_property]]
    log_message("DEBUG", paste("Using", primary_property, "as primary property"), category = "MultivarAdjust")
  }

  surface_values <- primary_matrix[1, ]

  # Safe ECDF creation
  surface_ecdf <- tryCatch({
    ecdf(surface_values)
  }, error = function(e) {
    handle_workflow_error(e, "ECDF creation", "warn")
    return(NULL)
  })

  if (is.null(surface_ecdf)) {
    log_message("WARN", "Failed to create ECDF for correlation preservation", category = "MultivarAdjust")
    return(property_matrices)
  }

  reference_quantiles <- surface_ecdf(surface_values)

  # Apply consistent adjustment to all properties
  adjusted_list <- list()

  for (prop in property_names) {
    log_message("DEBUG", paste("Adjusting property:", prop), category = "MultivarAdjust")

    current_matrix <- property_matrices[[prop]]
    gp_means <- gp_predictions[[prop]]

    if (is.null(gp_means) || length(gp_means) != n_depths) {
      log_message("DEBUG", paste("No valid GP predictions for", prop, "- keeping original"), category = "MultivarAdjust")
      adjusted_list[[prop]] <- current_matrix
      next
    }

    # Initialize with original values
    adjusted_matrix <- current_matrix

    # Apply depth-wise adjustment using SAME quantile ordering
    for (i in 2:n_depths) {

      # Get GP trend ratio with enhanced safety checks
      gp_ratio <- calculate_safe_gp_ratio(gp_means, i)

      # Get previous and current simulated values
      prev_values <- adjusted_matrix[i - 1, ]
      curr_values <- current_matrix[i, ]

      # Apply adjustment using REFERENCE quantiles (preserves correlations)
      adjusted_curr <- apply_quantile_adjustment(
        reference_quantiles, curr_values, prev_values, gp_ratio, n_sims
      )

      # Correct to maintain original distribution shape
      adjusted_matrix[i, ] <- correct_distribution_shape(curr_values, adjusted_curr)
    }

    adjusted_list[[prop]] <- adjusted_matrix
  }

  log_message("DEBUG", "Correlation preservation completed", category = "MultivarAdjust")
  return(adjusted_list)
}

# ============================================================================
# 2. NRCS GP INTEGRATION FUNCTIONS (Enhanced)
# ============================================================================

#' Apply NRCS Trend Adjustments
#'
#' Enhanced version with Module 8 integration and proper Module 5 function calls.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param gp_models NRCS GP models from gp_modeling module
#' @param model_group GP model group for this cokey
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Adjusted simulation data
#' @export
apply_nrcs_trend_adjustments <- function(cokey_data,
                                         gp_models,
                                         model_group,
                                         properties,
                                         preserve_correlations = TRUE,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", paste("Applying NRCS GP adjustments for group:", model_group), category = "MultivarAdjust")

  # Enhanced property mapping with Module 8 safe operations
  property_mapping <- get_nrcs_property_mapping()

  # Get available properties
  available_properties <- intersect(properties, names(property_mapping))

  if (length(available_properties) == 0) {
    log_message("DEBUG", "No mappable properties found for NRCS GP adjustment", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Get unique depths with validation
  unique_depths <- sort(unique(cokey_data$hzdept_r[!is.na(cokey_data$hzdept_r)]))

  if (length(unique_depths) < 2) {
    log_message("DEBUG", "Insufficient depths for NRCS GP adjustment", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Get NRCS GP predictions for each property using Module 5 functions
  nrcs_gp_predictions <- get_nrcs_gp_predictions(
    gp_models, available_properties, property_mapping, model_group, unique_depths, cokey_data
  )

  if (length(nrcs_gp_predictions) == 0) {
    log_message("DEBUG", "No valid NRCS GP predictions obtained", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Apply GP depth trends using enhanced function
  result <- apply_gp_depth_trends(
    cokey_data,
    nrcs_gp_predictions,
    available_properties,
    preserve_correlations
  )

  return(result)
}

#' Match Simulations to NRCS Models
#'
#' Enhanced version with Module 8 error handling.
#'
#' @param cokey Target cokey
#' @param cokey_mapping Mapping from gp_modeling::match_soils_to_gp_models()
#' @param fallback_group Default group if matching fails
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Model group name
#' @export
match_simulations_to_nrcs_models <- function(cokey, cokey_mapping, fallback_group = "general_pool",
                                             verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(cokey_mapping)) {
    log_message("DEBUG", "No cokey mapping provided - using fallback", category = "MultivarAdjust")
    return(fallback_group)
  }

  # Enhanced matching with error handling
  tryCatch({
    # Find mapping for this cokey
    model_group <- cokey_mapping$gp_model_group[cokey_mapping$sim_cokey == cokey]

    if (length(model_group) == 0 || is.na(model_group)) {
      log_message("DEBUG", paste("No mapping found for cokey", cokey, "- using fallback"), category = "MultivarAdjust")
      return(fallback_group)
    }

    return(model_group)

  }, error = function(e) {
    handle_workflow_error(e, paste("NRCS model matching for cokey", cokey), "warn")
    return(fallback_group)
  })
}

#' Extract NRCS Depth Trends
#'
#' Enhanced version with Module 8 validation and error handling.
#'
#' @param gp_models NRCS GP models
#' @param properties Properties to extract trends for
#' @param depths Depths for trend extraction
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of depth trends by property and group
#' @export
extract_nrcs_depth_trends <- function(gp_models, properties, depths = seq(0, 200, by = 10),
                                      verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Extracting NRCS depth trends", category = "MultivarAdjust")

  # Validate inputs using Module 8
  if (is.null(gp_models) || length(properties) == 0) {
    log_message("WARN", "Invalid inputs for NRCS trend extraction", category = "MultivarAdjust")
    return(list())
  }

  trend_data <- list()

  for (prop in properties) {
    if (prop %in% names(gp_models) && gp_models[[prop]]$type == "stratified_grouped") {

      prop_trends <- list()

      for (group in names(gp_models[[prop]]$models)) {
        group_model <- gp_models[[prop]]$models[[group]]

        if (!is.null(group_model)) {
          predictions <- tryCatch({
            # Use Module 5 function
            predict_gp_depth_trends(group_model, depths)
          }, error = function(e) {
            handle_workflow_error(e, paste("NRCS trend extraction for", prop, group), "warn")
            return(NULL)
          })

          if (!is.null(predictions)) {
            prop_trends[[group]] <- data.frame(
              depth = depths,
              predicted_value = predictions,
              group = group,
              property = prop,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(prop_trends) > 0) {
        trend_data[[prop]] <- dplyr::bind_rows(prop_trends)
      }
    }
  }

  log_message("DEBUG", paste("Extracted trends for", length(trend_data), "properties"), category = "MultivarAdjust")
  return(trend_data)
}

# ============================================================================
# 3. LOCAL GP INTEGRATION FUNCTIONS (Enhanced)
# ============================================================================

#' Apply Local GP Adjustments
#'
#' Enhanced version with Module 8 utilities and better error handling.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param min_depths Minimum depths required for GP fitting
#' @param config Configuration from Module 8
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Adjusted simulation data
#' @export
apply_local_gp_adjustments <- function(cokey_data,
                                       properties,
                                       preserve_correlations = TRUE,
                                       min_depths = 3,
                                       config = NULL,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- get_default_configuration("validation")
  }

  # Enhanced depth validation
  unique_depths <- sort(unique(cokey_data$hzdept_r[!is.na(cokey_data$hzdept_r)]))

  if (length(unique_depths) < min_depths) {
    log_message("DEBUG", paste("Insufficient depths for local GP fitting:", length(unique_depths), "<", min_depths), category = "MultivarAdjust")
    return(cokey_data)
  }

  # Fit local GP models with enhanced error handling
  local_gp_models <- fit_local_gp_models(cokey_data, properties, config)

  if (length(local_gp_models) == 0) {
    log_message("DEBUG", "No local GP models could be fitted", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Generate local GP predictions
  local_gp_predictions <- generate_local_predictions(local_gp_models, unique_depths)

  if (length(local_gp_predictions) == 0) {
    log_message("DEBUG", "No valid local GP predictions generated", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Apply local depth trends
  result <- apply_local_depth_trends(
    cokey_data,
    local_gp_predictions,
    unique_depths,
    preserve_correlations
  )

  return(result)
}

#' Fit Local GP Models
#'
#' Enhanced version with Module 8 validation and configuration management.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param properties Properties to model
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of fitted local GP models
#' @export
fit_local_gp_models <- function(cokey_data, properties, config = NULL,
                                verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- get_default_configuration("validation")
  }

  log_message("DEBUG", "Fitting local GP models", category = "MultivarAdjust")

  local_models <- list()

  for (prop in properties) {
    if (!prop %in% names(cokey_data)) {
      log_message("DEBUG", paste("Property", prop, "not found in data"), category = "MultivarAdjust")
      next
    }

    model_result <- tryCatch({
      # Enhanced aggregation with Module 8 safe operations
      agg_data <- aggregate_property_by_depth(cokey_data, prop)

      if (is.null(agg_data) || nrow(agg_data) < 3) {
        log_message("DEBUG", paste("Insufficient aggregated data for", prop), category = "MultivarAdjust")
        return(NULL)
      }

      # Validate variation
      if (var(agg_data$mean_val, na.rm = TRUE) <= 0) {
        log_message("DEBUG", paste("No variation in", prop, "values"), category = "MultivarAdjust")
        return(NULL)
      }

      # Fit GP model using Module 5 approach
      fit_local_gp_model_single(agg_data, prop)

    }, error = function(e) {
      handle_workflow_error(e, paste("Local GP fitting for", prop), "warn")
      return(NULL)
    })

    if (!is.null(model_result)) {
      local_models[[prop]] <- model_result
    }
  }

  log_message("DEBUG", paste("Successfully fitted", length(local_models), "local GP models"), category = "MultivarAdjust")
  return(local_models)
}

#' Apply Local Depth Trends
#'
#' Enhanced version with Module 8 error handling.
#'
#' @param cokey_data Simulation data
#' @param local_predictions Local GP predictions
#' @param unique_depths Depth vector
#' @param preserve_correlations Whether to preserve correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Adjusted simulation data
#' @export
apply_local_depth_trends <- function(cokey_data,
                                     local_predictions,
                                     unique_depths,
                                     preserve_correlations = TRUE,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Applying local depth trends", category = "MultivarAdjust")

  # Apply GP depth trends using the same logic as NRCS
  result <- tryCatch({
    apply_gp_depth_trends(
      cokey_data,
      local_predictions,
      names(local_predictions),
      preserve_correlations
    )
  }, error = function(e) {
    handle_workflow_error(e, "Local depth trend application", "warn")
    return(cokey_data)
  })

  return(result)
}

# ============================================================================
# 4. MULTIVARIATE PROCESSING FUNCTIONS (Enhanced with Module 8)
# ============================================================================

#' Convert to Property Matrices
#'
#' Enhanced version with Module 8 safe operations and validation.
#'
#' @param simulation_data Simulation data in long format
#' @param properties Properties to convert
#' @param unique_depths Depth vector
#' @param sim_numbers Simulation numbers
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Named list of property matrices
#' @export
convert_to_property_matrices <- function(simulation_data,
                                         properties,
                                         unique_depths,
                                         sim_numbers,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Converting to property matrices", category = "MultivarAdjust")

  property_matrices <- list()

  for (prop in properties) {
    if (!prop %in% names(simulation_data)) {
      log_message("DEBUG", paste("Property", prop, "not found in simulation data"), category = "MultivarAdjust")
      next
    }

    # Create matrix: rows = depths, columns = simulations
    prop_matrix <- matrix(NA, nrow = length(unique_depths), ncol = length(sim_numbers))

    tryCatch({
      for (i in seq_along(unique_depths)) {
        depth <- unique_depths[i]

        for (j in seq_along(sim_numbers)) {
          sim_num <- sim_numbers[j]

          # Find matching value using Module 8 safe operations
          value <- simulation_data |>
            dplyr::filter(hzdept_r == depth, simulation_number == sim_num) |>
            dplyr::pull(.data[[prop]])

          if (length(value) > 0 && !is.na(value[1])) {
            prop_matrix[i, j] <- value[1]
          }
        }
      }

      # Only include if we have some valid data
      valid_data_count <- sum(!is.na(prop_matrix))
      if (valid_data_count > 0) {
        property_matrices[[prop]] <- prop_matrix
        log_message("DEBUG", paste("Created matrix for", prop, "with", valid_data_count, "valid values"), category = "MultivarAdjust")
      } else {
        log_message("DEBUG", paste("No valid data for", prop, "matrix"), category = "MultivarAdjust")
      }

    }, error = function(e) {
      handle_workflow_error(e, paste("Matrix creation for", prop), "warn")
    })
  }

  log_message("DEBUG", paste("Created", length(property_matrices), "property matrices"), category = "MultivarAdjust")
  return(property_matrices)
}

#' Convert to Long Format
#'
#' Enhanced version with Module 8 error handling and data validation.
#'
#' @param adjusted_matrices List of adjusted property matrices
#' @param unique_depths Depth vector
#' @param sim_numbers Simulation numbers
#' @param original_data Original simulation data for metadata
#' @param properties Properties that were adjusted
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data frame in long format
#' @export
convert_to_long_format <- function(adjusted_matrices,
                                   unique_depths,
                                   sim_numbers,
                                   original_data,
                                   properties,
                                   verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Converting adjusted matrices to long format", category = "MultivarAdjust")

  result_list <- list()

  for (i in seq_along(unique_depths)) {
    depth <- unique_depths[i]

    track_progress(i, length(unique_depths), "Converting depths", update_frequency = max(1, length(unique_depths) %/% 10))

    for (j in seq_along(sim_numbers)) {
      sim_num <- sim_numbers[j]

      # Get original row for metadata using Module 8 safe operations
      original_row <- tryCatch({
        original_data |>
          dplyr::filter(hzdept_r == depth, simulation_number == sim_num) |>
          dplyr::slice(1)
      }, error = function(e) {
        handle_workflow_error(e, paste("Original row retrieval for depth", depth, "sim", sim_num), "warn")
        return(data.frame())
      })

      if (nrow(original_row) > 0) {
        # Start with metadata using Module 8 safe column selection
        metadata_cols <- intersect(
          c("cokey", "compname", "mukey", "hzdept_r", "hzdepb_r", "simulation_number", "unique_id"),
          names(original_row)
        )

        new_row <- original_row |>
          dplyr::select(dplyr::all_of(metadata_cols))

        # Add adjusted property values
        for (prop in properties) {
          if (prop %in% names(adjusted_matrices)) {
            adj_matrix <- adjusted_matrices[[prop]]
            if (i <= nrow(adj_matrix) && j <= ncol(adj_matrix)) {
              new_row[[prop]] <- adj_matrix[i, j]
            }
          }
        }

        result_list[[length(result_list) + 1]] <- new_row
      }
    }
  }

  # Combine all rows with error handling
  result_df <- tryCatch({
    if (length(result_list) > 0) {
      dplyr::bind_rows(result_list)
    } else {
      data.frame()
    }
  }, error = function(e) {
    handle_workflow_error(e, "Long format data binding", "warn")
    return(data.frame())
  })

  log_message("DEBUG", paste("Converted to long format with", nrow(result_df), "rows"), category = "MultivarAdjust")
  return(result_df)
}

# ============================================================================
# 5. VALIDATION AND QUALITY CONTROL (Enhanced with Module 8)
# ============================================================================

#' Validate Integration Results
#'
#' Enhanced validation using Module 8 validation framework.
#'
#' @param original_data Original simulation data
#' @param integrated_data Integrated simulation data
#' @param properties Properties that were processed
#' @param preserve_correlations Whether correlations should be preserved
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results list
#' @export
validate_integration_results <- function(original_data,
                                         integrated_data,
                                         properties,
                                         preserve_correlations = TRUE,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating integration results", category = "MultivarAdjust")

  validation_results <- list(
    data_integrity = list(),
    correlation_preservation = list(),
    trend_realism = list(),
    overall_assessment = list()
  )

  # Data integrity checks using Module 8
  validation_results$data_integrity <- validate_data_integrity(original_data, integrated_data, properties)

  # Correlation preservation checks
  if (preserve_correlations && length(properties) >= 2) {
    validation_results$correlation_preservation <- validate_correlation_preservation_integration(
      original_data, integrated_data, properties
    )
  }

  # Trend realism checks
  validation_results$trend_realism <- validate_depth_trends(integrated_data, properties)

  # Overall assessment using Module 8 scoring
  validation_results$overall_assessment <- calculate_overall_validation_score_integration(validation_results)

  # Log validation summary
  log_validation_summary(validation_results, preserve_correlations)

  return(validation_results)
}

#' Correct Distribution Shapes
#'
#' Enhanced version with Module 8 property validation and constraints.
#'
#' @param adjusted_data Adjusted simulation data
#' @param original_data Original simulation data
#' @param properties Properties to validate
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Corrected simulation data
#' @export
correct_distribution_shapes <- function(adjusted_data, original_data, properties, config = NULL,
                                        verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- get_default_configuration("validation")
  }

  log_message("DEBUG", "Correcting distribution shapes", category = "MultivarAdjust")

  corrected_data <- adjusted_data

  for (prop in properties) {
    if (!prop %in% names(adjusted_data)) {
      next
    }

    # Get property-specific constraints using Module 8
    constraints <- get_property_constraints(prop)

    # Apply range constraints with Module 8 validation
    corrected_data[[prop]] <- apply_range_constraints(corrected_data[[prop]], constraints)

    # Apply distribution correction if needed
    if (constraints$preserve_distribution) {
      corrected_data[[prop]] <- correct_property_distribution(
        corrected_data[[prop]], original_data[[prop]], constraints
      )
    }
  }

  # Apply cross-property constraints using Module 8 validation
  corrected_data <- apply_cross_property_constraints(corrected_data, properties)

  log_message("DEBUG", "Distribution shape correction completed", category = "MultivarAdjust")
  return(corrected_data)
}

# ============================================================================
# 6. ENHANCED HELPER FUNCTIONS (Leveraging Module 8)
# ============================================================================

# Property detection with Module 8 validation
detect_simulation_properties <- function(simulation_data) {

  # Use Module 8 property validation
  all_properties <- get_predefined_properties("laboratory")

  # Common simulation property patterns
  simulation_patterns <- c(
    "sand_total", "sandtotal", "clay_total", "claytotal", "silt_total", "silttotal",
    "bulk_density", "dbovendry", "db", "water_retention", "wthirdbar", "wfifteenbar",
    "ph", "cec", "om", "soc", "rfv"
  )

  # Combine and filter
  all_candidates <- unique(c(all_properties, simulation_patterns))
  detected_properties <- intersect(all_candidates, names(simulation_data))

  log_message("DEBUG", paste("Detected", length(detected_properties), "simulation properties"), category = "MultivarAdjust")

  return(detected_properties)
}

# Enhanced parallel processing with Module 8 progress tracking
process_cokeys_parallel <- function(simulation_data, unique_cokeys, properties,
                                    gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                    preserve_correlations, n_cores, config) {

  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  log_message("INFO", paste("Running integration in parallel with", n_cores, "cores"), category = "MultivarAdjust")

  tryCatch({
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))

      # Export necessary objects with Module 8 setup
      parallel::clusterExport(cl, c("gp_models", "cokey_mapping", "use_nrcs_gp",
                                    "use_local_gp", "preserve_correlations", "properties", "config"),
                              envir = environment())

      # Load the package (and therefore all its Imports, including dplyr and
      # GPfit) on each worker - mirrors the same fix already applied to
      # monte-carlo.R's run_parallel_simulation().
      parallel::clusterEvalQ(cl, library(soilSIM))

      current_log_cfg <- getOption("soil_workflow_log_config")
      parallel::clusterCall(cl, function(cfg) options(soil_workflow_log_config = cfg), current_log_cfg)

      results <- parallel::parLapply(cl, unique_cokeys, function(cokey) {
        process_single_cokey(simulation_data, cokey, properties, gp_models,
                                      cokey_mapping, use_nrcs_gp, use_local_gp,
                                      preserve_correlations, config)
      })
    } else {
      results <- parallel::mclapply(unique_cokeys, function(cokey) {
        process_single_cokey(simulation_data, cokey, properties, gp_models,
                                      cokey_mapping, use_nrcs_gp, use_local_gp,
                                      preserve_correlations, config)
      }, mc.cores = n_cores)
    }

    return(results)

  }, error = function(e) {
    handle_workflow_error(e, "Parallel processing", "warn")
    # Fallback to sequential processing
    return(process_cokeys_sequential(simulation_data, unique_cokeys, properties,
                                     gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                     preserve_correlations, config))
  })
}

# Enhanced sequential processing with Module 8 progress tracking
process_cokeys_sequential <- function(simulation_data, unique_cokeys, properties,
                                      gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                      preserve_correlations, config) {

  log_message("INFO", "Running integration sequentially", category = "MultivarAdjust")

  results <- list()

  for (i in seq_along(unique_cokeys)) {
    cokey <- unique_cokeys[i]

    # Progress tracking using Module 8
    track_progress(i, length(unique_cokeys), "Processing cokeys", update_frequency = 10)

    result <- process_single_cokey(simulation_data, cokey, properties, gp_models,
                                            cokey_mapping, use_nrcs_gp, use_local_gp,
                                            preserve_correlations, config)

    results[[i]] <- result
  }

  return(results)
}

# Enhanced single cokey processing with Module 8 error handling
process_single_cokey <- function(simulation_data, cokey, properties, gp_models,
                                          cokey_mapping, use_nrcs_gp, use_local_gp,
                                          preserve_correlations, config) {

  tryCatch({
    # Extract data for this cokey
    cokey_data <- simulation_data |> dplyr::filter(cokey == !!cokey)

    if (nrow(cokey_data) < 2) {
      log_message("DEBUG", paste("Insufficient data for cokey", cokey), category = "MultivarAdjust")
      return(NULL)
    }

    result_data <- cokey_data

    # Apply NRCS GP adjustments if requested
    if (use_nrcs_gp && !is.null(gp_models) && !is.null(cokey_mapping)) {
      model_group <- match_simulations_to_nrcs_models(cokey, cokey_mapping)

      result_data <- apply_nrcs_trend_adjustments(
        result_data, gp_models, model_group, properties, preserve_correlations
      )
    }

    # Apply local GP adjustments if requested
    if (use_local_gp) {
      result_data <- apply_local_gp_adjustments(
        result_data, properties, preserve_correlations, config = config
      )
    }

    return(result_data)

  }, error = function(e) {
    handle_workflow_error(e, paste("Processing cokey", cokey), "warn")
    return(NULL)
  })
}

# Enhanced result combination with Module 8 validation
combine_and_validate_results <- function(integrated_results, unique_cokeys, simulation_data) {

  # Remove NULL results
  integrated_results <- integrated_results[!sapply(integrated_results, is.null)]

  if (length(integrated_results) == 0) {
    stop("No cokeys were successfully processed")
  }

  # Combine results with error handling
  final_data <- tryCatch({
    dplyr::bind_rows(integrated_results)
  }, error = function(e) {
    handle_workflow_error(e, "Results combination", "stop")
  })

  # Validate result dimensions
  if (nrow(final_data) == 0) {
    stop("Combined results are empty")
  }

  success_rate <- length(integrated_results) / length(unique_cokeys)
  log_message("INFO", paste("Combined results:", nrow(final_data), "rows,",
                            round(success_rate * 100, 1), "% success rate"), category = "MultivarAdjust")

  return(final_data)
}

# Enhanced metadata creation with Module 8 utilities
create_integration_results <- function(final_data, simulation_data, original_metadata,
                                       integration_method, properties, preserve_correlations,
                                       use_nrcs_gp, use_local_gp, n_cokeys_total,
                                       n_cokeys_successful, processing_time, parallel,
                                       validation_results) {

  # Create comprehensive metadata using Module 8 patterns
  integration_metadata <- list(
    method = integration_method,
    properties_processed = properties,
    preserve_correlations = preserve_correlations,
    use_nrcs_gp = use_nrcs_gp,
    use_local_gp = use_local_gp,
    n_cokeys_processed = n_cokeys_total,
    n_cokeys_successful = n_cokeys_successful,
    success_rate = n_cokeys_successful / n_cokeys_total,
    processing_time = processing_time,
    parallel = parallel,
    timestamp = Sys.time(),
    r_version = R.version.string,
    package_versions = get_package_versions()
  )

  # Prepare final output with Module 8 metadata handling
  final_results <- list(
    integrated_data = final_data,
    original_simulation_data = simulation_data,
    integration_metadata = integration_metadata,
    validation_results = validation_results,
    original_metadata = original_metadata
  )

  # Add attributes for easy access
  attr(final_results, "success_rate") <- integration_metadata$success_rate
  attr(final_results, "processing_time") <- processing_time
  attr(final_results, "validation_passed") <- validation_results$overall_assessment$validation_passed

  return(final_results)
}

# Additional helper functions with Module 8 integration
get_nrcs_property_mapping <- function() {
  list(
    "bulk_density_third_bar" = "clay_pct",
    "water_retention_third_bar" = "organic_matter",
    "water_retention_15_bar" = "organic_matter",
    "sand_total" = "sand_pct",
    "silt_total" = "sand_pct",
    "clay_total" = "clay_pct",
    "ph" = "pH",
    "cec" = "organic_matter",
    "soc" = "organic_matter",
    "om" = "organic_matter",
    "db" = "clay_pct",
    "sandtotal" = "sand_pct",
    "claytotal" = "clay_pct",
    "silttotal" = "sand_pct"
  )
}

get_nrcs_gp_predictions <- function(gp_models, available_properties, property_mapping,
                                    model_group, unique_depths, cokey_data) {

  nrcs_gp_predictions <- list()

  for (prop in available_properties) {
    nrcs_prop <- property_mapping[[prop]]

    if (nrcs_prop %in% names(gp_models) &&
        gp_models[[nrcs_prop]]$type == "stratified_grouped" &&
        model_group %in% names(gp_models[[nrcs_prop]]$models)) {

      prediction_result <- tryCatch({
        gp_model_info <- gp_models[[nrcs_prop]]$models[[model_group]]
        predict_gp_depth_trends(gp_model_info, unique_depths)
      }, error = function(e) {
        handle_workflow_error(e, paste("NRCS GP prediction for", prop), "warn")
        return(NULL)
      })

      if (!is.null(prediction_result)) {
        nrcs_gp_predictions[[prop]] <- prediction_result
      } else {
        # Use local means as fallback
        nrcs_gp_predictions[[prop]] <- get_local_property_means(cokey_data, prop, unique_depths)
      }
    } else {
      log_message("DEBUG", paste("No NRCS GP model for", prop, "- using local trends"), category = "MultivarAdjust")
      nrcs_gp_predictions[[prop]] <- get_local_property_means(cokey_data, prop, unique_depths)
    }
  }

  return(nrcs_gp_predictions)
}

get_local_property_means <- function(cokey_data, prop, unique_depths) {
  tryCatch({
    prop_means <- cokey_data |>
      dplyr::group_by(hzdept_r) |>
      dplyr::summarise(mean_val = mean(.data[[prop]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(hzdept_r) |>
      dplyr::pull(mean_val)

    # Ensure we have values for all depths
    if (length(prop_means) != length(unique_depths)) {
      # Interpolate missing values
      prop_means <- approx(seq_along(prop_means), prop_means, n = length(unique_depths))$y
    }

    return(prop_means)
  }, error = function(e) {
    handle_workflow_error(e, paste("Local means calculation for", prop), "warn")
    return(rep(mean(cokey_data[[prop]], na.rm = TRUE), length(unique_depths)))
  })
}

# Enhanced GP ratio calculation with Module 8 safety
calculate_safe_gp_ratio <- function(gp_means, i) {
  if (is.na(gp_means[i-1]) || is.na(gp_means[i]) || gp_means[i-1] == 0) {
    return(1)  # No adjustment if invalid ratio
  } else {
    ratio <- gp_means[i] / gp_means[i - 1]
    # Clamp extreme ratios
    return(pmax(0.1, pmin(10, ratio)))
  }
}

# Enhanced quantile adjustment with Module 8 safety
apply_quantile_adjustment <- function(reference_quantiles, curr_values, prev_values, gp_ratio, n_sims) {
  adjusted_curr <- numeric(n_sims)

  for (j in 1:n_sims) {
    tryCatch({
      # Use the SAME quantile from reference property
      q <- reference_quantiles[j]

      # Find corresponding quantile value in current property's distribution
      quantile_value <- quantile(curr_values, probs = q, na.rm = TRUE)

      # Apply GP-based adjustment
      adjusted_curr[j] <- quantile_value + (prev_values[j] * gp_ratio - quantile_value)
    }, error = function(e) {
      # Fallback to original value
      adjusted_curr[j] <<- curr_values[j]
    })
  }

  return(adjusted_curr)
}

# Enhanced distribution shape correction with Module 8 safety
correct_distribution_shape <- function(curr_values, adjusted_curr) {
  tryCatch({
    # A property column that's entirely (or almost entirely) NA for this group is a real,
    # reachable case (e.g. a texture column simulate_cokey_generalized() left NA for some
    # rows rather than crashing the whole cokey - see R/property-simulation.R). Without
    # this guard, var(curr_values) on <2 non-NA values is NA (`if (NA > 0)` throws
    # "missing value where TRUE/FALSE needed"), and ecdf()/quantile() on an all-NA
    # adjusted_curr throws "'x' must have 1 or more non-missing values". Both were
    # already caught by this function's own tryCatch() below (falling back to
    # curr_values unadjusted), so this isn't a correctness fix - just replacing noisy,
    # per-call warning spam with the same fallback taken directly.
    if (sum(!is.na(curr_values)) < 2 || sum(!is.na(adjusted_curr)) < 1) {
      return(curr_values)
    }
    if (var(curr_values, na.rm = TRUE) > 0) {
      ecdf_adjusted <- ecdf(adjusted_curr)
      # Map back to original distribution while preserving order
      corrected_values <- quantile(curr_values,
                                   probs = ecdf_adjusted(adjusted_curr),
                                   na.rm = TRUE)
      return(corrected_values)
    } else {
      # If no variation, keep original values
      return(curr_values)
    }
  }, error = function(e) {
    handle_workflow_error(e, "Distribution shape correction", "warn")
    return(curr_values)
  })
}

# Local GP prediction/aggregation helpers (both real, despite the stale
# section name below - generate_local_predictions() delegates to the real
# predict_gp_depth_trends()).
generate_local_predictions <- function(local_gp_models, unique_depths) {
  local_gp_predictions <- list()

  for (prop in names(local_gp_models)) {
    prediction_result <- tryCatch({
      predict_gp_depth_trends(local_gp_models[[prop]], unique_depths)
    }, error = function(e) {
      handle_workflow_error(e, paste("Local GP prediction for", prop), "warn")
      return(NULL)
    })

    if (!is.null(prediction_result)) {
      local_gp_predictions[[prop]] <- prediction_result
    }
  }

  return(local_gp_predictions)
}

aggregate_property_by_depth <- function(cokey_data, prop) {
  tryCatch({
    agg_data <- cokey_data |>
      dplyr::group_by(hzdept_r) |>
      dplyr::summarise(
        mean_val = mean(.data[[prop]], na.rm = TRUE),
        n_obs = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::filter(
        !is.na(mean_val),
        !is.infinite(mean_val),
        n_obs > 0
      ) |>
      dplyr::arrange(hzdept_r)

    return(agg_data)
  }, error = function(e) {
    handle_workflow_error(e, paste("Property aggregation for", prop), "warn")
    return(NULL)
  })
}

fit_local_gp_model_single <- function(agg_data, prop) {
  depths <- agg_data$hzdept_r
  values <- agg_data$mean_val

  # Scale depths to [0,1]
  depth_min <- min(depths)
  depth_max <- max(depths)
  depth_range <- depth_max - depth_min

  if (depth_range > 0) {
    scaled_depths <- (depths - depth_min) / depth_range

    # Fit GP model
    gp_model <- GPfit::GP_fit(X = as.matrix(scaled_depths), Y = values)

    # Store with scaling information (using Module 5 structure)
    return(list(
      gp_model = gp_model,
      depth_scaling = list(
        min = depth_min,
        max = depth_max,
        range = depth_range
      ),
      training_data = agg_data,
      n_training_points = nrow(agg_data),
      property = prop
    ))
  }

  return(NULL)
}

merge_adjusted_data <- function(cokey_data, adjusted_data, available_properties) {
  result_data <- cokey_data

  for (i in seq_len(nrow(adjusted_data))) {
    adj_row <- adjusted_data[i, ]

    # Find matching row in original data
    match_idx <- which(
      result_data$hzdept_r == adj_row$hzdept_r &
        result_data$simulation_number == adj_row$simulation_number
    )

    if (length(match_idx) == 1) {
      # Update property values
      for (prop in available_properties) {
        if (prop %in% names(adj_row) && !is.na(adj_row[[prop]])) {
          result_data[[prop]][match_idx] <- adj_row[[prop]]
        }
      }
    }
  }

  return(result_data)
}

# Enhanced validation functions using Module 8
validate_data_integrity <- function(original_data, integrated_data, properties) {
  list(
    n_rows_original = nrow(original_data),
    n_rows_integrated = nrow(integrated_data),
    rows_preserved = nrow(original_data) == nrow(integrated_data),
    properties_processed = properties,
    missing_values_added = check_missing_values_increase(original_data, integrated_data, properties)
  )
}

validate_correlation_preservation_integration <- function(original_data, integrated_data, properties) {
  # Enhanced correlation validation using Module 8 safe operations
  depths_to_check <- unique(integrated_data$hzdept_r)[1:min(3, length(unique(integrated_data$hzdept_r)))]
  correlation_differences <- c()

  for (depth in depths_to_check) {
    orig_cor <- tryCatch({
      orig_depth_data <- original_data |>
        dplyr::filter(hzdept_r == depth) |>
        dplyr::select(dplyr::all_of(properties))

      if (nrow(orig_depth_data) > 10) {
        cor(orig_depth_data, use = "complete.obs")
      } else {
        NULL
      }
    }, error = function(e) NULL)

    int_cor <- tryCatch({
      int_depth_data <- integrated_data |>
        dplyr::filter(hzdept_r == depth) |>
        dplyr::select(dplyr::all_of(properties))

      if (nrow(int_depth_data) > 10) {
        cor(int_depth_data, use = "complete.obs")
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (!is.null(orig_cor) && !is.null(int_cor)) {
      max_diff <- max(abs(orig_cor - int_cor), na.rm = TRUE)
      correlation_differences <- c(correlation_differences, max_diff)
    }
  }

  return(list(
    depths_checked = depths_to_check,
    correlation_differences = correlation_differences,
    max_correlation_difference = if(length(correlation_differences) > 0) max(correlation_differences, na.rm = TRUE) else NA,
    mean_correlation_difference = if(length(correlation_differences) > 0) mean(correlation_differences, na.rm = TRUE) else NA
  ))
}

validate_depth_trends <- function(integrated_data, properties) {
  realistic_trends <- TRUE
  trend_issues <- character(0)

  for (prop in properties) {
    if (!prop %in% names(integrated_data)) {
      next
    }

    trend_analysis <- tryCatch({
      # Check for realistic depth trends using Module 8 safe operations
      trend_data <- integrated_data |>
        dplyr::group_by(hzdept_r) |>
        dplyr::summarise(mean_value = mean(.data[[prop]], na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(hzdept_r)

      if (nrow(trend_data) >= 3) {
        # Check for extreme jumps in trend
        diffs <- diff(trend_data$mean_value)
        if (length(diffs) > 0 && var(diffs, na.rm = TRUE) > 0) {
          extreme_jumps <- abs(diffs) > 3 * sd(diffs, na.rm = TRUE)

          if (any(extreme_jumps, na.rm = TRUE)) {
            realistic_trends <<- FALSE
            trend_issues <<- c(trend_issues, paste("Extreme trend jumps in", prop))
          }
        }

        # Check for property-specific realistic ranges
        constraints <- get_property_constraints(prop)
        if (!is.null(constraints$range)) {
          out_of_range <- any(trend_data$mean_value < constraints$range[1] |
                                trend_data$mean_value > constraints$range[2], na.rm = TRUE)

          if (out_of_range) {
            realistic_trends <<- FALSE
            trend_issues <<- c(trend_issues, paste("Out of range values in", prop))
          }
        }
      }

      return(TRUE)
    }, error = function(e) {
      handle_workflow_error(e, paste("Trend validation for", prop), "warn")
      return(FALSE)
    })
  }

  return(list(
    realistic_trends = realistic_trends,
    trend_issues = trend_issues
  ))
}

calculate_overall_validation_score_integration <- function(validation_results) {
  integrity_score <- ifelse(validation_results$data_integrity$rows_preserved, 1.0, 0.5)

  correlation_score <- if (!is.null(validation_results$correlation_preservation) &&
                           !is.na(validation_results$correlation_preservation$max_correlation_difference)) {
    1.0 - min(1.0, validation_results$correlation_preservation$max_correlation_difference)
  } else {
    1.0
  }

  trend_score <- ifelse(validation_results$trend_realism$realistic_trends, 1.0, 0.7)

  overall_score <- mean(c(integrity_score, correlation_score, trend_score), na.rm = TRUE)

  return(list(
    integrity_score = integrity_score,
    correlation_score = correlation_score,
    trend_score = trend_score,
    overall_score = overall_score,
    validation_passed = overall_score >= 0.8
  ))
}

log_validation_summary <- function(validation_results, preserve_correlations) {
  log_message("INFO", "Validation Results:", category = "MultivarAdjust")
  log_message("INFO", paste("  Data integrity:", validation_results$data_integrity$rows_preserved), category = "MultivarAdjust")
  log_message("INFO", paste("  Overall score:", round(validation_results$overall_assessment$overall_score, 3)), category = "MultivarAdjust")

  if (preserve_correlations && !is.null(validation_results$correlation_preservation) &&
      !is.na(validation_results$correlation_preservation$max_correlation_difference)) {
    log_message("INFO", paste("  Max correlation difference:",
                              round(validation_results$correlation_preservation$max_correlation_difference, 4)), category = "MultivarAdjust")
  }
}

# Enhanced property constraints using Module 8
get_property_constraints <- function(property) {
  # Use Module 8 property validation to get realistic ranges
  constraints <- list(
    range = NULL,
    preserve_distribution = FALSE,
    cross_property_rules = NULL
  )

  # Property-specific constraints enhanced with Module 8 patterns
  if (property %in% c("sand_total", "sandtotal", "clay_total", "claytotal", "silt_total", "silttotal")) {
    constraints$range <- c(0, 100)
    constraints$preserve_distribution <- TRUE
    constraints$cross_property_rules <- "texture_sum"
  } else if (property %in% c("bulk_density", "dbovendry", "db")) {
    constraints$range <- c(0.3, 3.0)
    constraints$preserve_distribution <- TRUE
  } else if (property %in% c("ph")) {
    constraints$range <- c(2.5, 11.0)
    constraints$preserve_distribution <- TRUE
  } else if (property %in% c("wthirdbar", "wfifteenbar", "water_retention_third_bar", "water_retention_15_bar")) {
    constraints$range <- c(0, 70)
    constraints$preserve_distribution <- TRUE
  } else if (property %in% c("rfv")) {
    constraints$range <- c(0, 95)
    constraints$preserve_distribution <- TRUE
  }

  return(constraints)
}

apply_range_constraints <- function(values, constraints) {
  if (!is.null(constraints$range)) {
    values <- pmax(values, constraints$range[1])
    values <- pmin(values, constraints$range[2])
  }
  return(values)
}

correct_property_distribution <- function(adjusted_values, original_values, constraints) {
  if (all(is.na(adjusted_values)) || all(is.na(original_values))) {
    return(adjusted_values)
  }

  # Enhanced quantile mapping approach using Module 8 safe operations
  tryCatch({
    original_quantiles <- ecdf(original_values)
    adjusted_quantiles <- ecdf(adjusted_values)

    # Map adjusted values back to original distribution
    corrected_values <- quantile(original_values,
                                 probs = adjusted_quantiles(adjusted_values),
                                 na.rm = TRUE)

    return(corrected_values)
  }, error = function(e) {
    handle_workflow_error(e, "Property distribution correction", "warn")
    return(adjusted_values)
  })
}

apply_cross_property_constraints <- function(data, properties) {
  # Enhanced texture sum constraint using Module 8 validation
  texture_props <- intersect(c("sand_total", "sandtotal", "clay_total", "claytotal", "silt_total", "silttotal"),
                             properties)

  if (length(texture_props) >= 2) {
    tryCatch({
      # Normalize texture properties to sum to 100%
      for (i in seq_len(nrow(data))) {
        texture_values <- unlist(data[i, texture_props])
        texture_sum <- sum(texture_values, na.rm = TRUE)

        if (texture_sum > 0 && !is.na(texture_sum)) {
          scaling_factor <- 100 / texture_sum
          data[i, texture_props] <- texture_values * scaling_factor
        }
      }
    }, error = function(e) {
      handle_workflow_error(e, "Cross-property constraint application", "warn")
    })
  }

  return(data)
}

check_missing_values_increase <- function(original_data, integrated_data, properties) {
  results <- list()

  for (prop in properties) {
    if (prop %in% names(original_data) && prop %in% names(integrated_data)) {
      tryCatch({
        original_missing <- sum(is.na(original_data[[prop]]))
        integrated_missing <- sum(is.na(integrated_data[[prop]]))

        results[[prop]] <- list(
          original_missing = original_missing,
          integrated_missing = integrated_missing,
          increase = integrated_missing - original_missing
        )
      }, error = function(e) {
        handle_workflow_error(e, paste("Missing value check for", prop), "warn")
      })
    }
  }

  return(results)
}

# Individual adjustments (fallback approach)
apply_individual_adjustments <- function(property_matrices, gp_predictions, depths) {
  adjusted_matrices <- list()

  for (prop in names(property_matrices)) {
    current_matrix <- property_matrices[[prop]]
    gp_means <- gp_predictions[[prop]]

    if (is.null(gp_means)) {
      adjusted_matrices[[prop]] <- current_matrix
      next
    }

    # Simple scaling approach with Module 8 safety
    adjusted_matrix <- current_matrix

    for (i in 2:nrow(current_matrix)) {
      if (!is.na(gp_means[i-1]) && !is.na(gp_means[i]) && gp_means[i-1] != 0) {
        scaling_factor <- gp_means[i] / gp_means[i-1]
        # Clamp extreme scaling factors
        scaling_factor <- pmax(0.1, pmin(10, scaling_factor))
        adjusted_matrix[i, ] <- adjusted_matrix[i-1, ] * scaling_factor
      }
    }

    adjusted_matrices[[prop]] <- adjusted_matrix
  }

  return(adjusted_matrices)
}
