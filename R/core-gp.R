#' @title Gaussian Process Depth Modeling
#' @description Functions for fitting Gaussian Process models of soil property
#'   depth trends and applying GP-informed, correlation-preserving adjustments
#'   to simulated soil profiles.
#' @name gp_modeling
NULL

# ============================================================================
# 1. NRCS DATA PREPARATION FUNCTIONS
# ============================================================================

#' Prepare NRCS Training Data for GP Model Building
#'
#' @param nrcs_combined_data Combined NRCS data with horizon and component information
#' @param grouping_strategy How to group soils ("auto" for intelligent selection)
#' @param min_profiles_per_group Minimum profiles required per group (default = 3)
#' @param min_observations_per_group Minimum total observations per group (default = 15)
#' @param target_min_groups Minimum number of adequate groups desired (default = 3)
#' @param max_depth Maximum depth to include in training (default = 250 cm)
#' @param validation_config Validation configuration from the shared utilities
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Processed data frame ready for GP model building
#' @family gp-modeling
#' @export
prepare_nrcs_training_data <- function(nrcs_combined_data,
                                       grouping_strategy = "auto",
                                       min_profiles_per_group = 3,
                                       min_observations_per_group = 15,
                                       target_min_groups = 3,
                                       max_depth = DEFAULT_MAX_DEPTH_CM,
                                       validation_config = NULL,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== NRCS TRAINING DATA PREPARATION ===", category = "GPModeling")
  log_message("INFO", paste("Grouping strategy:", grouping_strategy), category = "GPModeling")
  log_message("INFO", paste("Min profiles per group:", min_profiles_per_group), category = "GPModeling")
  log_message("INFO", paste("Max depth:", max_depth, "cm"), category = "GPModeling")

  # Initial data validation
  if (is.null(validation_config)) {
    validation_config <- default_config("validation")
  }

  required_columns <- c("cokey", "hzdept_r")
  validation_results <- validate_data_quality(
    nrcs_combined_data,
    required_columns = required_columns,
    quality_thresholds = validation_config
  )

  if (!validation_results$overall_quality$validation_passed) {
    log_message("WARN", "Input data quality issues detected", category = "GPModeling")
  }

  # Standardize property names
  processed_data <- standardize_property_names(
    nrcs_combined_data,
    target_standard = "ssurgo"
  )

  # Apply unsuitable-horizon detection from the shared utilities. is_unsuitable() needs the
  # data frame passed as a named argument (not just piped in), so this is written out rather
  # than chained:
  # computing the flags against `processed_data` directly beforehand,
  # rather than leaving a now-unbound bare `.` inside mutate().
  unsuitable_horizon_flags <- is_unsuitable(processed_data, strict_mode = TRUE)
  processed_data <- processed_data |>
    dplyr::mutate(
      unsuitable_horizon = unsuitable_horizon_flags
    ) |>
    dplyr::filter(!unsuitable_horizon)

  log_message("DEBUG", paste("Filtered out", sum(processed_data$unsuitable_horizon, na.rm = TRUE), "unsuitable horizons"), category = "GPModeling")

  # Auto-select optimal grouping strategy if requested
  if (grouping_strategy == "auto") {
    log_message("INFO", "Auto-selecting optimal grouping strategy...", category = "GPModeling")
    grouping_strategy <- select_optimal_grouping(
      processed_data,
      min_profiles_per_group,
      min_observations_per_group,
      target_min_groups
    )
    log_message("INFO", paste("Selected strategy:", grouping_strategy), category = "GPModeling")
  }

  # Process depth information
  processed_data <- processed_data |>
    dplyr::mutate(
      hzdept_r = dplyr::case_when(
        !is.na(hzdept_r) ~ hzdept_r,
        !is.na(hzdepb_r) & !is.na(hzdept_l) ~ (hzdept_l + hzdepb_r) / 2,
        !is.na(hzdept_l) & !is.na(hzdepb_l) ~ (hzdept_l + hzdepb_l) / 2,
        !is.na(hzdept_h) & !is.na(hzdepb_h) ~ (hzdept_h + hzdepb_h) / 2,
        TRUE ~ as.numeric(NA)
      ),
      hzdepb_r = dplyr::if_else(is.na(hzdepb_r), hzdept_r + 10, hzdepb_r),
      depth_midpoint = (hzdept_r + hzdepb_r) / 2
    ) |>
    dplyr::filter(hzdept_r >= 0, hzdepb_r <= max_depth)

  # Use the shared safe-coalesce helper for property standardization
  property_mappings <- list(
    clay_pct = c("claytotal_r", "clay_r", "clay"),
    sand_pct = c("sandtotal_r", "sand_r", "sand"),
    silt_pct = c("silttotal_r", "silt_r", "silt"),
    pH = c("ph1to1h2o_r", "pH_r", "pH", "ph_r"),
    organic_matter = c("om_r", "organic_matter_r", "organic_matter", "om"),
    bulk_density = c("dbovendry_r", "dbthirdbar_r", "bulk_density_r", "bulk_density", "db_r"),
    cec = c("cec7_r", "cec_r", "cec"),
    awc = c("awc_r", "awc")
  )

  for (prop_name in names(property_mappings)) {
    source_cols <- property_mappings[[prop_name]]
    processed_data[[prop_name]] <- safe_coalesce(
      processed_data,
      source_cols,
      target_type = "numeric"
    )
  }

  # Create soil grouping
  processed_data <- create_soil_groups(processed_data, grouping_strategy)

  # Final data quality assessment
  processed_data <- processed_data |>
    dplyr::filter(
      !is.na(soil_group),
      soil_group != "unknown",
      !is.na(hzdept_r),
      hzdept_r >= 0
    ) |>
    dplyr::group_by(soil_group) |>
    dplyr::filter(dplyr::n() >= min_observations_per_group) |>
    dplyr::ungroup() |>
    dplyr::arrange(cokey, hzdept_r)

  # Generate comprehensive summary
  generate_gp_processing_summary(processed_data, grouping_strategy)

  return(processed_data)
}

#' Select Optimal Grouping Strategy
#'
#'
#' @param data Input NRCS data
#' @param min_profiles Minimum profiles per group
#' @param min_obs Minimum observations per group
#' @param target_groups Minimum number of adequate groups desired
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Selected grouping strategy name
#' @family gp-modeling
#' @export
select_optimal_grouping <- function(data, min_profiles, min_obs, target_groups,
                                    verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Define grouping hierarchy with validation
  grouping_strategies <- list(
    soil_series = list(strategy = "soil_series", required_col = "compname"),
    taxonomic_class = list(strategy = "taxonomic_class", required_col = "taxclname"),
    particle_size = list(strategy = "particle_size", required_col = "taxpartsize"),
    soil_grtgroup = list(strategy = "soil_grtgroup", required_col = "taxgrtgroup"),
    soil_suborder = list(strategy = "soil_suborder", required_col = "taxsuborder"),
    soil_order = list(strategy = "soil_order", required_col = "taxorder"),
    none = list(strategy = "none", required_col = NA)
  )

  # Use shared column validation
  available_strategies <- list()
  for (strat_name in names(grouping_strategies)) {
    required_col <- grouping_strategies[[strat_name]]$required_col
    if (is.na(required_col) || required_col %in% names(data)) {
      available_strategies[[strat_name]] <- grouping_strategies[[strat_name]]
    }
  }

  if (length(available_strategies) == 0) {
    stop("No grouping strategies available - missing all required columns")
  }

  log_message("DEBUG", paste("Testing", length(available_strategies), "available grouping strategies"), category = "GPModeling")

  # Test each strategy
  strategy_results <- list()
  for (strat_name in names(available_strategies)) {
    strategy <- available_strategies[[strat_name]]$strategy

    log_message("DEBUG", paste("Testing strategy:", strategy), category = "GPModeling")

    test_data <- create_test_grouping(data, strategy)

    if (nrow(test_data) == 0) {
      log_message("DEBUG", paste("Strategy", strategy, "produced no data"), category = "GPModeling")
      next
    }

    # Calculate group statistics
    group_stats <- test_data |>
      dplyr::group_by(test_group) |>
      dplyr::summarise(
        n_profiles = dplyr::n_distinct(cokey),
        n_observations = dplyr::n(),
        depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
        .groups = "drop"
      )

    # Identify adequate groups
    adequate_groups <- group_stats |>
      dplyr::filter(
        n_profiles >= min_profiles,
        n_observations >= min_obs,
        depth_range >= 20
      )

    n_adequate <- nrow(adequate_groups)
    pct_data_covered <- if (n_adequate > 0) {
      sum(adequate_groups$n_observations) / sum(group_stats$n_observations) * 100
    } else {
      0
    }

    log_message("DEBUG", paste("Strategy", strategy, ":", n_adequate, "adequate groups,",
                               round(pct_data_covered, 1), "% data coverage"), category = "GPModeling")

    strategy_results[[strategy]] <- list(
      n_adequate_groups = n_adequate,
      pct_data_covered = pct_data_covered,
      meets_criteria = n_adequate >= target_groups
    )

    # Return first strategy that meets criteria
    if (n_adequate >= target_groups) {
      log_message("INFO", paste("Strategy", strategy, "meets criteria"), category = "GPModeling")
      return(strategy)
    }
  }

  # Find best available option
  best_n_groups <- max(sapply(strategy_results, function(x) x$n_adequate_groups))
  best_strategy <- names(strategy_results)[which.max(sapply(strategy_results, function(x) x$n_adequate_groups))]

  log_message("WARN", paste("Using", best_strategy, "as best available option (", best_n_groups, "adequate groups)"), category = "GPModeling")

  return(best_strategy)
}

#' Validate Training Groups
#'
#'
#' @param processed_data Processed NRCS data
#' @param properties Properties to validate
#' @param min_profiles Minimum profiles per group
#' @param min_observations Minimum observations per group
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results list
#' @keywords internal
validate_training_groups <- function(processed_data,
                                     properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
                                     min_profiles = 3,
                                     min_observations = 15,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== TRAINING GROUP VALIDATION ===", category = "GPModeling")

  # Use the shared validation helpers
  validation_results <- list(
    overall_adequacy = TRUE,
    group_validation = list(),
    property_coverage = list(),
    recommendations = character(0)
  )

  # Overall data assessment
  total_groups <- dplyr::n_distinct(processed_data$soil_group)
  total_profiles <- dplyr::n_distinct(processed_data$cokey)
  total_observations <- nrow(processed_data)

  log_message("INFO", paste("Overall data - Groups:", total_groups, "Profiles:", total_profiles, "Observations:", total_observations), category = "GPModeling")

  # Group adequacy assessment
  group_summary <- processed_data |>
    dplyr::group_by(soil_group) |>
    dplyr::summarise(
      n_profiles = dplyr::n_distinct(cokey),
      n_observations = dplyr::n(),
      depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
      min_depth = min(hzdept_r, na.rm = TRUE),
      max_depth = max(hzdept_r, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      adequate_profiles = n_profiles >= min_profiles,
      adequate_observations = n_observations >= min_observations,
      adequate_depth_range = depth_range >= 20,
      overall_adequate = adequate_profiles & adequate_observations & adequate_depth_range
    )

  adequate_groups <- sum(group_summary$overall_adequate)

  log_message("INFO", paste("Group adequacy:", adequate_groups, "/", total_groups, "groups adequate"), category = "GPModeling")

  if (adequate_groups < 2) {
    validation_results$overall_adequacy <- FALSE
    validation_results$recommendations <- c(
      validation_results$recommendations,
      "Insufficient adequate groups for reliable GP modeling"
    )
  }

  validation_results$group_validation <- group_summary

  # Property coverage assessment
  available_properties <- intersect(properties, names(processed_data))

  if (length(available_properties) > 0) {
    for (prop in available_properties) {
      prop_coverage <- processed_data |>
        dplyr::filter(soil_group %in% group_summary$soil_group[group_summary$overall_adequate]) |>
        dplyr::group_by(soil_group) |>
        dplyr::summarise(
          total_obs = dplyr::n(),
          prop_available = sum(!is.na(.data[[prop]])),
          coverage_rate = prop_available / total_obs,
          .groups = "drop"
        )

      overall_coverage <- sum(prop_coverage$prop_available) / sum(prop_coverage$total_obs)
      groups_with_good_coverage <- sum(prop_coverage$coverage_rate >= 0.5)

      log_message("DEBUG", paste("Property", prop, ":", round(overall_coverage * 100, 1),
                                 "% coverage,", groups_with_good_coverage, "groups >50%"), category = "GPModeling")

      validation_results$property_coverage[[prop]] <- prop_coverage

      if (overall_coverage < 0.3) {
        validation_results$recommendations <- c(
          validation_results$recommendations,
          paste("Low coverage for", prop, "- consider data enrichment")
        )
      }
    }
  }

  # Log recommendations
  if (length(validation_results$recommendations) > 0) {
    log_message("WARN", "Validation recommendations available", category = "GPModeling")
    for (rec in validation_results$recommendations) {
      log_message("WARN", paste("-", rec), category = "GPModeling")
    }
  } else {
    log_message("INFO", "Training groups meet adequacy criteria", category = "GPModeling")
  }

  return(validation_results)
}

# ============================================================================
# 2. GP MODEL BUILDING FUNCTIONS# ============================================================================

#' Build Stratified GP Models
#'
#'
#' @param processed_nrcs_data Processed NRCS data from prepare_nrcs_training_data()
#' @param properties Character vector of properties to model
#' @param min_profiles_per_group Minimum profiles per group (default = 3)
#' @param min_observations_per_group Minimum observations per group (default = 15)
#' @param optimize_hyperparameters Whether to optimize GP hyperparameters (default = TRUE)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of GP models organized by property and group
#' @family gp-modeling
#' @export
fit_depth_gp_models <- function(processed_nrcs_data,
                                       properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
                                       min_profiles_per_group = 3,
                                       min_observations_per_group = 15,
                                       optimize_hyperparameters = TRUE,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== BUILDING STRATIFIED GP MODELS ===", category = "GPModeling")
  log_message("INFO", paste("Properties:", paste(properties, collapse = ", ")), category = "GPModeling")

  gp_models <- list()

  for (prop in properties) {
    log_message("INFO", paste("Building GP models for", prop), category = "GPModeling")

    if (!prop %in% names(processed_nrcs_data)) {
      log_message("WARN", paste("Property", prop, "not found in data"), category = "GPModeling")
      next
    }

    # Apply hierarchical grouping with intelligent fallback
    grouped_data <- tryCatch({
      apply_hierarchical_grouping(
        processed_nrcs_data,
        prop,
        min_profiles_per_group,
        min_observations_per_group
      )
    }, error = function(e) {
      handle_workflow_error(e, paste("Hierarchical grouping for", prop), "warn")
      return(NULL)
    })

    if (is.null(grouped_data)) {
      log_message("ERROR", paste("Failed to create groups for", prop), category = "GPModeling")
      next
    }

    # Analyze final grouping strategy
    group_summary <- grouped_data |>
      dplyr::group_by(final_group) |>
      dplyr::summarise(
        n_profiles = dplyr::n_distinct(cokey),
        n_observations = dplyr::n(),
        depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
        mean_value = mean(.data[[prop]], na.rm = TRUE),
        original_groups = dplyr::n_distinct(soil_group),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(n_observations))

    log_message("DEBUG", paste("Final grouping for", prop, ":", nrow(group_summary), "groups"), category = "GPModeling")

    # Build GP models for each final group
    group_models <- list()
    model_diagnostics <- list()

    for (group in group_summary$final_group) {
      group_data <- grouped_data |>
        dplyr::filter(final_group == group, !is.na(.data[[prop]]))

      log_message("DEBUG", paste("Building model for group:", group, "(", nrow(group_data), "observations)"), category = "GPModeling")

      # Fit individual GP model with error handling
      model_result <- tryCatch({
        fit_individual_gp_model(
          group_data,
          prop,
          optimize_hyperparameters = optimize_hyperparameters
        )
      }, error = function(e) {
        handle_workflow_error(e, paste("GP fitting for", prop, "group", group), "warn")
        return(NULL)
      })

      if (!is.null(model_result)) {
        group_models[[group]] <- model_result$model
        model_diagnostics[[group]] <- model_result$diagnostics
      }
    }

    # Store models and metadata
    gp_models[[prop]] <- list(
      type = "stratified_grouped",
      models = group_models,
      grouping_summary = group_summary,
      model_diagnostics = model_diagnostics,
      property = prop,
      n_groups = length(group_models),
      total_training_points = sum(group_summary$n_observations)
    )

    log_message("INFO", paste("Completed", prop, "- built", length(group_models), "group models"), category = "GPModeling")
  }

  # Add overall model summary
  gp_models$model_summary <- create_model_summary(gp_models, processed_nrcs_data)

  log_message("INFO", "=== STRATIFIED GP MODEL BUILDING COMPLETE ===", category = "GPModeling")

  return(new_soilSIM_gp_models(gp_models))
}

#' Fit Individual GP Model
#'
#'
#' @param data Input data for the group
#' @param property Property name to model
#' @param optimize_hyperparameters Whether to optimize hyperparameters
#' @param gp_control Passed through to `GPfit::GP_fit()`'s `control` argument (population size /
#'   iteration counts for its internal hyperparameter search). Defaults to `c(20, 10, 2)`, far
#'   below `GP_fit()`'s own default `c(200*d, 80*d, 2*d)` - see
#'   `fit_local_gp_model_single()`'s docs (`core-gp.R`) for the empirical
#'   justification (identical fitted results, ~5x faster per call). Only used when
#'   `optimize_hyperparameters = TRUE` calls through to `optimize_gp_hyperparameters()`/
#'   `k_fold_gp_cv()` (many `GP_fit()` calls per fit - 5-fold CV x 3 correlation candidates + a
#'   refit), or by the direct `GP_fit()` call when `optimize_hyperparameters = FALSE`.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List containing GP model and diagnostics
#' @family gp-modeling
#' @export
fit_individual_gp_model <- function(data, property, optimize_hyperparameters = TRUE,
                                    gp_control = c(20, 10, 2),
                                    verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  tryCatch({
    # Aggregate by depth to get smooth trends
    aggregated_data <- data |>
      dplyr::group_by(cokey, hzdept_r) |>
      dplyr::summarise(
        prop_value = mean(.data[[property]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::group_by(hzdept_r) |>
      dplyr::summarise(
        mean_value = mean(prop_value, na.rm = TRUE),
        n_profiles = dplyr::n_distinct(cokey),
        sd_value = sd(prop_value, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::filter(
        !is.na(mean_value),
        !is.infinite(mean_value),
        n_profiles >= 1
      ) |>
      dplyr::arrange(hzdept_r)

    if (nrow(aggregated_data) < 3) {
      log_message("WARN", paste("Insufficient data after aggregation for", property), category = "GPModeling")
      return(NULL)
    }

    # Prepare for GP fitting
    depths <- aggregated_data$hzdept_r
    values <- aggregated_data$mean_value

    # Check for variation
    if (var(values, na.rm = TRUE) <= 0) {
      log_message("WARN", paste("No variation in", property, "values"), category = "GPModeling")
      return(NULL)
    }

    # Standardize depths to [0,1]
    depth_min <- min(depths)
    depth_max <- max(depths)
    depth_range <- depth_max - depth_min

    if (depth_range <= 0) {
      log_message("WARN", paste("No depth variation for", property), category = "GPModeling")
      return(NULL)
    }

    scaled_depths <- (depths - depth_min) / depth_range

    # Fit GP model with optimization
    if (optimize_hyperparameters) {
      gp_model <- optimize_gp_hyperparameters(scaled_depths, values, gp_control = gp_control)
    } else {
      gp_model <- GPfit::GP_fit(X = as.matrix(scaled_depths), Y = values, control = gp_control)
    }

    # Calculate diagnostics
    diagnostics <- calculate_model_diagnostics(gp_model, aggregated_data, property)

    # Return model with metadata
    result <- list(
      model = list(
        gp_model = gp_model,
        depth_scaling = list(
          min = depth_min,
          max = depth_max,
          range = depth_range
        ),
        training_data = aggregated_data,
        n_training_points = nrow(aggregated_data),
        property = property
      ),
      diagnostics = diagnostics
    )

    return(result)

  }, error = function(e) {
    handle_workflow_error(e, paste("GP fitting for", property), "warn")
    return(NULL)
  })
}

# ============================================================================
# 3. MULTIVARIATE CORRELATION PRESERVATION
# ============================================================================

#' Adjust Multiple Soil Properties While Preserving Correlations
#'
#' Adjusts multiple simulated soil properties
#' simultaneously to follow GP-predicted trends while preserving within-depth correlations.
#'
#' @param simulated_list A named list of matrices, where each matrix contains simulated
#'   values for one property (rows = depths, columns = simulations)
#' @param gp_models A named list of fitted GP models for each property
#' @param depths A numeric vector of depth values
#' @param primary_property Character string specifying which property to use as the
#'   "reference" for correlation preservation (default: first property in list)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#'
#' @return A named list of adjusted matrices with preserved correlations
#'
#' @family gp-modeling
#' @export
adjust_simulation_depthwise <- function(simulated_list, gp_models, depths,
                                             primary_property = NULL,
                                             verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== MULTIVARIATE GP DEPTH ADJUSTMENT ===", category = "GPModeling")

  # Input validation
  if (length(simulated_list) == 0) {
    stop("simulated_list cannot be empty")
  }

  property_names <- names(simulated_list)
  if (is.null(property_names) || any(property_names == "")) {
    stop("simulated_list must have named elements")
  }

  # Check GP models match properties
  if (!all(property_names %in% names(gp_models))) {
    stop("All properties in simulated_list must have corresponding GP models")
  }

  # Get dimensions from first property
  n_depths <- nrow(simulated_list[[1]])
  n_sims <- ncol(simulated_list[[1]])

  # Validate all properties have same dimensions
  for (prop in property_names) {
    if (nrow(simulated_list[[prop]]) != n_depths ||
        ncol(simulated_list[[prop]]) != n_sims) {
      stop(paste("Property", prop, "has inconsistent dimensions"))
    }
  }

  # Check depths match
  if (length(depths) != n_depths) {
    stop("Length of depths must match number of rows in property matrices")
  }

  # Set primary property
  if (is.null(primary_property)) {
    primary_property <- property_names[1]
    log_message("INFO", paste("Using", primary_property, "as primary property for correlation preservation"), category = "GPModeling")
  }

  if (!primary_property %in% property_names) {
    stop("primary_property must be one of the properties in simulated_list")
  }

  # Get GP predictions for all properties
  gp_predictions <- list()
  for (prop in property_names) {
    tryCatch({
      if (is.list(gp_models[[prop]]) && "gp_model" %in% names(gp_models[[prop]])) {
        # Handle nested GP model structure
        gp_predictions[[prop]] <- predict_gp_depth_trends(gp_models[[prop]], depths)
      } else {
        # Handle direct GP model
        gp_predictions[[prop]] <- GPfit::predict.GP(gp_models[[prop]],
                                                    xnew = as.matrix(depths))$Y_hat
      }
    }, error = function(e) {
      handle_workflow_error(e, paste("GP prediction for", prop), "warn")
      # Use mean as fallback
      gp_predictions[[prop]] <<- rowMeans(simulated_list[[prop]], na.rm = TRUE)
    })
  }

  # Establish reference quantile ordering
  primary_matrix <- simulated_list[[primary_property]]
  surface_values <- primary_matrix[1, ]

  # Create ECDF and get quantiles for surface (reference depth)
  surface_ecdf <- ecdf(surface_values)
  reference_quantiles <- surface_ecdf(surface_values)

  # Apply consistent adjustment to all properties
  adjusted_list <- list()

  for (prop in property_names) {

    log_message("DEBUG", paste("Adjusting property:", prop), category = "GPModeling")

    current_matrix <- simulated_list[[prop]]
    gp_means <- gp_predictions[[prop]]

    # Initialize with original values
    adjusted_matrix <- current_matrix

    # Apply depth-wise adjustment using SAME quantile ordering for all properties
    for (i in 2:n_depths) {

      # Get GP trend ratio with safety checks
      if (is.na(gp_means[i-1]) || is.na(gp_means[i]) || gp_means[i-1] == 0) {
        gp_ratio <- 1  # No adjustment if invalid
      } else {
        gp_ratio <- gp_means[i] / gp_means[i - 1]
      }

      # Get previous and current simulated values
      prev_values <- adjusted_matrix[i - 1, ]
      curr_values <- current_matrix[i, ]  # Use original current values

      # quantile() accepts a probs VECTOR, so this is one call for all replicates at once.
      # unname() drops quantile()'s "37%"-style names.
      quantile_values <- unname(stats::quantile(curr_values, probs = reference_quantiles, na.rm = TRUE))
      adjusted_curr <- quantile_values + (prev_values * gp_ratio - quantile_values)

      # Correct to maintain original distribution shape
      if (var(curr_values, na.rm = TRUE) > 0) {
        ecdf_adjusted <- ecdf(adjusted_curr)
        # Map back to original distribution while preserving order
        corrected_values <- quantile(curr_values,
                                     probs = ecdf_adjusted(adjusted_curr),
                                     na.rm = TRUE)
        adjusted_matrix[i, ] <- corrected_values
      } else {
        # If no variation, keep original values
        adjusted_matrix[i, ] <- curr_values
      }
    }

    adjusted_list[[prop]] <- adjusted_matrix
  }

  # Verify correlation preservation
  log_message("INFO", "Verifying correlation preservation...", category = "GPModeling")

  for (depth_idx in 1:min(5, n_depths)) {  # Check first 5 depths
    original_cors <- tryCatch({
      cor(sapply(simulated_list, function(x) x[depth_idx, ]), use = "complete.obs")
    }, error = function(e) matrix(1, nrow = length(property_names), ncol = length(property_names)))

    adjusted_cors <- tryCatch({
      cor(sapply(adjusted_list, function(x) x[depth_idx, ]), use = "complete.obs")
    }, error = function(e) matrix(1, nrow = length(property_names), ncol = length(property_names)))

    max_cor_diff <- max(abs(original_cors - adjusted_cors), na.rm = TRUE)

    log_message("DEBUG", paste("Depth", depth_idx, "- Max correlation difference:",
                               round(max_cor_diff, 4)), category = "GPModeling")

    if (max_cor_diff > 0.1) {
      log_message("WARN", paste("Large correlation change at depth", depth_idx), category = "GPModeling")
    }
  }

  log_message("INFO", "Multivariate GP adjustment completed", category = "GPModeling")

  return(adjusted_list)
}

#' Validate Correlation Preservation
#'
#' Validates that correlation structure is preserved after depth-trend adjustment.
#'
#' @param original_list Original simulated property list
#' @param adjusted_list Adjusted property list
#' @param depths Depth vector
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#'
#' @keywords internal
validate_correlation_preservation <- function(original_list, adjusted_list, depths,
                                              verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== CORRELATION PRESERVATION VALIDATION ===", category = "GPModeling")

  property_names <- names(original_list)
  n_depths <- length(depths)

  correlation_summary <- data.frame(
    depth_idx = 1:n_depths,
    depth_value = depths,
    max_cor_difference = NA,
    mean_cor_difference = NA
  )

  for (i in 1:n_depths) {
    # Original correlations
    orig_data <- tryCatch({
      sapply(original_list, function(x) x[i, ])
    }, error = function(e) {
      log_message("WARN", paste("Error extracting original data at depth", i), category = "GPModeling")
      return(NULL)
    })

    if (is.null(orig_data)) {
      next
    }

    orig_cors <- tryCatch({
      cor(orig_data, use = "complete.obs")
    }, error = function(e) {
      log_message("WARN", paste("Error calculating original correlations at depth", i), category = "GPModeling")
      return(NULL)
    })

    # Adjusted correlations
    adj_data <- tryCatch({
      sapply(adjusted_list, function(x) x[i, ])
    }, error = function(e) {
      log_message("WARN", paste("Error extracting adjusted data at depth", i), category = "GPModeling")
      return(NULL)
    })

    if (is.null(adj_data)) {
      next
    }

    adj_cors <- tryCatch({
      cor(adj_data, use = "complete.obs")
    }, error = function(e) {
      log_message("WARN", paste("Error calculating adjusted correlations at depth", i), category = "GPModeling")
      return(NULL)
    })

    if (!is.null(orig_cors) && !is.null(adj_cors)) {
      # Calculate differences
      cor_diffs <- abs(orig_cors - adj_cors)

      correlation_summary$max_cor_difference[i] <- max(cor_diffs, na.rm = TRUE)
      correlation_summary$mean_cor_difference[i] <- mean(cor_diffs[upper.tri(cor_diffs)], na.rm = TRUE)
    }
  }

  # Summary statistics
  overall_max_diff <- max(correlation_summary$max_cor_difference, na.rm = TRUE)
  overall_mean_diff <- mean(correlation_summary$mean_cor_difference, na.rm = TRUE)

  log_message("INFO", paste("Overall correlation preservation - Max diff:", round(overall_max_diff, 4),
                            "Mean diff:", round(overall_mean_diff, 4)), category = "GPModeling")

  return(correlation_summary)
}

#' Validate Joint Depth x Property Correlation Structure
#'
#' Diagnostic helper for the joint depth x property correlation structure. Unlike
#' `validate_correlation_preservation()`, which only checks cross-property correlation and (via
#' its caller `adjust_simulation_depthwise()`) only at the first 5 depths, this reports BOTH
#' halves of the joint structure - cross-property correlation at every depth, and depth-lag
#' (across-depth, within-property) correlation for every property - against explicit target
#' matrices when supplied. This is the acceptance test any vertical-correlation method
#' (the existing GP-quantile-retrofit path or a future joint-copula replacement) should be judged
#' against: does the SAME simulated output actually achieve both the intended property correlation
#' matrix and the intended depth correlation structure, simultaneously.
#'
#' @param simulated_list A named list of matrices (rows = depths, columns = simulations), the same
#'   shape consumed/produced by `preserve_correlation_structure()`/
#'   `adjust_simulation_depthwise()`.
#' @param target_property_corr Optional k x k target property-correlation matrix (row/column names
#'   matching `names(simulated_list)`, or unnamed matching order). If supplied, per-depth and
#'   overall max absolute deviation from this target is reported.
#' @param target_depth_corr Optional n_depths x n_depths target depth-correlation matrix. If
#'   supplied, per-property and overall max absolute deviation from this target is reported.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE}
#'   - quiet). See \code{set_verbose_logging()}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{achieved_property_correlation}{List (one per depth) of achieved k x k property
#'       correlation matrices.}
#'     \item{property_correlation_max_diff}{Numeric vector (one per depth) of max absolute
#'       deviation from `target_property_corr`, or all `NA` if no target supplied.}
#'     \item{achieved_depth_correlation}{Named list (one per property) of achieved n_depths x
#'       n_depths depth correlation matrices.}
#'     \item{depth_correlation_max_diff}{Named numeric vector (one per property) of max absolute
#'       deviation from `target_depth_corr`, or all `NA` if no target supplied.}
#'     \item{overall_property_max_diff, overall_depth_max_diff}{Single worst-case summary values
#'       across all depths/properties respectively.}
#'   }
#'
#' @export
validate_joint_correlation_structure <- function(simulated_list,
                                                  target_property_corr = NULL,
                                                  target_depth_corr = NULL,
                                                  verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (length(simulated_list) == 0) {
    stop("simulated_list cannot be empty")
  }

  property_names <- names(simulated_list)
  n_depths <- nrow(simulated_list[[1]])

  log_message("INFO", "=== JOINT DEPTH x PROPERTY CORRELATION VALIDATION ===", category = "GPModeling")

  # --- Cross-property correlation, at EVERY depth (not just the first 5) ---
  achieved_property_correlation <- vector("list", n_depths)
  property_correlation_max_diff <- rep(NA_real_, n_depths)

  for (i in seq_len(n_depths)) {
    depth_data <- tryCatch({
      sapply(simulated_list, function(x) x[i, ])
    }, error = function(e) NULL)

    cors <- if (!is.null(depth_data)) {
      tryCatch(cor(depth_data, use = "complete.obs"), error = function(e) NULL)
    } else {
      NULL
    }

    achieved_property_correlation[[i]] <- cors

    if (!is.null(cors) && !is.null(target_property_corr)) {
      property_correlation_max_diff[i] <- max(abs(cors - target_property_corr), na.rm = TRUE)
    }
  }

  # --- Depth-lag correlation, for EVERY property (not checked by
  # validate_correlation_preservation() at all) ---
  achieved_depth_correlation <- vector("list", length(property_names))
  names(achieved_depth_correlation) <- property_names
  depth_correlation_max_diff <- rep(NA_real_, length(property_names))
  names(depth_correlation_max_diff) <- property_names

  for (prop in property_names) {
    m <- simulated_list[[prop]]
    dc <- tryCatch(cor(t(m), use = "complete.obs"), error = function(e) NULL)
    achieved_depth_correlation[[prop]] <- dc

    if (!is.null(dc) && !is.null(target_depth_corr)) {
      depth_correlation_max_diff[prop] <- max(abs(dc - target_depth_corr), na.rm = TRUE)
    }
  }

  overall_property_max_diff <- if (all(is.na(property_correlation_max_diff))) {
    NA_real_
  } else {
    max(property_correlation_max_diff, na.rm = TRUE)
  }

  overall_depth_max_diff <- if (all(is.na(depth_correlation_max_diff))) {
    NA_real_
  } else {
    max(depth_correlation_max_diff, na.rm = TRUE)
  }

  log_message("INFO", paste("Property correlation - overall max diff:", round(overall_property_max_diff, 4),
                            "| Depth correlation - overall max diff:", round(overall_depth_max_diff, 4)),
              category = "GPModeling")

  list(
    achieved_property_correlation = achieved_property_correlation,
    property_correlation_max_diff = property_correlation_max_diff,
    achieved_depth_correlation = achieved_depth_correlation,
    depth_correlation_max_diff = depth_correlation_max_diff,
    overall_property_max_diff = overall_property_max_diff,
    overall_depth_max_diff = overall_depth_max_diff
  )
}

# ============================================================================
# 3b. DEPTH CORRELATION KERNEL
# ============================================================================
#
# These two functions build an explicit depth-correlation matrix R_depth from a real-distance
# kernel whose length-scale is taken from the already-fitted, already-cross-validated GPfit
# model for each property (fit_local_gp_model_single()). R_depth is consumed by the
# joint copula sampler as one half of a Kronecker-separable (depth x property) covariance.

#' Extract a Real-Units Depth Length-Scale from a Fitted GP Model
#'
#' `GPfit::GP_fit()` fits a power-exponential (or Matern) correlation function
#' `R(h) = corr(0, h; beta, correlation_param)` over depths already rescaled to `[0,1]` by
#' `fit_local_gp_model_single()`. This reuses that fitted correlation function exactly (via
#' `GPfit::corr_matrix()`, the same internal machinery `GP_fit()`/`predict.GP()` use) to solve for
#' the distance `h` (in scaled `[0,1]` units) at which the correlation first drops to `target_corr`
#' (default `exp(-1)`, the conventional "range" definition for both exponential and Matern
#' kernels), then rescales `h` back to real depth units (cm) using the same `depth_scaling`
#' `fit_local_gp_model_single()` stored alongside the model - so no separate/new length-scale
#' estimation step is introduced; the value is entirely derived from the already-validated GP fit.
#'
#' @param gp_model_list Either the full list returned by `fit_local_gp_model_single()`
#'   (`list(gp_model, depth_scaling, ...)`) or a raw `GPfit`-classed `"GP"` object (in which case
#'   `depth_scaling` must be supplied separately, or the result stays in scaled `[0,1]` units).
#' @param depth_scaling Optional `list(min=, max=, range=)` as stored by
#'   `fit_local_gp_model_single()`; only needed when `gp_model_list` is a raw `GP` object without
#'   its own `$depth_scaling`. If neither is available, the returned length-scale stays in scaled
#'   `[0,1]` units.
#' @param target_corr The correlation value defining "range" (default `exp(-1) ~= 0.368`).
#'
#' @return A single numeric length-scale in real depth units (or scaled `[0,1]` units if no
#'   `depth_scaling` is available), or `NA_real_` if `gp_model_list` doesn't contain a usable
#'   fitted model (mirrors `fit_local_gp_model_single()`'s own `NULL`-on-failure contract - callers
#'   should treat `NA_real_` the same way they'd treat a `NULL` GP model).
#'
#' @family gp-modeling
#' @export
extract_depth_length_scale <- function(gp_model_list, depth_scaling = NULL, target_corr = exp(-1)) {

  gp_model <- if (is.list(gp_model_list) && "gp_model" %in% names(gp_model_list)) {
    gp_model_list$gp_model
  } else {
    gp_model_list
  }

  if (is.null(depth_scaling) && is.list(gp_model_list) && !is.null(gp_model_list$depth_scaling)) {
    depth_scaling <- gp_model_list$depth_scaling
  }

  if (is.null(gp_model) || is.null(gp_model$beta) || is.null(gp_model$correlation_param)) {
    return(NA_real_)
  }

  corr_at <- function(h) {
    if (h <= 0) return(1)
    GPfit::corr_matrix(X = matrix(c(0, h), ncol = 1), beta = gp_model$beta,
                       corr = gp_model$correlation_param)[1, 2]
  }

  # GPfit's supported correlation families (power-exponential, Matern) are monotonically
  # non-increasing in distance, so a bracket search + uniroot() reliably finds the crossing point
  # without needing to hand-derive a closed-form inverse per kernel family (which would need to be
  # re-derived if GPfit's correlation() options ever change).
  upper <- 1
  tries <- 0
  while (tries < 10 && is.finite(corr_at(upper)) && corr_at(upper) > target_corr) {
    upper <- upper * 2
    tries <- tries + 1
  }

  scaled_length_scale <- tryCatch({
    stats::uniroot(function(h) corr_at(h) - target_corr, interval = c(1e-8, upper))$root
  }, error = function(e) NA_real_)

  if (is.na(scaled_length_scale)) {
    return(NA_real_)
  }

  if (!is.null(depth_scaling) && !is.null(depth_scaling$range) &&
      is.finite(depth_scaling$range) && depth_scaling$range > 0) {
    return(scaled_length_scale * depth_scaling$range)
  }

  scaled_length_scale
}

#' Build a Depth Correlation Kernel Matrix
#'
#' Constructs an `n_depths x n_depths` correlation matrix from real depth *distances* (not just
#' adjacent-lag steps), via an exponential or Matern kernel with a given `length_scale`.
#' This is the depth half of a Kronecker-separable `R_depth (x) R_property` joint covariance
#' (see `sample_joint_depth_property_copula()`).
#'
#' @param depths Numeric vector of depth values (real units, e.g. cm; need not be evenly spaced or
#'   sorted - internally sorted for the `boundary_distinctness` gating pass, then mapped back to
#'   `depths`' original order for the returned matrix).
#' @param length_scale A single positive numeric length-scale (real units matching `depths`),
#'   typically from `extract_depth_length_scale()`. A non-finite or non-positive value degrades
#'   gracefully to the identity matrix (no depth correlation) rather than erroring, matching this
#'   package's established fallback pattern for missing/invalid correlation inputs (see
#'   `simulate_cokey_generalized()`'s pooled-matrix fallback).
#' @param kernel `"exponential"` (default) or `"matern"`.
#' @param nu Matern smoothness parameter (only used when `kernel = "matern"`); default `1.5`.
#' @param boundary_distinctness Optional numeric vector, same length as `depths`, of OSD-derived
#'   `bound_sd` values (see `attach_osd_boundary_distinctness()`) - one value per depth, interpreted as
#'   the distinctness of the boundary immediately ABOVE that depth (in sorted-ascending order; the
#'   value for the shallowest depth is unused, since there is no boundary above the top of the
#'   profile). Smaller `bound_sd` (e.g. `aqp::hzDistinctnessCodeToOffset()`'s `abrupt` ~= 1,
#'   `clear` ~= 2.5) means a sharper real discontinuity and more strongly suppresses correlation
#'   across that boundary; larger `bound_sd` (`gradual` ~= 7.5, `diffuse` ~= 10) leaves the plain
#'   distance-decay kernel almost untouched. `NULL` (default) or all-`NA` skips gating entirely,
#'   reproducing the plain distance-decay kernel exactly.
#' @param distinctness_range Named `c(min=, max=)` giving the `bound_sd` values that map to the
#'   strongest (`min_gate_weight`) and weakest (`1`, no suppression) gating respectively. Default
#'   `c(min=1, max=10)` spans `aqp::hzDistinctnessCodeToOffset()`'s actual `abrupt`-to-`diffuse`
#'   range for the codes `infill_missing_distinctness()` ever assigns.
#' @param min_gate_weight The gating floor (default `0.05`) - even the sharpest boundary leaves a
#'   small residual correlation rather than exactly zero, both physically (a truly zero-correlation
#'   entry is rarely justified) and numerically (helps `ensure_positive_definite_matrix()`'s
#'   eigenvalue repair stay well-conditioned).
#'
#' @return A symmetric, unit-diagonal, positive-definite `length(depths) x length(depths)`
#'   correlation matrix (repaired via `ensure_positive_definite_matrix()` if numerical
#'   floating-point asymmetry from the distance-matrix/gating construction pushes it slightly off).
#'
#' @keywords internal
build_depth_correlation_kernel <- function(depths, length_scale,
                                           kernel = c("exponential", "matern"), nu = 1.5,
                                           boundary_distinctness = NULL,
                                           distinctness_range = c(min = 1, max = 10),
                                           min_gate_weight = 0.05) {
  kernel <- match.arg(kernel)
  n <- length(depths)

  if (n < 1) {
    stop("depths must have at least one element")
  }

  R <- if (!is.finite(length_scale) || length_scale <= 0) {
    diag(n)
  } else {
    dist_mat <- as.matrix(stats::dist(depths, method = "euclidean"))
    dimnames(dist_mat) <- NULL  # stats::dist() labels dims "1","2",... - not meaningful here
    scaled_dist <- dist_mat / length_scale

    if (kernel == "exponential") {
      exp(-scaled_dist)
    } else {
      # Matern correlation, parameterized directly by length_scale (not GPfit's internal
      # beta/power-exponential form - a standalone, more standard geostatistical kernel since
      # this function's `length_scale` is already in real depth units by the time it's called).
      d <- scaled_dist
      d[d == 0] <- .Machine$double.eps  # avoid 0/0 in the Bessel term; true limit at d=0 is 1
      const <- (2^(1 - nu)) / gamma(nu)
      val <- const * (sqrt(2 * nu) * d)^nu * besselK(sqrt(2 * nu) * d, nu)
      diag(val) <- 1
      val
    }
  }

  if (!is.null(boundary_distinctness) && n > 1 && !all(is.na(boundary_distinctness))) {
    if (length(boundary_distinctness) != n) {
      stop("boundary_distinctness must have the same length as depths")
    }

    ord <- order(depths)
    sorted_bound_sd <- boundary_distinctness[ord]

    # One gate weight per adjacent boundary in sorted-depth order (n-1 of them): bound_sd for
    # sorted depth k (k = 2..n) is the distinctness of the boundary immediately above it, i.e.
    # between sorted depths k-1 and k. Missing/NA bound_sd degrades to weight 1 (no gating for
    # that boundary), the same graceful-fallback contract established throughout this package for
    # missing correlation-relevant inputs.
    raw <- sorted_bound_sd[-1]
    range_span <- distinctness_range["max"] - distinctness_range["min"]
    gate_weights <- (raw - distinctness_range["min"]) / range_span
    gate_weights[is.na(gate_weights)] <- 1
    gate_weights <- pmin(1, pmax(min_gate_weight, gate_weights))

    # Compound ALL boundary weights strictly between two depths (not just adjacent pairs) via
    # cumulative log-sums, so a sharp boundary anywhere along the path suppresses correlation
    # between the depths on either side of it, however many other depths lie between them.
    # log_cum[m] = sum of log(gate_weights) for every boundary at or above sorted depth m
    # (log_cum[1] = 0, the shallowest depth has no boundary above it); the gate between sorted
    # depths i and j is then exp(-|log_cum[i] - log_cum[j]|) = product of the weights strictly
    # between them.
    log_cum <- c(0, cumsum(log(gate_weights)))
    gate_matrix <- exp(-abs(outer(log_cum, log_cum, `-`)))

    R_sorted <- R[ord, ord, drop = FALSE] * gate_matrix
    R[ord, ord] <- R_sorted
  }

  R <- (R + t(R)) / 2
  diag(R) <- 1

  ensure_positive_definite_matrix(R)
}

# ============================================================================
# 4. SOIL PROPERTY SIMULATION
# ============================================================================

#' Select a Correlation Matrix and Property Set for Single-Cokey Simulation
#'
#' Chooses a single flat correlation matrix (and its associated property set)
#' to hand to [simulate_monte_carlo()] for one cokey's data,
#' from either a plain matrix or a genhz-keyed list of matrices.
#'
#' @param cokey_data Simulation data for a single cokey.
#' @param correlation_matrices A matrix, a genhz-keyed named list of
#'   matrices, or `NULL`.
#' @param txt_correlation_matrices Same shape as `correlation_matrices`, for
#'   texture properties; used as a fallback when `correlation_matrices` is
#'   unusable.
#'
#' @return List with `properties` (character vector) and `correlation_matrix`
#'   (a matrix or `NULL`, letting [simulate_monte_carlo()]
#'   estimate one empirically).
select_simulation_correlation_matrix <- function(cokey_data, correlation_matrices, txt_correlation_matrices) {

  pick_flat_matrix <- function(candidate) {
    if (is.null(candidate)) return(NULL)
    if (is.matrix(candidate)) return(candidate)
    if (is.list(candidate) && length(candidate) > 0) {
      if ("genhz" %in% names(cokey_data) && any(!is.na(cokey_data$genhz))) {
        genhz_counts <- table(cokey_data$genhz)
        modal_genhz <- names(genhz_counts)[which.max(genhz_counts)]
        if (modal_genhz %in% names(candidate) && is.matrix(candidate[[modal_genhz]])) {
          return(candidate[[modal_genhz]])
        }
      }
      first_matrix <- Filter(is.matrix, candidate)
      if (length(first_matrix) > 0) return(first_matrix[[1]])
    }
    NULL
  }

  chosen_matrix <- pick_flat_matrix(correlation_matrices)
  if (is.null(chosen_matrix)) {
    chosen_matrix <- pick_flat_matrix(txt_correlation_matrices)
  }

  if (!is.null(chosen_matrix) && !is.null(rownames(chosen_matrix))) {
    properties <- intersect(rownames(chosen_matrix), auto_detect_soil_properties(cokey_data))
  } else {
    properties <- auto_detect_soil_properties(cokey_data)
    chosen_matrix <- NULL
  }

  if (length(properties) > 0 && !is.null(chosen_matrix)) {
    chosen_matrix <- chosen_matrix[properties, properties, drop = FALSE]
  }

  list(properties = properties, correlation_matrix = chosen_matrix)
}

#' Flatten a Monte Carlo Simulation Array to Long Format
#'
#' Converts the `[horizon, property, realization]` array returned by
#' [simulate_monte_carlo()] (as `result$simulation_data`) into a
#' long-format data frame with one row per horizon-realization combination,
#' matching the shape [apply_local_gp_adjustments()], [apply_nrcs_trend_adjustments()],
#' and [convert_to_property_matrices()] expect (a `simulation_number` column
#' plus one column per property, alongside the horizon metadata).
#'
#' @param sim_array A `[horizon, property, realization]` array as produced by
#'   [simulate_correlated_properties()].
#' @param cokey_data The original per-horizon data frame passed in as
#'   `soil_data` - must have exactly one row per `sim_array` horizon index,
#'   in the same order.
#'
#' @return Long-format data frame.
flatten_simulation_array_to_long <- function(sim_array, cokey_data) {

  dims <- dim(sim_array)
  n_horizons <- dims[1]
  n_realizations <- dims[3]
  properties <- dimnames(sim_array)[[2]]

  if (nrow(cokey_data) != n_horizons) {
    stop("cokey_data must have exactly one row per horizon in sim_array (",
         n_horizons, " expected, got ", nrow(cokey_data), ")")
  }

  horizon_idx <- rep(seq_len(n_horizons), times = n_realizations)
  realization_idx <- rep(seq_len(n_realizations), each = n_horizons)

  long_df <- cokey_data[horizon_idx, , drop = FALSE]
  rownames(long_df) <- NULL
  long_df$simulation_number <- realization_idx

  for (p in seq_along(properties)) {
    # sim_array[, p, ] is an [n_horizons x n_realizations] matrix; as.vector()
    # flattens it column-major (horizon varies fastest within a realization),
    # which is exactly the order horizon_idx/realization_idx were built in.
    long_df[[properties[p]]] <- as.vector(sim_array[, p, ])
  }

  long_df
}

#' Soil Property Simulation with NRCS GP Models and Cholesky Correlations
#'
#' Combines triangular-distribution + Cholesky-correlated sampling with NRCS GP
#' models for realistic depth trends.
#'
#' @param target_cokey Character string of the cokey to simulate
#' @param nrcs_gp_models Your fitted NRCS GP models (optional)
#' @param cokey_mapping Mapping from match_simulated_soils_to_gp_models (optional)
#' @param sim_data Your processed simulation data for this cokey
#' @param correlation_matrices Your existing correlation matrices by genetic horizon
#' @param txt_correlation_matrices Your existing texture correlation matrices
#' @param n_simulations Number of Monte Carlo realizations
#' @param use_nrcs_gp Whether to use NRCS GP models (TRUE) or local GP models (FALSE)
#' @param preserve_correlations Whether to preserve within-depth correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#'
#' @return simulation results with realistic depth trends and preserved correlations
#' @keywords internal
simulate_soil_properties <- function(target_cokey,
                                              nrcs_gp_models = NULL,
                                              cokey_mapping = NULL,
                                              sim_data,
                                              correlation_matrices,
                                              txt_correlation_matrices,
                                              n_simulations = 100,
                                              use_nrcs_gp = !is.null(nrcs_gp_models),
                                              preserve_correlations = TRUE,
                                              verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== SOIL PROPERTY SIMULATION ===", category = "GPModeling")
  log_message("INFO", paste("Target cokey:", target_cokey), category = "GPModeling")
  log_message("INFO", paste("Simulations:", n_simulations, "Use NRCS GP:", use_nrcs_gp, "Preserve correlations:", preserve_correlations), category = "GPModeling")

  # Step 1: Generate initial correlated simulations using existing approach
  log_message("INFO", "Step 1: Generating correlated simulations...", category = "GPModeling")

  # Filter data for target cokey
  cokey_data <- sim_data |> dplyr::filter(cokey == target_cokey)

  if (nrow(cokey_data) == 0) {
    log_message("ERROR", paste("No data found for target cokey:", target_cokey), category = "GPModeling")
    return(data.frame())
  }

  # Generate initial simulations using the package's real Monte Carlo engine.
  # Choose a single correlation matrix to hand to simulate_monte_carlo():
  # correlation_matrices/txt_correlation_matrices are documented as
  # "existing correlation matrices by genetic horizon" (genhz-keyed lists),
  # but simulate_monte_carlo() takes one flat matrix per call -
  # use the matrix for cokey_data's most common genhz as a reasonable
  # single-cokey simplification (a cokey's horizons are usually dominated by
  # one or two genhz groups). Falls back to auto-detected properties and no
  # explicit correlation matrix (empirical estimation) if no usable
  # genhz-keyed or flat matrix is supplied.
  initial_simulations <- tryCatch({
    chosen <- select_simulation_correlation_matrix(cokey_data, correlation_matrices, txt_correlation_matrices)

    mc_result <- simulate_monte_carlo(
      soil_data = cokey_data,
      properties = chosen$properties,
      correlation_matrix = chosen$correlation_matrix,
      n_realizations = n_simulations
    )

    flatten_simulation_array_to_long(mc_result$simulation_data, cokey_data)
  }, error = function(e) {
    handle_workflow_error(e, "Initial simulation generation", "stop")
  })

  if (nrow(initial_simulations) == 0) {
    log_message("WARN", paste("No simulations generated for cokey:", target_cokey), category = "GPModeling")
    return(data.frame())
  }

  log_message("INFO", paste("Generated", nrow(initial_simulations), "initial simulation rows"), category = "GPModeling")

  # Step 2: Apply GP-based depth adjustments
  enhanced_simulations <- if (use_nrcs_gp && !is.null(nrcs_gp_models) && !is.null(cokey_mapping)) {
    log_message("INFO", "Step 2: Applying NRCS GP depth adjustments...", category = "GPModeling")

    # Find which NRCS GP model group to use
    model_group <- cokey_mapping$gp_model_group[cokey_mapping$sim_cokey == target_cokey]
    if (length(model_group) == 0 || is.na(model_group)) {
      model_group <- "fallback_general"
      log_message("INFO", paste("Using fallback model group:", model_group), category = "GPModeling")
    } else {
      log_message("INFO", paste("Using NRCS model group:", model_group), category = "GPModeling")
    }

    # Apply NRCS GP adjustments with correlation preservation
    apply_nrcs_gp_adjustments_with_correlations(
      initial_simulations,
      nrcs_gp_models,
      model_group,
      preserve_correlations
    )
  } else {
    log_message("INFO", "Step 2: Applying local GP depth adjustments...", category = "GPModeling")

    # Use local GP adjustment approach
    apply_local_gp_adjustments_with_correlations(
      initial_simulations,
      preserve_correlations
    )
  }

  log_message("INFO", paste("Simulations completed with", nrow(enhanced_simulations), "rows"), category = "GPModeling")

  # Step 3: Add validation and summary statistics
  log_message("INFO", "Step 3: Validating results...", category = "GPModeling")

  validation_results <- validate_enhanced_simulations(
    initial_simulations,
    enhanced_simulations,
    preserve_correlations
  )

  # Add validation results as attributes
  attr(enhanced_simulations, "validation") <- validation_results
  attr(enhanced_simulations, "n_simulations") <- n_simulations
  attr(enhanced_simulations, "use_nrcs_gp") <- use_nrcs_gp
  attr(enhanced_simulations, "preserve_correlations") <- preserve_correlations

  log_message("INFO", "=== SIMULATION COMPLETE ===", category = "GPModeling")

  return(enhanced_simulations)
}

# ============================================================================
# 5. SOIL MATCHING AND MODEL APPLICATION FUNCTIONS
# ============================================================================

#' Match Soils to GP Models
#'
#'
#' @param simulated_cokeys Vector of cokeys from simulation data
#' @param nrcs_combined_data Original NRCS data used for GP training
#' @param gp_models Fitted GP models from fit_depth_gp_models()
#' @param property Property name for matching (default = "clay_pct")
#' @param matching_strategy Matching approach (default = "exact_cokey")
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data frame mapping cokeys to GP model groups
#' @keywords internal
match_soils_to_gp_models <- function(simulated_cokeys,
                                     nrcs_combined_data,
                                     gp_models,
                                     property = "clay_pct",
                                     matching_strategy = "exact_cokey",
                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== MATCHING SIMULATED SOILS TO GP MODELS ===", category = "GPModeling")
  log_message("INFO", paste("Property:", property, "Strategy:", matching_strategy, "Simulated cokeys:", length(simulated_cokeys)), category = "GPModeling")

  if (!property %in% names(gp_models)) {
    stop("Property ", property, " not found in GP models")
  }

  # Get available model groups
  if (gp_models[[property]]$type == "stratified_grouped") {
    available_groups <- names(gp_models[[property]]$models)
  } else {
    # For pooled models, all soils use the same model
    return(data.frame(
      sim_cokey = simulated_cokeys,
      gp_model_group = "pooled_model",
      stringsAsFactors = FALSE
    ))
  }

  log_message("INFO", paste("Available model groups:", length(available_groups)), category = "GPModeling")

  # Initialize mapping
  cokey_to_model_map <- data.frame(
    sim_cokey = simulated_cokeys,
    gp_model_group = NA_character_,
    stringsAsFactors = FALSE
  )

  if (matching_strategy == "exact_cokey") {
    # Get component information for matching with validation
    available_cols <- intersect(
      c("cokey", "compname", "taxclname", "taxpartsize", "taxgrtgroup", "taxsuborder", "taxorder"),
      names(nrcs_combined_data)
    )

    if (length(available_cols) <= 1) {
      log_message("WARN", "Insufficient taxonomic columns for matching. Using default model.", category = "GPModeling")
      cokey_to_model_map$gp_model_group <- available_groups[1]
      return(cokey_to_model_map)
    }

    # Create component lookup
    component_lookup <- nrcs_combined_data[, available_cols, drop = FALSE] |>
      dplyr::distinct()

    # Apply hierarchical matching
    cokey_to_model_map <- cokey_to_model_map |>
      dplyr::left_join(
        component_lookup |> dplyr::rename(sim_cokey = cokey),
        by = "sim_cokey"
      ) |>
      dplyr::mutate(gp_model_group = NA_character_)

    # Apply matching hierarchy step by step
    cokey_to_model_map <- apply_matching_hierarchy(cokey_to_model_map, available_groups)
  }

  # Report matching statistics
  log_message("INFO", "=== MATCHING SUMMARY ===", category = "GPModeling")
  matching_summary <- cokey_to_model_map |>
    dplyr::count(gp_model_group, sort = TRUE)

  for (i in 1:nrow(matching_summary)) {
    log_message("INFO", paste("Group", matching_summary$gp_model_group[i], ":", matching_summary$n[i], "soils"), category = "GPModeling")
  }

  unmatched <- sum(is.na(cokey_to_model_map$gp_model_group))
  if (unmatched > 0) {
    log_message("WARN", paste(unmatched, "simulated soils could not be matched"), category = "GPModeling")
  }

  return(cokey_to_model_map |> dplyr::select(sim_cokey, gp_model_group))
}

#' Predict GP Depth Trends
#'
#'
#' @param gp_model_info GP model information from fit_individual_gp_model()
#' @param new_depths Vector of depths for prediction
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Vector of predictions
#' @family gp-modeling
#' @export
predict_gp_depth_trends <- function(gp_model_info, new_depths,
                                    verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(gp_model_info) || is.null(gp_model_info$gp_model)) {
    log_message("WARN", "Invalid GP model provided for prediction", category = "GPModeling")
    return(rep(NA, length(new_depths)))
  }

  tryCatch({
    # Scale new depths using training scaling
    scaling <- gp_model_info$depth_scaling
    scaled_depths <- (new_depths - scaling$min) / scaling$range
    scaled_depths <- pmax(0, pmin(1, scaled_depths))  # Clamp to [0,1]

    # Predict
    predictions <- GPfit::predict.GP(
      gp_model_info$gp_model,
      xnew = as.matrix(scaled_depths)
    )

    return(predictions$Y_hat)

  }, error = function(e) {
    handle_workflow_error(e, "GP depth prediction", "warn")
    return(rep(NA, length(new_depths)))
  })
}

# ============================================================================
# 6. MODEL VALIDATION AND DIAGNOSTICS# ============================================================================

#' Validate GP Models
#'
#'
#' @param gp_models List of GP models from fit_depth_gp_models()
#' @param nrcs_data Original NRCS training data
#' @param validation_depths Depths for prediction testing
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results
#' @family gp-modeling
#' @export
validate_gp_models <- function(gp_models, nrcs_data, validation_depths = seq(0, 200, by = 10),
                               verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== GP MODEL VALIDATION ===", category = "GPModeling")

  validation_results <- list(
    model_validation = list(),
    prediction_validation = list(),
    overall_assessment = list(),
    recommendations = character(0)
  )

  properties <- names(gp_models)[names(gp_models) != "model_summary"]

  for (prop in properties) {
    log_message("INFO", paste("Validating", prop, "models..."), category = "GPModeling")

    model_info <- gp_models[[prop]]
    property_validation <- list()

    if (model_info$type == "stratified_grouped") {
      # Validate each group model
      group_validations <- list()

      for (group in names(model_info$models)) {
        group_model <- model_info$models[[group]]

        if (!is.null(group_model)) {
          # Test predictions at standard depths
          predictions <- tryCatch({
            predict_gp_depth_trends(group_model, validation_depths)
          }, error = function(e) {
            handle_workflow_error(e, paste("Validation prediction for", prop, group), "warn")
            rep(NA, length(validation_depths))
          })

          # Calculate validation metrics
          group_validation <- list(
            depths = validation_depths,
            predictions = predictions,
            model_valid = !any(is.na(predictions)),
            prediction_range = range(predictions, na.rm = TRUE),
            trend_monotonicity = assess_trend_monotonicity(predictions),
            realistic_values = assess_realistic_values(predictions, prop)
          )

          group_validations[[group]] <- group_validation
        }
      }

      property_validation <- group_validations

    } else {
      # Validate single model
      if (!is.null(model_info$model)) {
        predictions <- tryCatch({
          predict_gp_depth_trends(model_info$model, validation_depths)
        }, error = function(e) {
          handle_workflow_error(e, paste("Validation prediction for", prop), "warn")
          rep(NA, length(validation_depths))
        })

        property_validation <- list(
          depths = validation_depths,
          predictions = predictions,
          model_valid = !any(is.na(predictions)),
          prediction_range = range(predictions, na.rm = TRUE),
          trend_monotonicity = assess_trend_monotonicity(predictions),
          realistic_values = assess_realistic_values(predictions, prop)
        )
      }
    }

    validation_results$model_validation[[prop]] <- property_validation

    # Assess property-specific issues
    if (is.list(property_validation) && length(property_validation) > 0) {
      # Count valid models
      valid_models <- sum(sapply(property_validation, function(x) x$model_valid))
      total_models <- length(property_validation)

      log_message("INFO", paste(prop, ":", valid_models, "/", total_models, "models valid"), category = "GPModeling")

      if (valid_models < total_models) {
        validation_results$recommendations <- c(
          validation_results$recommendations,
          paste("Some", prop, "models failed validation - check training data quality")
        )
      }

      # Check for unrealistic predictions
      unrealistic_models <- sum(sapply(property_validation, function(x) {
        !x$realistic_values
      }))

      if (unrealistic_models > 0) {
        validation_results$recommendations <- c(
          validation_results$recommendations,
          paste(unrealistic_models, prop, "models produce unrealistic predictions")
        )
      }
    }
  }

  # Overall assessment
  total_models <- sum(sapply(gp_models[properties], function(x) x$n_groups))
  valid_models <- sum(sapply(validation_results$model_validation, function(prop_val) {
    sum(sapply(prop_val, function(x) x$model_valid))
  }))

  validation_results$overall_assessment <- list(
    total_models = total_models,
    valid_models = valid_models,
    success_rate = valid_models / total_models,
    validation_passed = valid_models / total_models >= 0.8
  )

  log_message("INFO", paste("Overall validation - Total:", total_models, "Valid:", valid_models,
                            "Success rate:", round(validation_results$overall_assessment$success_rate * 100, 1), "%"), category = "GPModeling")

  # Log recommendations
  if (length(validation_results$recommendations) > 0) {
    log_message("WARN", "Validation recommendations:", category = "GPModeling")
    for (rec in validation_results$recommendations) {
      log_message("WARN", paste("-", rec), category = "GPModeling")
    }
  } else {
    log_message("INFO", "All models passed validation", category = "GPModeling")
  }

  return(validation_results)
}

# ============================================================================
# 7. SUPPORTING HELPER FUNCTIONS
# ============================================================================

# Create soil groups based on strategy
create_soil_groups <- function(data, strategy) {
  available_cols <- names(data)

  if (strategy == "soil_series" && "compname" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$compname), data$compname, "unknown")
  } else if (strategy == "taxonomic_class" && "taxclname" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$taxclname), data$taxclname, "unknown")
  } else if (strategy == "particle_size" && "taxpartsize" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$taxpartsize), data$taxpartsize, "unknown")
  } else if (strategy == "soil_grtgroup" && "taxgrtgroup" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$taxgrtgroup), data$taxgrtgroup, "unknown")
  } else if (strategy == "soil_suborder" && "taxsuborder" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$taxsuborder), data$taxsuborder, "unknown")
  } else if (strategy == "soil_order" && "taxorder" %in% available_cols) {
    data$soil_group <- ifelse(!is.na(data$taxorder), data$taxorder, "unknown")
  } else if (strategy == "none") {
    data$soil_group <- "all_soils"
  } else {
    # Fallback hierarchy
    if ("compname" %in% available_cols) {
      data$soil_group <- ifelse(!is.na(data$compname), data$compname, "unknown")
    } else if ("taxclname" %in% available_cols) {
      data$soil_group <- ifelse(!is.na(data$taxclname), data$taxclname, "unknown")
    } else {
      data$soil_group <- "unknown"
    }
  }

  return(data)
}

# Generate processing summary
generate_gp_processing_summary <- function(processed_data, grouping_strategy) {
  log_message("INFO", "=== TRAINING DATA PROCESSING SUMMARY ===", category = "GPModeling")
  log_message("INFO", paste("Total observations:", nrow(processed_data)), category = "GPModeling")
  log_message("INFO", paste("Unique soil profiles:", dplyr::n_distinct(processed_data$cokey)), category = "GPModeling")
  log_message("INFO", paste("Grouping strategy:", grouping_strategy), category = "GPModeling")
  log_message("INFO", paste("Number of soil groups:", dplyr::n_distinct(processed_data$soil_group)), category = "GPModeling")

  # Show group summary
  group_summary <- processed_data |>
    dplyr::group_by(soil_group) |>
    dplyr::summarise(
      n_profiles = dplyr::n_distinct(cokey),
      n_observations = dplyr::n(),
      depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_observations))

  log_message("DEBUG", "Top soil groups by observation count:", category = "GPModeling")
  for (i in 1:min(5, nrow(group_summary))) {
    log_message("DEBUG", paste("  ", group_summary$soil_group[i], ":", group_summary$n_observations[i], "obs"), category = "GPModeling")
  }

  # Property data availability
  property_cols <- c("clay_pct", "sand_pct", "silt_pct", "pH", "organic_matter", "bulk_density")
  available_props <- intersect(property_cols, names(processed_data))

  if (length(available_props) > 0) {
    log_message("DEBUG", "Property data availability:", category = "GPModeling")
    for (prop in available_props) {
      n_available <- sum(!is.na(processed_data[[prop]]))
      pct_available <- round(n_available / nrow(processed_data) * 100, 1)
      log_message("DEBUG", paste("  ", prop, ":", n_available, "(", pct_available, "%)"), category = "GPModeling")
    }
  }
}

# Create test grouping for strategy evaluation
create_test_grouping <- function(data, strategy) {
  if (strategy == "soil_series" && "compname" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(compname), compname, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "taxonomic_class" && "taxclname" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(taxclname), taxclname, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "particle_size" && "taxpartsize" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(taxpartsize), taxpartsize, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "soil_grtgroup" && "taxgrtgroup" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(taxgrtgroup), taxgrtgroup, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "soil_suborder" && "taxsuborder" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(taxsuborder), taxsuborder, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "soil_order" && "taxorder" %in% names(data)) {
    test_data <- data |>
      dplyr::mutate(test_group = ifelse(!is.na(taxorder), taxorder, NA_character_)) |>
      dplyr::filter(!is.na(test_group))
  } else if (strategy == "none") {
    test_data <- data |>
      dplyr::mutate(test_group = "all_soils")
  } else {
    test_data <- data.frame()
  }

  return(test_data)
}

# Apply hierarchical grouping with fallbacks
apply_hierarchical_grouping <- function(data, property, min_profiles, min_obs) {

  # Check original grouping adequacy
  original_summary <- data |>
    dplyr::filter(!is.na(.data[[property]])) |>
    dplyr::group_by(soil_group) |>
    dplyr::summarise(
      n_profiles = dplyr::n_distinct(cokey),
      n_observations = dplyr::n(),
      depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
      .groups = "drop"
    )

  # Identify adequate and inadequate groups
  adequate_groups <- original_summary |>
    dplyr::filter(
      n_profiles >= min_profiles,
      n_observations >= min_obs,
      depth_range >= 20
    ) |>
    dplyr::pull(soil_group)

  inadequate_groups <- original_summary |>
    dplyr::filter(!soil_group %in% adequate_groups) |>
    dplyr::pull(soil_group)

  log_message("DEBUG", paste("Adequate groups:", length(adequate_groups), "Small groups needing fallback:", length(inadequate_groups)), category = "GPModeling")

  # Apply hierarchical fallback
  data_with_fallback <- data |>
    dplyr::mutate(
      final_group = dplyr::if_else(
        soil_group %in% adequate_groups,
        soil_group,
        "needs_fallback"
      )
    )

  # Apply fallbacks step by step
  fallback_hierarchy <- list(
    list(col = "taxpartsize", prefix = "fallback_particle_"),
    list(col = "taxgrtgroup", prefix = "fallback_grtgroup_"),
    list(col = "taxsuborder", prefix = "fallback_suborder_"),
    list(col = "taxorder", prefix = "fallback_order_")
  )

  for (fallback in fallback_hierarchy) {
    col_name <- fallback$col
    prefix <- fallback$prefix

    if (col_name %in% names(data_with_fallback)) {
      data_with_fallback <- data_with_fallback |>
        dplyr::mutate(
          final_group = dplyr::if_else(
            final_group == "needs_fallback" & !is.na(.data[[col_name]]),
            paste0(prefix, .data[[col_name]]),
            final_group
          )
        )
    }
  }

  # Final fallback
  data_with_fallback <- data_with_fallback |>
    dplyr::mutate(
      final_group = dplyr::if_else(
        final_group == "needs_fallback",
        "fallback_general",
        final_group
      )
    )

  # Check if fallback groups are adequate and merge small ones
  fallback_summary <- data_with_fallback |>
    dplyr::filter(!is.na(.data[[property]])) |>
    dplyr::group_by(final_group) |>
    dplyr::summarise(
      n_profiles = dplyr::n_distinct(cokey),
      n_observations = dplyr::n(),
      depth_range = max(hzdept_r, na.rm = TRUE) - min(hzdept_r, na.rm = TRUE),
      .groups = "drop"
    )

  still_inadequate <- fallback_summary |>
    dplyr::filter(
      n_profiles < min_profiles |
        n_observations < min_obs |
        depth_range < 20
    ) |>
    dplyr::pull(final_group)

  if (length(still_inadequate) > 0) {
    log_message("DEBUG", paste("Merging", length(still_inadequate), "remaining small groups into general pool"), category = "GPModeling")

    data_with_fallback <- data_with_fallback |>
      dplyr::mutate(
        final_group = dplyr::if_else(
          final_group %in% still_inadequate,
          "general_pool",
          final_group
        )
      )
  }

  return(data_with_fallback)
}

# Apply matching hierarchy for soil-to-model matching
apply_matching_hierarchy <- function(cokey_to_model_map, available_groups) {

  # Direct matches first
  if ("compname" %in% names(cokey_to_model_map)) {
    cokey_to_model_map <- cokey_to_model_map |>
      dplyr::mutate(
        gp_model_group = dplyr::if_else(
          is.na(gp_model_group) & compname %in% available_groups,
          compname,
          gp_model_group
        )
      )
  }

  if ("taxclname" %in% names(cokey_to_model_map)) {
    cokey_to_model_map <- cokey_to_model_map |>
      dplyr::mutate(
        gp_model_group = dplyr::if_else(
          is.na(gp_model_group) & taxclname %in% available_groups,
          taxclname,
          gp_model_group
        )
      )
  }

  # Fallback group matching
  fallback_mappings <- list(
    list(col = "taxpartsize", prefix = "fallback_particle_"),
    list(col = "taxgrtgroup", prefix = "fallback_grtgroup_"),
    list(col = "taxsuborder", prefix = "fallback_suborder_"),
    list(col = "taxorder", prefix = "fallback_order_")
  )

  for (mapping in fallback_mappings) {
    col_name <- mapping$col
    prefix <- mapping$prefix

    if (col_name %in% names(cokey_to_model_map)) {
      cokey_to_model_map <- cokey_to_model_map |>
        dplyr::mutate(
          fallback_group = paste0(prefix, .data[[col_name]]),
          gp_model_group = dplyr::if_else(
            is.na(gp_model_group) & fallback_group %in% available_groups,
            fallback_group,
            gp_model_group
          )
        ) |>
        dplyr::select(-fallback_group)
    }
  }

  # Final fallbacks
  if ("fallback_general" %in% available_groups) {
    cokey_to_model_map <- cokey_to_model_map |>
      dplyr::mutate(
        gp_model_group = dplyr::if_else(
          is.na(gp_model_group),
          "fallback_general",
          gp_model_group
        )
      )
  }

  if ("general_pool" %in% available_groups) {
    cokey_to_model_map <- cokey_to_model_map |>
      dplyr::mutate(
        gp_model_group = dplyr::if_else(
          is.na(gp_model_group),
          "general_pool",
          gp_model_group
        )
      )
  }

  # Last resort
  cokey_to_model_map <- cokey_to_model_map |>
    dplyr::mutate(
      gp_model_group = dplyr::if_else(
        is.na(gp_model_group),
        available_groups[1],
        gp_model_group
      )
    )

  return(cokey_to_model_map)
}

# Create model summary
create_model_summary <- function(gp_models, processed_nrcs_data) {
  properties_modeled <- names(gp_models)[names(gp_models) != "model_summary"]

  list(
    properties_modeled = properties_modeled,
    total_models = sum(sapply(gp_models[properties_modeled], function(x) x$n_groups)),
    timestamp = Sys.time(),
    training_data_summary = list(
      n_observations = nrow(processed_nrcs_data),
      n_profiles = dplyr::n_distinct(processed_nrcs_data$cokey),
      depth_range = c(min(processed_nrcs_data$hzdept_r, na.rm = TRUE),
                      max(processed_nrcs_data$hzdept_r, na.rm = TRUE))
    )
  )
}

#' K-Fold Cross-Validation for GP Correlation Structure Selection
#'
#' Runs `n_folds`-fold cross-validation of [GPfit::GP_fit()] for each supplied
#' correlation-family candidate, returning the mean held-out RMSE per
#' candidate. Shared helper used by [optimize_gp_hyperparameters()] and
#' `perform_gp_cross_validation()` (diagnostics.R) so the
#' fold-splitting logic is not duplicated across files.
#'
#' @param X Predictor matrix (one row per observation).
#' @param Y Response vector.
#' @param n_folds Number of cross-validation folds.
#' @param corr_candidates Named list of `corr` specs to pass to
#'   [GPfit::GP_fit()] (e.g. `list(type = "exponential", power = 1.95)`).
#'
#' @param gp_control Passed through to each fold's `GPfit::GP_fit()` `control` argument - see
#'   `fit_individual_gp_model()`'s docs for why the default is much smaller than `GP_fit()`'s own
#'   default. This is the highest-multiplier call site in the package: `n_folds x
#'   length(corr_candidates)` `GP_fit()` calls per invocation.
#' @return List with `mean_rmse_by_candidate`, a named list of mean held-out
#'   RMSE per candidate (`NA` for a candidate that failed on every fold).
k_fold_gp_cv <- function(X, Y, n_folds, corr_candidates, gp_control = c(20, 10, 2)) {
  X <- as.matrix(X)
  n <- length(Y)
  n_folds <- max(2, min(n_folds, n))
  fold_id <- sample(rep(seq_len(n_folds), length.out = n))

  mean_rmse_by_candidate <- list()

  for (cand_name in names(corr_candidates)) {
    corr_spec <- corr_candidates[[cand_name]]
    fold_rmses <- numeric(0)

    for (k in seq_len(n_folds)) {
      test_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k)

      if (length(train_idx) < 2 || length(test_idx) < 1) next

      fold_rmse <- tryCatch({
        model <- GPfit::GP_fit(X = X[train_idx, , drop = FALSE], Y = Y[train_idx], corr = corr_spec,
                                control = gp_control)
        pred <- GPfit::predict.GP(model, xnew = X[test_idx, , drop = FALSE])
        sqrt(mean((pred$Y_hat - Y[test_idx])^2))
      }, error = function(e) NA_real_)

      if (is.finite(fold_rmse)) {
        fold_rmses <- c(fold_rmses, fold_rmse)
      }
    }

    mean_rmse_by_candidate[[cand_name]] <- if (length(fold_rmses) > 0) mean(fold_rmses) else NA_real_
  }

  list(mean_rmse_by_candidate = mean_rmse_by_candidate, folds = n_folds)
}

#' Optimize GP Hyperparameters
#'
#' Performs real `n_folds`-fold cross-validation comparing a small set of
#' [GPfit::GP_fit()] correlation-family specs (exponential power=1.95,
#' Matern nu=3/2, Matern nu=5/2 - the practically tunable surface GPfit
#' exposes), picks the lowest mean cross-validated RMSE, and refits on the
#' full data with the winning spec. CV results (folds requested, candidates
#' compared, winning candidate, per-candidate mean RMSE) are attached as the
#' `cv_results` attribute on the returned model so the search performed is
#' inspectable without changing the model's class/shape.
#'
#' @param X Scaled predictor matrix (depths).
#' @param Y Response vector.
#' @param n_folds Number of cross-validation folds.
#' @param gp_control Passed through to every `GPfit::GP_fit()` call this function makes
#'   (`k_fold_gp_cv()`'s per-fold fits, the winning-candidate refit, and the fallback/baseline
#'   fits) - see `fit_individual_gp_model()`'s docs for why the default is much smaller than
#'   `GP_fit()`'s own default. With the default 5 folds x 3 correlation candidates + a refit,
#'   this is ~16 `GP_fit()` calls per invocation.
#'
#' @return A `"GP"`-classed [GPfit::GP_fit()] model, with a `cv_results`
#'   attribute describing the cross-validation that selected it.
optimize_gp_hyperparameters <- function(X, Y, n_folds = 5, gp_control = c(20, 10, 2)) {
  tryCatch({
    X <- as.matrix(X)
    n <- length(Y)

    corr_candidates <- list(
      exponential_power1.95 = list(type = "exponential", power = 1.95),
      matern_nu1.5 = list(type = "matern", nu = 3 / 2),
      matern_nu2.5 = list(type = "matern", nu = 5 / 2)
    )

    effective_folds <- max(2, min(n_folds, n))

    # Cross-validation needs at least 2 training points per fold to be meaningful; fall back
    # to a single baseline model otherwise.
    cv <- if (n >= 2 * effective_folds) {
      k_fold_gp_cv(X, Y, effective_folds, corr_candidates, gp_control = gp_control)
    } else {
      NULL
    }

    rmses <- if (!is.null(cv)) unlist(cv$mean_rmse_by_candidate) else numeric(0)

    if (length(rmses) > 0 && any(is.finite(rmses))) {
      best_name <- names(which.min(rmses))
      best_corr <- corr_candidates[[best_name]]

      final_model <- tryCatch({
        GPfit::GP_fit(X = X, Y = Y, corr = best_corr, control = gp_control)
      }, error = function(e) {
        GPfit::GP_fit(X = X, Y = Y, control = gp_control)
      })

      attr(final_model, "cv_results") <- list(
        folds = effective_folds,
        corr_candidates = names(corr_candidates),
        best_candidate = best_name,
        mean_rmse_by_candidate = cv$mean_rmse_by_candidate
      )

      return(final_model)
    }

    # Not enough data for meaningful cross-validation - single baseline fit.
    baseline_model <- GPfit::GP_fit(X = X, Y = Y, control = gp_control)
    attr(baseline_model, "cv_results") <- list(
      folds = effective_folds,
      corr_candidates = names(corr_candidates),
      best_candidate = NA_character_,
      mean_rmse_by_candidate = list()
    )
    return(baseline_model)

  }, error = function(e) {
    log_message("WARN", paste("Hyperparameter optimization failed, using default:", e$message), category = "GPModeling")
    return(GPfit::GP_fit(X = as.matrix(X), Y = Y, control = gp_control))
  })
}

# Calculate model diagnostics
calculate_model_diagnostics <- function(gp_model, training_data, property) {
  diagnostics <- list(
    property = property,
    n_training_points = nrow(training_data),
    depth_range = range(training_data$hzdept_r),
    value_range = range(training_data$mean_value, na.rm = TRUE),
    training_rmse = NA,
    log_likelihood = NA
  )

  tryCatch({
    # Calculate training RMSE
    scaled_depths <- (training_data$hzdept_r - min(training_data$hzdept_r)) /
      (max(training_data$hzdept_r) - min(training_data$hzdept_r))

    predictions <- GPfit::predict.GP(gp_model, xnew = as.matrix(scaled_depths))
    diagnostics$training_rmse <- sqrt(mean((predictions$Y_hat - training_data$mean_value)^2))

    # Get log likelihood if available
    if ("loglikelihood" %in% names(gp_model)) {
      diagnostics$log_likelihood <- gp_model$loglikelihood
    }

  }, error = function(e) {
    log_message("WARN", paste("Error calculating diagnostics for", property, ":", e$message), category = "GPModeling")
  })

  return(diagnostics)
}

# Trend assessment functions (using shared safe operations)
assess_trend_monotonicity <- function(predictions) {
  if (length(predictions) < 3 || all(is.na(predictions))) {
    return(NA)
  }

  # Check if trend is generally monotonic (allowing some variation)
  diffs <- diff(predictions)
  pos_diffs <- sum(diffs > 0, na.rm = TRUE)
  neg_diffs <- sum(diffs < 0, na.rm = TRUE)

  # Consider monotonic if 70% of changes are in the same direction
  monotonicity_score <- max(pos_diffs, neg_diffs) / length(diffs)

  return(monotonicity_score >= 0.7)
}

assess_realistic_values <- function(predictions, property) {
  if (all(is.na(predictions))) {
    return(FALSE)
  }

  # Property-specific realistic ranges
  realistic_ranges <- list(
    clay_pct = c(0, 100),
    sand_pct = c(0, 100),
    silt_pct = c(0, 100),
    pH = c(2.5, 11.0),
    organic_matter = c(0, 50),
    bulk_density = c(0.3, 3.0),
    cec = c(0, 200),
    awc = c(0, 0.5)
  )

  if (property %in% names(realistic_ranges)) {
    range_vals <- realistic_ranges[[property]]
    return(all(predictions >= range_vals[1] & predictions <= range_vals[2], na.rm = TRUE))
  }

  # Generic check - no negative values, reasonable magnitude
  return(all(predictions >= 0, na.rm = TRUE) && all(predictions < 1000, na.rm = TRUE))
}

#' Infer Simulation Property Columns
#'
#' Infers which columns of a long-format simulation data frame represent
#' soil properties to adjust, as the numeric columns remaining after
#' excluding known structural/metadata columns.
#'
#' @param simulation_data Long-format simulation data frame.
#'
#' @return Character vector of inferred property column names.
infer_simulation_properties <- function(simulation_data) {
  structural_cols <- c("cokey", "hzdept_r", "hzdepb_r", "hzname", "genhz", "mukey",
                       "simulation_number", "infill_method", "unsuitable_horizon",
                       "property_data_complete")
  numeric_cols <- names(simulation_data)[vapply(simulation_data, is.numeric, logical(1))]
  setdiff(numeric_cols, structural_cols)
}

#' Apply NRCS GP Adjustments With Correlation Preservation
#'
#' Delegates to the real `apply_nrcs_trend_adjustments()`
#' (`core-gp.R`), which already implements NRCS GP-model
#' depth-trend adjustment with correlation preservation. This wrapper's only
#' job is to infer the `properties` argument that function requires (as the
#' numeric, non-structural columns of `simulation_data`), since this
#' function's own signature does not accept one explicitly.
#'
#' @param simulation_data Simulated soil property data for one cokey.
#' @param nrcs_gp_models Fitted NRCS GP models (see [fit_depth_gp_models()]).
#' @param model_group GP model group to apply.
#' @param preserve_correlations Whether to preserve within-depth correlations.
#'
#' @return Adjusted simulation data.
apply_nrcs_gp_adjustments_with_correlations <- function(simulation_data, nrcs_gp_models, model_group, preserve_correlations = TRUE) {
  log_message("INFO", "Applying NRCS GP adjustments with correlation preservation", category = "GPModeling")

  properties <- infer_simulation_properties(simulation_data)
  if (length(properties) == 0) {
    log_message("WARN", "No numeric property columns found for NRCS GP adjustment", category = "GPModeling")
    return(simulation_data)
  }

  apply_nrcs_trend_adjustments(
    cokey_data = simulation_data,
    gp_models = nrcs_gp_models,
    model_group = model_group,
    properties = properties,
    preserve_correlations = preserve_correlations
  )
}

#' Apply Local GP Adjustments With Correlation Preservation
#'
#' Delegates to [apply_local_gp_adjustments()], inferring the `properties` argument that
#' function requires (as the numeric, non-structural columns of `simulation_data`),
#' since this function's own signature does not accept one explicitly.
#'
#' @param simulation_data Simulated soil property data for one cokey.
#' @param preserve_correlations Whether to preserve within-depth correlations.
#'
#' @return Adjusted simulation data.
apply_local_gp_adjustments_with_correlations <- function(simulation_data, preserve_correlations = TRUE) {
  log_message("INFO", "Applying local GP adjustments with correlation preservation", category = "GPModeling")

  properties <- infer_simulation_properties(simulation_data)
  if (length(properties) == 0) {
    log_message("WARN", "No numeric property columns found for local GP adjustment", category = "GPModeling")
    return(simulation_data)
  }

  apply_local_gp_adjustments(
    cokey_data = simulation_data,
    properties = properties,
    preserve_correlations = preserve_correlations
  )
}

validate_enhanced_simulations <- function(original_data, enhanced_data, preserve_correlations) {
  # Basic validation structure
  validation <- list(
    n_rows_original = nrow(original_data),
    n_rows_enhanced = nrow(enhanced_data),
    dimensions_match = nrow(original_data) == nrow(enhanced_data)
  )

  log_message("INFO", paste("Validation - Dimensions match:", validation$dimensions_match), category = "GPModeling")

  return(validation)
}

# Export GP models function
export_gp_models <- function(gp_models, output_path, include_training_data = FALSE) {
  log_message("INFO", paste("Exporting GP models to:", output_path), category = "GPModeling")

  tryCatch({
    export_data <- list(
      models = gp_models,
      export_metadata = list(
        timestamp = Sys.time(),
        r_version = R.version.string,
        include_training_data = include_training_data
      )
    )

    # Optionally remove training data to reduce file size
    if (!include_training_data) {
      properties <- names(gp_models)[names(gp_models) != "model_summary"]
      for (prop in properties) {
        if (gp_models[[prop]]$type == "stratified_grouped") {
          for (group in names(gp_models[[prop]]$models)) {
            if (!is.null(gp_models[[prop]]$models[[group]]$training_data)) {
              export_data$models[[prop]]$models[[group]]$training_data <- NULL
            }
          }
        }
      }
      log_message("DEBUG", "Training data excluded from export", category = "GPModeling")
    }

    # Save to file
    saveRDS(export_data, file = output_path)

    if (file.exists(output_path)) {
      file_size <- file.size(output_path)
      log_message("INFO", paste("Models exported successfully, file size:", round(file_size / 1024^2, 2), "MB"), category = "GPModeling")
      return(TRUE)
    } else {
      stop("Export file was not created")
    }

  }, error = function(e) {
    handle_workflow_error(e, "GP model export", "warn")
    return(FALSE)
  })
}

# Load exported models function
load_gp_models <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("GP model file not found: ", input_path)
  }

  tryCatch({
    loaded_data <- readRDS(input_path)

    log_message("INFO", "=== LOADING GP MODELS ===", category = "GPModeling")
    log_message("INFO", paste("Export timestamp:", as.character(loaded_data$export_metadata$timestamp)), category = "GPModeling")

    if (!is.null(loaded_data$models$model_summary)) {
      log_message("INFO", paste("Properties modeled:", length(loaded_data$models$model_summary$properties_modeled)), category = "GPModeling")
      log_message("INFO", paste("Total models:", loaded_data$models$model_summary$total_models), category = "GPModeling")
    }

    return(loaded_data$models)

  }, error = function(e) {
    stop("Failed to load GP models: ", e$message)
  })
}


# ============================================================================
# merged from core-gp.R (P1 file reorg)
# ============================================================================

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
#'
#' @param simulation_results Results from simulate_monte_carlo()
#' @param gp_models Optional fitted NRCS GP depth models (from fit_depth_gp_models())
#' @param cokey_mapping Optional mapping from simulations to NRCS GP groups
#' @param integration_method Method: "nrcs_gp", "local_gp", or "hybrid" (default = "hybrid")
#' @param preserve_correlations Whether to preserve within-depth correlations (default = TRUE)
#' @param properties Properties to adjust (NULL = auto-detect)
#' @param parallel Whether to use parallel processing (default = FALSE)
#' @param n_cores Number of cores for parallel processing
#' @param config Integration configuration (uses the package default configuration if NULL)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Integrated simulation results with realistic depth trends
#' @keywords internal
apply_depth_gp_to_simulation <- function(simulation_results,
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

  # Load configuration
  if (is.null(config)) {
    config <- default_config("full")
  }

  # Validate parameters
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

  # Input validation
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

  # Validate properties
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

  # Split the simulation data by cokey once here and pass each component its own subset down,
  # rather than filtering the full multi-cokey table once per cokey inside process_single_cokey()
  # (which would be O(rows x cokeys)). Same pattern as
  # maybe_adjust_soil_data_depth_trend()/run_fusion_group() elsewhere in the package.
  cokey_groups <- split(simulation_data, simulation_data$cokey)

  # Process cokeys with progress tracking
  if (parallel && length(unique_cokeys) > 1) {
    integrated_results <- process_cokeys_parallel(
      cokey_groups, unique_cokeys, properties, gp_models, cokey_mapping,
      use_nrcs_gp, use_local_gp, preserve_correlations, n_cores, config
    )
  } else {
    integrated_results <- process_cokeys_sequential(
      cokey_groups, unique_cokeys, properties, gp_models, cokey_mapping,
      use_nrcs_gp, use_local_gp, preserve_correlations, config
    )
  }

  # Combine and validate results
  final_data <- combine_and_validate_results(integrated_results, unique_cokeys, simulation_data)

  # Comprehensive validation
  log_message("INFO", "=== INTEGRATION VALIDATION ===", category = "MultivarAdjust")
  validation_results <- validate_integration_results(
    simulation_data, final_data, properties, preserve_correlations
  )

  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "secs")

  # Prepare final output 
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
#' Applies GP-predicted depth trends to one component's simulation data, optionally
#' preserving within-depth cross-property correlations.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param gp_predictions Named list of GP predictions by property
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param primary_property Reference property for correlation preservation
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @param config Optional Monte Carlo config (as from \code{default_monte_carlo_config()}) whose
#'   \code{monte_carlo$vertical_correlation_method} selects between \code{"joint_copula"}
#'   (the default; dispatches to
#'   \code{preserve_correlation_structure_joint()}, drawing depth correlation and property
#'   correlation simultaneously) and \code{"gp_quantile_retrofit"} (an explicit opt-out that dispatches to
#'   \code{preserve_correlation_structure()}). \code{NULL} (default) resolves to
#'   \code{"joint_copula"}, matching \code{default_monte_carlo_config()}'s own default - set
#'   \code{config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit"} explicitly to
#'   select the quantile-retrofit method. Under \code{"joint_copula"},
#'   \code{config$monte_carlo$vertical_correlation_gating} (default \code{FALSE}) separately
#'   controls whether \code{bound_sd}-based discontinuity gating is applied. Its numeric defaults are conservative and kept independent of the core
#'   method choice.
#' @param gp_models Optional named list of fitted GP models (as `fit_local_gp_model_single()`
#'   returns), keyed by property - passed through to \code{preserve_correlation_structure_joint()}
#'   when \code{vertical_correlation_method = "joint_copula"}, so its depth kernel can reuse each
#'   property's already-fitted, already-cross-validated length-scale (see
#'   \code{extract_depth_length_scale()}) instead of falling back to a full-depth-range default.
#'   Ignored entirely under \code{"gp_quantile_retrofit"}.
#' @return Adjusted simulation data
#' @family gp-modeling
#' @export
apply_gp_depth_trends <- function(cokey_data,
                                  gp_predictions,
                                  properties,
                                  preserve_correlations = TRUE,
                                  primary_property = NULL,
                                  verbose = getOption("ssurgo.verbose", FALSE),
                                  config = NULL,
                                  gp_models = NULL) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # validation
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

  # Vertical-correlation method. Default "joint_copula"; kept in sync with
  # default_monte_carlo_config() so "no config passed" means the same thing whether or not a caller
  # goes through default_monte_carlo_config() first. Set
  # config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit" for the
  # quantile-retrofit method instead.
  vertical_correlation_method <- config$monte_carlo$vertical_correlation_method %||% "joint_copula"

  # Discontinuity gating (build_depth_correlation_kernel()'s boundary_distinctness suppression) is
  # a separate opt-in from the joint_copula method itself: bound_sd is attached unconditionally
  # upstream (attach_osd_boundary_distinctness() in simulate_ssurgo_mapunit_draws()), so this flag
  # is what allows joint_copula to run without gating when OSD lookup succeeds. Its numeric
  # defaults are conservative; defaults to FALSE so method and gating are chosen independently.
  vertical_correlation_gating <- isTRUE(config$monte_carlo$vertical_correlation_gating)

  # Apply multivariate adjustment with error handling
  adjusted_matrices <- tryCatch({
    if (preserve_correlations && length(property_matrices) >= 2) {
      if (identical(vertical_correlation_method, "joint_copula")) {
        # bound_sd (OSD boundary distinctness) is broadcast identically across every simulation realization at a given depth
        # (simulate_cokey_generalized()), so the first non-NA value per unique depth is that
        # depth's boundary_distinctness for build_depth_correlation_kernel()'s discontinuity
        # gating - NULL (no gating) when the column isn't present, OR when
        # vertical_correlation_gating is not explicitly enabled.
        boundary_distinctness <- if (vertical_correlation_gating && "bound_sd" %in% names(valid_data)) {
          vapply(unique_depths, function(d) {
            vals <- valid_data$bound_sd[valid_data$hzdept_r == d]
            vals <- vals[!is.na(vals)]
            if (length(vals) > 0) vals[1] else NA_real_
          }, numeric(1))
        } else {
          NULL
        }

        preserve_correlation_structure_joint(
          property_matrices, gp_predictions, unique_depths, primary_property,
          gp_models = gp_models, boundary_distinctness = boundary_distinctness
        )
      } else {
        preserve_correlation_structure(
          property_matrices, gp_predictions, unique_depths, primary_property
        )
      }
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

  # Merge with original data
  result_data <- merge_adjusted_data(cokey_data, adjusted_data, available_properties)

  return(result_data)
}

# ============================================================================
# 1b. JOINT DEPTH x PROPERTY COPULA
# ============================================================================
#
# `preserve_correlation_structure()` below imposes vertical correlation on
# independent-across-depth draws via a sequential gp_ratio nudge keyed off one "primary
# property"'s rank. The two functions here instead draw depth correlation (R_depth, from
# build_depth_correlation_kernel()) and property correlation (R_prop, the same flat matrix
# simulate_correlated_triangular() uses) simultaneously from a single Kronecker-separable joint
# distribution, so both hold by construction rather than by approximation.
# preserve_correlation_structure_joint() wires these two functions into a drop-in alternative
# with the same signature and contract as preserve_correlation_structure().

#' Draw a Joint Depth x Property Gaussian Copula Sample
#'
#' Draws `n_sims` realizations of an `n_depths x k` standard-normal field whose ROWS (depths) are
#' correlated according to `R_depth` and whose COLUMNS (properties) are correlated according to
#' `R_prop`, SIMULTANEOUSLY - the standard separable/Kronecker-structured multivariate-normal
#' sampling identity `Z = L_depth %*% Eps %*% t(L_prop)` (where `L_depth`/`L_prop` are Cholesky
#' factors of `R_depth`/`R_prop` and `Eps` is iid standard normal), which achieves
#' `Cov(vec(Z)) = R_prop (x) R_depth` (a Kronecker product) without ever materializing that full
#' `(n_depths*k) x (n_depths*k)` matrix. Vectorized across all `n_sims` realizations at once via
#' array reshaping (no explicit per-realization loop).
#'
#' @param R_depth An `n_depths x n_depths` correlation matrix (e.g. from
#'   `build_depth_correlation_kernel()`).
#' @param R_prop A `k x k` correlation matrix (the same flat property-correlation matrix
#'   `simulate_correlated_triangular()` already uses for within-horizon correlation).
#' @param n_sims Number of joint realizations to draw.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A `c(n_depths, k, n_sims)` array of standard-normal values (mean 0, variance 1
#'   marginally), jointly correlated across both the depth and property dimensions as specified.
#'
#' @keywords internal
sample_joint_depth_property_copula <- function(R_depth, R_prop, n_sims, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_depths <- nrow(R_depth)
  k <- nrow(R_prop)

  if (n_sims < 1) {
    stop("n_sims must be at least 1")
  }

  L_depth <- t(chol(R_depth))
  L_prop <- t(chol(R_prop))

  eps <- array(stats::rnorm(n_depths * k * n_sims), dim = c(n_depths, k, n_sims))

  # Left-multiply by L_depth (depth correlation) across every (property, realization) slice at
  # once: matrix(eps, nrow = n_depths) flattens the array's trailing dimensions column-major,
  # exactly matching R's array storage order, so a single n_depths x n_depths matrix multiply
  # applies to all k * n_sims columns simultaneously.
  step1 <- L_depth %*% matrix(eps, nrow = n_depths)
  dim(step1) <- c(n_depths, k, n_sims)

  # Right-multiply by t(L_prop) (property correlation) across every (depth, realization) slice at
  # once - reshape to bring the property dimension first so the same "one matrix multiply over a
  # flattened array" trick applies again, then reshape back.
  step1_prop_first <- aperm(step1, c(2, 1, 3))
  step2 <- L_prop %*% matrix(step1_prop_first, nrow = k)
  dim(step2) <- c(k, n_depths, n_sims)

  aperm(step2, c(2, 1, 3))
}

#' Map a Joint Copula Sample onto Existing Per-Depth Marginal Distributions
#'
#' Converts the standard-normal joint sample from `sample_joint_depth_property_copula()` into
#' actual property values, by probability-integral-transforming each depth/property/realization
#' cell (`pnorm()`) and then inverse-transforming (`quantile()`) against that depth/property's
#' OWN already-simulated marginal distribution (`property_matrices[[prop]][depth, ]` - the same
#' per-horizon triangular-fit values `simulate_correlated_triangular()` produced) - a Gaussian
#' copula, in the standard statistical sense. This is what actually delivers the achieved
#' correlation structure from `sample_joint_depth_property_copula()` onto real property values
#' while preserving each depth's own fitted marginal shape exactly (for large `n_sims`).
#'
#' @param Z A `c(n_depths, k, n_sims)` array as returned by
#'   `sample_joint_depth_property_copula()`.
#' @param property_matrices Named list of `k` matrices (rows = depths, columns = simulations,
#'   `dim(Z)[1]`/`dim(Z)[3]` must match), in the same property order as `Z`'s second dimension -
#'   the same shape `preserve_correlation_structure()` already consumes.
#' @param gp_predictions Optional named list of per-property depth-trend mean vectors (length
#'   `n_depths`), e.g. from `predict_gp_depth_trends()`. When supplied for a property, that
#'   property's marginal at each depth is re-centered (a location shift only - shape/spread
#'   preserved) onto the GP-predicted mean before the quantile mapping, replacing the old
#'   sequential `gp_ratio` nudge with a single direct shift.
#'
#' @return A named list of `n_depths x n_sims` matrices, same shape/names as `property_matrices`.
#'
#' @keywords internal
apply_copula_to_marginals <- function(Z, property_matrices, gp_predictions = NULL) {
  property_names <- names(property_matrices)

  if (length(property_names) == 0) {
    stop("property_matrices must have named elements")
  }
  if (dim(Z)[2] != length(property_names)) {
    stop("Z's property dimension (dim(Z)[2]) must match length(property_matrices)")
  }

  n_depths <- dim(Z)[1]
  n_sims <- dim(Z)[3]
  adjusted <- vector("list", length(property_names))
  names(adjusted) <- property_names

  for (p in seq_along(property_names)) {
    prop <- property_names[p]
    current_matrix <- property_matrices[[prop]]
    adjusted_matrix <- matrix(NA_real_, nrow = n_depths, ncol = n_sims)

    gp_means <- if (!is.null(gp_predictions) && !is.null(gp_predictions[[prop]]) &&
                     length(gp_predictions[[prop]]) == n_depths) {
      gp_predictions[[prop]]
    } else {
      NULL
    }

    for (i in seq_len(n_depths)) {
      curr_values <- current_matrix[i, ]

      # Optional GP-mean recentering: shift the marginal's LOCATION to the GP-predicted mean
      # while preserving its already-fitted shape/spread exactly - the direct replacement for
      # preserve_correlation_structure()'s sequential gp_ratio nudge.
      #
      # isTRUE()-wrapped: stats::var(x, na.rm = TRUE) returns NA (not FALSE) when x is entirely
      # NA (a real, non-rare case on messy field data - e.g. a cokey row whose texture triplet
      # was incomplete, see simulate_cokey_generalized()'s own per-row texture tryCatch()) -
      # `if (NA > 0)` errors with "missing value where TRUE/FALSE needed" rather than falling
      # through to the else branch. isTRUE() treats that NA as FALSE, matching this function's
      # intended "no variation (or no data) - keep original values" contract. (a real case on messy field data - e.g. a cokey whose texture triplet was incomplete).
      target_values <- if (!is.null(gp_means) && is.finite(gp_means[i]) &&
                            isTRUE(stats::var(curr_values, na.rm = TRUE) > 0)) {
        curr_values + (gp_means[i] - mean(curr_values, na.rm = TRUE))
      } else {
        curr_values
      }

      if (isTRUE(stats::var(target_values, na.rm = TRUE) > 0)) {
        u <- stats::pnorm(Z[i, p, ])
        adjusted_matrix[i, ] <- unname(stats::quantile(target_values, probs = u, na.rm = TRUE))
      } else {
        adjusted_matrix[i, ] <- target_values
      }
    }

    adjusted[[prop]] <- adjusted_matrix
  }

  adjusted
}

#' Preserve Correlation Structure During GP Adjustment
#'
#' Imposes vertical correlation on the property matrices by re-ordering each property's
#' per-depth values to a shared surface-referenced rank, so the depth trend is applied without
#' disturbing within-depth cross-property correlation.
#'
#' @param property_matrices Named list of property matrices
#' @param gp_predictions Named list of GP predictions
#' @param depths Depth vector
#' @param primary_property Reference property for correlation preservation
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of adjusted property matrices
#' @keywords internal
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

      # Get GP trend ratio
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

#' Preserve Correlation Structure via a Joint Depth x Property Copula
#'
#' Drop-in alternative to `preserve_correlation_structure()` - same required parameters, in the
#' same order, so existing call sites work unchanged - that replaces its sequential
#' `gp_ratio`/single-"primary-property" retrofit with the joint Kronecker-copula sampler from
#' the joint Kronecker-copula sampler (sample_joint_depth_property_copula() +
#' `apply_copula_to_marginals()`), so depth correlation and property correlation are satisfied
#' SIMULTANEOUSLY by construction rather than approximated by a rank-copying retrofit.
#' `primary_property` is accepted (for signature compatibility with
#' `preserve_correlation_structure()`) but unused - the joint method has no need for a single
#' reference property, since every property's correlation to every other is modeled directly via
#' `R_prop`.
#'
#' @param property_matrices Named list of property matrices (rows = depths, columns =
#'   simulations) - same shape `preserve_correlation_structure()` consumes.
#' @param gp_predictions Named list of GP-predicted per-depth mean vectors - passed through to
#'   `apply_copula_to_marginals()`'s optional GP-mean recentering.
#' @param depths Depth vector (real units, e.g. cm).
#' @param primary_property Accepted for signature compatibility with
#'   `preserve_correlation_structure()`; unused by the joint method.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @param gp_models Optional named list of fitted GP models (as returned by
#'   `fit_local_gp_model_single()`, or raw `GPfit`-classed objects), keyed by property. When
#'   supplied, the depth kernel's length-scale is derived from `extract_depth_length_scale()`
#'   applied to every property with a usable model, averaged across them (reusing each property's
#'   already-fitted, already-cross-validated GP fit). When `NULL`/empty (e.g. this function is called
#'   standalone, without the fitted models available), falls back to a length-scale spanning the
#'   full depth range (`diff(range(depths))`) - a conservative "moderate smooth correlation across
#'   the whole profile" default that keeps this function usable without requiring GP models.
#' @param boundary_distinctness Optional per-depth `bound_sd` vector, passed through to
#'   `build_depth_correlation_kernel()` for discontinuity gating. `NULL` (default)
#'   skips gating.
#' @param kernel `"exponential"` (default) or `"matern"` - passed through to
#'   `build_depth_correlation_kernel()`.
#'
#' @return List of adjusted property matrices - same shape/contract as
#'   `preserve_correlation_structure()`'s return value. Degrades gracefully to the
#'   unadjusted `property_matrices` (with a warning) on insufficient dimensions or a sampling/
#'   mapping failure, matching `preserve_correlation_structure()`'s own graceful-failure contract.
#' @keywords internal
preserve_correlation_structure_joint <- function(property_matrices,
                                                 gp_predictions,
                                                 depths,
                                                 primary_property,
                                                 verbose = getOption("ssurgo.verbose", FALSE),
                                                 gp_models = NULL,
                                                 boundary_distinctness = NULL,
                                                 kernel = c("exponential", "matern")) {
  kernel <- match.arg(kernel)

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Preserving correlation structure via joint depth x property copula",
              category = "MultivarAdjust")

  property_names <- names(property_matrices)
  n_depths <- length(depths)
  n_sims <- ncol(property_matrices[[1]])
  k <- length(property_names)

  if (n_depths < 2 || n_sims < 1) {
    log_message("WARN", "Insufficient dimensions for correlation preservation", category = "MultivarAdjust")
    return(property_matrices)
  }

  # --- Property correlation (R_prop): empirical, from the surface (first) depth's already-drawn
  # values - the same source adjust_simulation_depthwise()'s own verification step already
  # uses. Falls back to the identity (no property correlation) on a single property or an
  # estimation failure, matching this package's established graceful-degradation pattern.
  if (k < 2) {
    R_prop <- matrix(1, 1, 1, dimnames = list(property_names, property_names))
  } else {
    R_prop <- tryCatch({
      surface_data <- sapply(property_matrices, function(x) x[1, ])
      ensure_positive_definite_matrix(stats::cor(surface_data, use = "complete.obs"))
    }, error = function(e) {
      handle_workflow_error(e, "Property correlation estimation for joint copula", "warn")
      NULL
    })

    if (is.null(R_prop)) {
      R_prop <- diag(k)
      dimnames(R_prop) <- list(property_names, property_names)
    }
  }

  # --- Depth correlation (R_depth): length-scale reused from already-fitted GP models when
  # supplied (see @param gp_models above), else a full-depth-range fallback.
  length_scale <- NA_real_
  if (!is.null(gp_models) && length(gp_models) > 0) {
    scales <- vapply(property_names, function(prop) {
      if (prop %in% names(gp_models)) extract_depth_length_scale(gp_models[[prop]]) else NA_real_
    }, numeric(1))
    scales <- scales[is.finite(scales) & scales > 0]
    if (length(scales) > 0) {
      length_scale <- mean(scales)
    }
  }
  if (!is.finite(length_scale) || length_scale <= 0) {
    length_scale <- diff(range(depths))
    if (!is.finite(length_scale) || length_scale <= 0) {
      length_scale <- 1
    }
  }

  R_depth <- build_depth_correlation_kernel(
    depths, length_scale, kernel = kernel, boundary_distinctness = boundary_distinctness
  )

  # --- Joint draw + marginal mapping. property_matrices' name order drives both R_prop's
  # dimnames (built above from the same sapply()) and Z's property axis (via
  # sample_joint_depth_property_copula(R_depth, R_prop, ...)), so no reordering is needed before
  # handing Z to apply_copula_to_marginals().
  Z <- tryCatch({
    sample_joint_depth_property_copula(R_depth, R_prop, n_sims)
  }, error = function(e) {
    handle_workflow_error(e, "Joint copula sampling", "warn")
    NULL
  })

  if (is.null(Z)) {
    log_message("WARN", "Joint copula sampling failed - returning original property matrices",
                category = "MultivarAdjust")
    return(property_matrices)
  }

  adjusted_list <- tryCatch({
    apply_copula_to_marginals(Z, property_matrices, gp_predictions = gp_predictions)
  }, error = function(e) {
    handle_workflow_error(e, "Copula-to-marginal mapping", "warn")
    NULL
  })

  if (is.null(adjusted_list)) {
    log_message("WARN", "Copula-to-marginal mapping failed - returning original property matrices",
                category = "MultivarAdjust")
    return(property_matrices)
  }

  log_message("DEBUG", "Joint correlation preservation completed", category = "MultivarAdjust")
  adjusted_list
}

# ============================================================================
# 2. NRCS GP INTEGRATION FUNCTIONS
# ============================================================================

#' Apply NRCS Trend Adjustments
#'
#' Maps each requested property to its NRCS regional GP model, predicts that model's
#' depth trend for the component's depths, and applies it via apply_gp_depth_trends().
#'
#' @param cokey_data Simulation data for a single cokey
#' @param gp_models Fitted NRCS GP depth models
#' @param model_group GP model group for this cokey
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @param config Optional Monte Carlo config, passed through to `apply_gp_depth_trends()` -
#'   `config$monte_carlo$vertical_correlation_method` (`"joint_copula"` default, or
#'   `"gp_quantile_retrofit"`) reaches the NRCS/regional GP path as well as the local-GP path.
#'   `NULL` (default) resolves to `"joint_copula"`, matching `default_monte_carlo_config()`.
#' @return Adjusted simulation data
#' @family gp-modeling
#' @export
apply_nrcs_trend_adjustments <- function(cokey_data,
                                         gp_models,
                                         model_group,
                                         properties,
                                         preserve_correlations = TRUE,
                                         verbose = getOption("ssurgo.verbose", FALSE),
                                         config = NULL) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", paste("Applying NRCS GP adjustments for group:", model_group), category = "MultivarAdjust")

  # property mapping
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

  # Get NRCS GP predictions for each property using the GP depth-modeling functions
  nrcs_gp_predictions <- get_nrcs_gp_predictions(
    gp_models, available_properties, property_mapping, model_group, unique_depths, cokey_data
  )

  if (length(nrcs_gp_predictions) == 0) {
    log_message("DEBUG", "No valid NRCS GP predictions obtained", category = "MultivarAdjust")
    return(cokey_data)
  }

  # Extract the actual fitted NRCS GP model objects (not just their predictions) for the joint
  # method's depth-kernel length-scale reuse (extract_depth_length_scale()) - same
  # gp_models[[nrcs_prop]]$models[[model_group]] lookup get_nrcs_gp_predictions() already does,
  # keyed here by `prop` (the cokey_data/property_matrices name) to match
  # apply_gp_depth_trends()'s gp_models contract. Ignored entirely under the default
  # "gp_quantile_retrofit" method.
  nrcs_fitted_models <- list()
  for (prop in available_properties) {
    nrcs_prop <- property_mapping[[prop]]
    if (nrcs_prop %in% names(gp_models) &&
        isTRUE(gp_models[[nrcs_prop]]$type == "stratified_grouped") &&
        model_group %in% names(gp_models[[nrcs_prop]]$models)) {
      nrcs_fitted_models[[prop]] <- gp_models[[nrcs_prop]]$models[[model_group]]
    }
  }

  # Apply GP depth trends
  result <- apply_gp_depth_trends(
    cokey_data,
    nrcs_gp_predictions,
    available_properties,
    preserve_correlations,
    config = config,
    gp_models = nrcs_fitted_models
  )

  return(result)
}

#' Match Simulations to NRCS Models
#'
#' Looks up the NRCS GP model group for a cokey in cokey_mapping, returning
#' fallback_group when no mapping is found.
#'
#' @param cokey Target cokey
#' @param cokey_mapping Mapping from match_soils_to_gp_models()
#' @param fallback_group Default group if matching fails
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Model group name
#' @keywords internal
match_simulations_to_nrcs_models <- function(cokey, cokey_mapping, fallback_group = "general_pool",
                                             verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(cokey_mapping)) {
    log_message("DEBUG", "No cokey mapping provided - using fallback", category = "MultivarAdjust")
    return(fallback_group)
  }

  # matching with error handling
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
#' Predicts each stratified-grouped NRCS GP model's depth trend on a common depth
#' grid, returning one data frame of predicted values per property.
#'
#' @param gp_models NRCS GP models
#' @param properties Properties to extract trends for
#' @param depths Depths for trend extraction
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of depth trends by property and group
#' @keywords internal
extract_nrcs_depth_trends <- function(gp_models, properties, depths = seq(0, 200, by = 10),
                                      verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Extracting NRCS depth trends", category = "MultivarAdjust")

  # Validate inputs
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
            # call the GP depth-trend predictor
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
# 3. LOCAL GP INTEGRATION FUNCTIONS
# ============================================================================

#' Apply Local GP Adjustments
#'
#' Fits per-component Gaussian-process depth models for the requested properties and
#' applies their predicted trends to the component's simulation data.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param properties Properties to adjust
#' @param preserve_correlations Whether to preserve correlations
#' @param min_depths Minimum depths required for GP fitting
#' @param config Optional configuration list; defaults to the package validation configuration
#' @param gp_control Passed through to `fit_local_gp_models()`/`fit_local_gp_model_single()`'s
#'   `gp_control` - see `fit_local_gp_model_single()`'s docs for why the default is much smaller
#'   than `GPfit::GP_fit()`'s own default.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Adjusted simulation data
#' @family gp-modeling
#' @export
apply_local_gp_adjustments <- function(cokey_data,
                                       properties,
                                       preserve_correlations = TRUE,
                                       min_depths = 3,
                                       config = NULL,
                                       gp_control = c(20, 10, 2),
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- default_config("validation")
  }

  # depth validation
  unique_depths <- sort(unique(cokey_data$hzdept_r[!is.na(cokey_data$hzdept_r)]))

  if (length(unique_depths) < min_depths) {
    log_message("DEBUG", paste("Insufficient depths for local GP fitting:", length(unique_depths), "<", min_depths), category = "MultivarAdjust")
    return(cokey_data)
  }

  # Fit local GP models
  local_gp_models <- fit_local_gp_models(cokey_data, properties, config, gp_control = gp_control)

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

  # Apply local depth trends - config/local_gp_models threaded through so
  # config$monte_carlo$vertical_correlation_method = "joint_copula" is reachable from this function's own already-fitted local_gp_models.
  result <- apply_local_depth_trends(
    cokey_data,
    local_gp_predictions,
    unique_depths,
    preserve_correlations,
    config = config,
    gp_models = local_gp_models
  )

  return(result)
}

#' Fit Local GP Models
#'
#' Fits one Gaussian-process depth model per property from a single component's
#' depth-aggregated values.
#'
#' @param cokey_data Simulation data for a single cokey
#' @param properties Properties to model
#' @param config Configuration settings
#' @param gp_control Passed through to `fit_local_gp_model_single()`'s `gp_control` - see its
#'   docs for why the default is much smaller than `GPfit::GP_fit()`'s own default.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of fitted local GP models
#' @keywords internal
fit_local_gp_models <- function(cokey_data, properties, config = NULL, gp_control = c(20, 10, 2),
                                verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- default_config("validation")
  }

  log_message("DEBUG", "Fitting local GP models", category = "MultivarAdjust")

  local_models <- list()

  for (prop in properties) {
    if (!prop %in% names(cokey_data)) {
      log_message("DEBUG", paste("Property", prop, "not found in data"), category = "MultivarAdjust")
      next
    }

    model_result <- tryCatch({
      # aggregation
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

      # Fit GP model 
      fit_local_gp_model_single(agg_data, prop, gp_control = gp_control)

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
#' Applies locally fitted GP depth-trend predictions to a component's simulation
#' data via apply_gp_depth_trends().
#'
#' @param cokey_data Simulation data
#' @param local_predictions Local GP predictions
#' @param unique_depths Depth vector
#' @param preserve_correlations Whether to preserve correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @param config Optional config, passed straight through to `apply_gp_depth_trends()` - lets
#'   `config$monte_carlo$vertical_correlation_method` (`"joint_copula"` default, or
#'   `"gp_quantile_retrofit"`) reach this call site. `NULL` (default) resolves to `"joint_copula"`,
#'   matching `default_monte_carlo_config()`.
#' @param gp_models Optional named list of fitted local GP models (as `apply_local_gp_adjustments()`
#'   already has in scope via `fit_local_gp_models()`), passed straight through to
#'   `apply_gp_depth_trends()` so the joint-copula depth kernel can reuse their fitted
#'   length-scales. Ignored under `"gp_quantile_retrofit"`.
#' @return Adjusted simulation data
#' @keywords internal
apply_local_depth_trends <- function(cokey_data,
                                     local_predictions,
                                     unique_depths,
                                     preserve_correlations = TRUE,
                                     verbose = getOption("ssurgo.verbose", FALSE),
                                     config = NULL,
                                     gp_models = NULL) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Applying local depth trends", category = "MultivarAdjust")

  # Apply GP depth trends using the same logic as NRCS
  result <- tryCatch({
    apply_gp_depth_trends(
      cokey_data,
      local_predictions,
      names(local_predictions),
      preserve_correlations,
      config = config,
      gp_models = gp_models
    )
  }, error = function(e) {
    handle_workflow_error(e, "Local depth trend application", "warn")
    return(cokey_data)
  })

  return(result)
}

# ============================================================================
# 4. MULTIVARIATE PROCESSING FUNCTIONS
# ============================================================================

#' Convert to Property Matrices
#'
#' Reshapes long-format simulation data into a named list of depth-by-simulation
#' matrices, one per property.
#'
#' @param simulation_data Simulation data in long format
#' @param properties Properties to convert
#' @param unique_depths Depth vector
#' @param sim_numbers Simulation numbers
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Named list of property matrices
#' @family gp-modeling
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
      # Vectorized lookup: match() against a combined depth/simulation_number key takes the
      # first matching row for each (depth, sim_num) cell; unmatched or NA values collapse to
      # NA. `target_key` order matches matrix()'s column-major fill (depth varies fastest,
      # matching the row index; sim_number slowest, matching the column index).
      row_key <- paste(simulation_data$hzdept_r, simulation_data$simulation_number, sep = "\r")
      target_key <- paste(
        rep(unique_depths, times = length(sim_numbers)),
        rep(sim_numbers, each = length(unique_depths)),
        sep = "\r"
      )
      match_idx <- match(target_key, row_key)
      prop_matrix[] <- simulation_data[[prop]][match_idx]

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
#' Reshapes a named list of adjusted depth-by-simulation property matrices back into
#' long format, re-attaching component metadata.
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
#' @keywords internal
convert_to_long_format <- function(adjusted_matrices,
                                   unique_depths,
                                   sim_numbers,
                                   original_data,
                                   properties,
                                   verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Converting adjusted matrices to long format", category = "MultivarAdjust")

  metadata_cols <- intersect(
    c("cokey", "compname", "mukey", "hzdept_r", "hzdepb_r", "simulation_number", "unique_id"),
    names(original_data)
  )

  # Vectorized "first matching row" lookup, same approach as convert_to_property_matrices().
  # Grid order is depth outer/slower, sim_num inner/faster.
  result_df <- tryCatch({
    grid_depth <- rep(unique_depths, each = length(sim_numbers))
    grid_sim <- rep(sim_numbers, times = length(unique_depths))

    row_key <- paste(original_data$hzdept_r, original_data$simulation_number, sep = "\r")
    target_key <- paste(grid_depth, grid_sim, sep = "\r")
    match_idx <- match(target_key, row_key)
    keep <- !is.na(match_idx)

    if (!any(keep)) {
      data.frame()
    } else {
      out <- original_data[match_idx[keep], metadata_cols, drop = FALSE]

      # Matrix row/col indices (into unique_depths/sim_numbers) for each kept grid cell, to
      # pull the corresponding adjusted value out of each property's matrix.
      depth_idx <- rep(seq_along(unique_depths), each = length(sim_numbers))[keep]
      sim_idx <- rep(seq_along(sim_numbers), times = length(unique_depths))[keep]

      for (prop in properties) {
        if (prop %in% names(adjusted_matrices)) {
          adj_matrix <- adjusted_matrices[[prop]]
          valid <- depth_idx <= nrow(adj_matrix) & sim_idx <= ncol(adj_matrix)
          vals <- rep(NA_real_, length(depth_idx))
          vals[valid] <- adj_matrix[cbind(depth_idx[valid], sim_idx[valid])]
          out[[prop]] <- vals
        }
      }
      out
    }
  }, error = function(e) {
    handle_workflow_error(e, "Long format data binding", "warn")
    return(data.frame())
  })

  log_message("DEBUG", paste("Converted to long format with", nrow(result_df), "rows"), category = "MultivarAdjust")
  return(result_df)
}

# ============================================================================
# 5. VALIDATION AND QUALITY CONTROL
# ============================================================================

#' Validate Integration Results
#'
#' Checks row-count preservation, within-depth correlation preservation, and
#' depth-trend realism of an integration result, and returns an overall score.
#'
#' @param original_data Original simulation data
#' @param integrated_data Integrated simulation data
#' @param properties Properties that were processed
#' @param preserve_correlations Whether correlations should be preserved
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results list
#' @keywords internal
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

  # Data integrity checks
  validation_results$data_integrity <- validate_data_integrity(original_data, integrated_data, properties)

  # Correlation preservation checks
  if (preserve_correlations && length(properties) >= 2) {
    validation_results$correlation_preservation <- validate_correlation_preservation_integration(
      original_data, integrated_data, properties
    )
  }

  # Trend realism checks
  validation_results$trend_realism <- validate_depth_trends(integrated_data, properties)

  # Overall assessment
  validation_results$overall_assessment <- calculate_overall_validation_score_integration(validation_results)

  # Log validation summary
  log_validation_summary(validation_results, preserve_correlations)

  return(validation_results)
}

#' Correct Distribution Shapes
#'
#' Clamps each adjusted property to its plausible range and, where requested, remaps
#' it back onto the shape of its pre-adjustment distribution.
#'
#' @param adjusted_data Adjusted simulation data
#' @param original_data Original simulation data
#' @param properties Properties to validate
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Corrected simulation data
#' @keywords internal
correct_distribution_shapes <- function(adjusted_data, original_data, properties, config = NULL,
                                        verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(config)) {
    config <- default_config("validation")
  }

  log_message("DEBUG", "Correcting distribution shapes", category = "MultivarAdjust")

  corrected_data <- adjusted_data

  for (prop in properties) {
    if (!prop %in% names(adjusted_data)) {
      next
    }

    # Get property-specific constraints
    constraints <- get_property_constraints(prop)

    # Apply range constraints
    corrected_data[[prop]] <- apply_range_constraints(corrected_data[[prop]], constraints)

    # Apply distribution correction if needed
    if (constraints$preserve_distribution) {
      corrected_data[[prop]] <- correct_property_distribution(
        corrected_data[[prop]], original_data[[prop]], constraints
      )
    }
  }

  # Apply cross-property constraints
  corrected_data <- apply_cross_property_constraints(corrected_data, properties)

  log_message("DEBUG", "Distribution shape correction completed", category = "MultivarAdjust")
  return(corrected_data)
}

# ============================================================================
# 6. HELPER FUNCTIONS
# ============================================================================

# Property detection
detect_simulation_properties <- function(simulation_data) {

  # Use property validation
  all_properties <- predefined_properties("laboratory")

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

# parallel processing
process_cokeys_parallel <- function(cokey_groups, unique_cokeys, properties,
                                    gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                    preserve_correlations, n_cores, config) {

  run_parallel_lapply(
    unique_cokeys,
    function(cokey) {
      process_single_cokey(cokey_groups[[as.character(cokey)]], cokey, properties, gp_models,
                            cokey_mapping, use_nrcs_gp, use_local_gp,
                            preserve_correlations, config)
    },
    n_cores = n_cores,
    # process_single_cokey() can call apply_local_gp_adjustments()/apply_nrcs_trend_adjustments(),
    # which fit GPfit models - GPfit::GP_fit()'s internal hyperparameter search is a genetic
    # algorithm that genuinely uses R's RNG (confirmed empirically via future_lapply's
    # "UNRELIABLE VALUE" warning when this was FALSE), despite the fitted result itself being
    # highly stable across runs (see fit_local_gp_model_single()'s gp_control docs).
    future_seed = TRUE,
    op_name = "cokey integration",
    sequential_fallback = function() {
      process_cokeys_sequential(cokey_groups, unique_cokeys, properties,
                                 gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                 preserve_correlations, config)
    }
  )
}

# sequential processing
process_cokeys_sequential <- function(cokey_groups, unique_cokeys, properties,
                                      gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                      preserve_correlations, config) {

  log_message("INFO", "Running integration sequentially", category = "MultivarAdjust")

  results <- list()

  for (i in seq_along(unique_cokeys)) {
    cokey <- unique_cokeys[i]

    # Progress tracking
    track_progress(i, length(unique_cokeys), "Processing cokeys", update_frequency = 10)

    result <- process_single_cokey(cokey_groups[[as.character(cokey)]], cokey, properties, gp_models,
                                            cokey_mapping, use_nrcs_gp, use_local_gp,
                                            preserve_correlations, config)

    results[[i]] <- result
  }

  return(results)
}

# single cokey processing
process_single_cokey <- function(cokey_data, cokey, properties, gp_models,
                                          cokey_mapping, use_nrcs_gp, use_local_gp,
                                          preserve_correlations, config) {

  tryCatch({
    if (nrow(cokey_data) < 2) {
      log_message("DEBUG", paste("Insufficient data for cokey", cokey), category = "MultivarAdjust")
      return(NULL)
    }

    result_data <- cokey_data

    # Apply NRCS GP adjustments if requested
    if (use_nrcs_gp && !is.null(gp_models) && !is.null(cokey_mapping)) {
      model_group <- match_simulations_to_nrcs_models(cokey, cokey_mapping)

      # config is threaded through so "joint_copula" is reachable via the NRCS/regional GP path,

      # not only the local-GP branch below.
      result_data <- apply_nrcs_trend_adjustments(
        result_data, gp_models, model_group, properties, preserve_correlations, config = config
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

# result combination
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

# metadata creation
create_integration_results <- function(final_data, simulation_data, original_metadata,
                                       integration_method, properties, preserve_correlations,
                                       use_nrcs_gp, use_local_gp, n_cokeys_total,
                                       n_cokeys_successful, processing_time, parallel,
                                       validation_results) {

  # Create comprehensive metadata
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

  # Prepare final output 
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

# Additional helper functions
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

# GP ratio calculation
calculate_safe_gp_ratio <- function(gp_means, i) {
  if (is.na(gp_means[i-1]) || is.na(gp_means[i]) || gp_means[i-1] == 0) {
    return(1)  # No adjustment if invalid ratio
  } else {
    ratio <- gp_means[i] / gp_means[i - 1]
    # Clamp extreme ratios
    return(pmax(0.1, pmin(10, ratio)))
  }
}

# quantile adjustment
apply_quantile_adjustment <- function(reference_quantiles, curr_values, prev_values, gp_ratio, n_sims) {
  # Vectorized: quantile() accepts a vector of probs and computes every requested quantile
  # from a single sort of curr_values. quantile()'s failure modes (e.g. curr_values all-NA)
  # depend on curr_values as a whole, not on the individual prob requested, so one tryCatch
  # around the vectorized call suffices; on failure the row's curr_values are returned
  # unadjusted.
  tryCatch({
    quantile_values <- stats::quantile(curr_values, probs = reference_quantiles, na.rm = TRUE, names = FALSE)
    quantile_values + (prev_values * gp_ratio - quantile_values)
  }, error = function(e) {
    curr_values
  })
}

# distribution shape correction
correct_distribution_shape <- function(curr_values, adjusted_curr) {
  tryCatch({
    # A property column that's entirely (or almost entirely) NA for this group is a real,
    # reachable case (e.g. a texture column simulate_cokey_generalized() left NA for some
    # rows rather than crashing the whole cokey - see R/core-simulation.R). Without
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

#' Fit a Single Local GP Model
#'
#' @param agg_data One row per depth (`hzdept_r`, `mean_val`), as returned by
#'   `aggregate_property_by_depth()` - always a handful of points (typically <10, one per
#'   unique depth in a single cokey), never one row per Monte Carlo replicate.
#' @param prop Property name, stored on the returned model for reference.
#' @param gp_control `GPfit::GP_fit()`'s `control` argument (population size / iteration counts
#'   for its internal hyperparameter search). Defaults to `c(20, 10, 2)`, far below `GP_fit()`'s
#'   own default `c(200*d, 80*d, 2*d)` (`d` = input dimensionality, always 1 here - depth is the
#'   only predictor). `GP_fit()`'s default search effort scales with `d`, not with the number of
#'   training points, and a 1-D, single-hyperparameter (`beta`) fit on this few points has a
#'   simple enough likelihood surface that the smaller search converges to the identical
#'   optimum every time - verified empirically across several representative depth/value series
#'   (identical fitted `beta` and identical predictions vs. the default, ~5x faster per call).
#'   Pass `c(200, 80, 2)` (or larger) to restore `GP_fit()`'s own default search effort if a
#'   future property/dataset needs a more thorough search.
#' @return A list with the fitted GP model, depth scaling info, training data, and `prop`.
#' @family gp-modeling
fit_local_gp_model_single <- function(agg_data, prop, gp_control = c(20, 10, 2)) {
  depths <- agg_data$hzdept_r
  values <- agg_data$mean_val

  # Scale depths to [0,1]
  depth_min <- min(depths)
  depth_max <- max(depths)
  depth_range <- depth_max - depth_min

  if (depth_range > 0) {
    scaled_depths <- (depths - depth_min) / depth_range

    # Fit GP model
    gp_model <- GPfit::GP_fit(X = as.matrix(scaled_depths), Y = values, control = gp_control)

    # Store with scaling information
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

#' @section Behavior:
#' A row of `adjusted_data` updates `result_data` only when its (`hzdept_r`,
#' `simulation_number`) key matches **exactly one** `result_data` row (zero or multiple
#' matches are skipped, not an error) and its own value for that property is non-`NA`. When
#' multiple `adjusted_data` rows share a key, the assignment applies them in order and the last
#' one wins. This is done with a single vectorized key match rather than a per-row table scan.
merge_adjusted_data <- function(cokey_data, adjusted_data, available_properties) {
  result_data <- cokey_data

  key_result <- paste(result_data$hzdept_r, result_data$simulation_number, sep = "\a")
  key_adj <- paste(adjusted_data$hzdept_r, adjusted_data$simulation_number, sep = "\a")

  # A key must be unique within result_data (exactly one candidate row): match() alone only
  # finds the first occurrence and can't tell a duplicate from a genuine single match.
  key_counts <- table(key_result)
  unique_keys <- names(key_counts)[key_counts == 1]
  match_idx <- match(key_adj, key_result)
  match_idx[!(key_adj %in% unique_keys)] <- NA_integer_

  for (prop in available_properties) {
    if (!(prop %in% names(adjusted_data))) next
    valid <- !is.na(match_idx) & !is.na(adjusted_data[[prop]])
    if (any(valid)) {
      result_data[[prop]][match_idx[valid]] <- adjusted_data[[prop]][valid]
    }
  }

  return(result_data)
}

# validation functions
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
  # correlation validation
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
      # Check for realistic depth trends
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

# property constraints
get_property_constraints <- function(property) {
  # Use property validation to get realistic ranges
  constraints <- list(
    range = NULL,
    preserve_distribution = FALSE,
    cross_property_rules = NULL
  )

  # Property-specific constraints
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

  # quantile mapping approach
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

#' @section Behavior:
#' Rescales the texture columns so they sum to 100 where that sum is positive and non-`NA`,
#' operating on the whole texture-column matrix at once with `rowSums()`. An `NA` texture sum
#' skips rescaling for that row.
apply_cross_property_constraints <- function(data, properties) {
  # texture sum constraint
  texture_props <- intersect(c("sand_total", "sandtotal", "clay_total", "claytotal", "silt_total", "silttotal"),
                             properties)

  if (length(texture_props) >= 2) {
    tryCatch({
      # Normalize texture properties to sum to 100%
      texture_mat <- as.matrix(data[, texture_props, drop = FALSE])
      texture_sum <- rowSums(texture_mat, na.rm = TRUE)
      needs_scaling <- texture_sum > 0 & !is.na(texture_sum)
      if (any(needs_scaling)) {
        scaling_factor <- 100 / texture_sum[needs_scaling]
        data[needs_scaling, texture_props] <- texture_mat[needs_scaling, , drop = FALSE] * scaling_factor
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

    # Simple scaling approach
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
