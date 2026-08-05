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
#' @param validation_config Validation configuration from Module 0
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Processed data frame ready for GP model building
#' @export
prepare_nrcs_training_data <- function(nrcs_combined_data,
                                       grouping_strategy = "auto",
                                       min_profiles_per_group = 3,
                                       min_observations_per_group = 15,
                                       target_min_groups = 3,
                                       max_depth = 250,
                                       validation_config = NULL,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== NRCS TRAINING DATA PREPARATION ===", category = "GPModeling")
  log_message("INFO", paste("Grouping strategy:", grouping_strategy), category = "GPModeling")
  log_message("INFO", paste("Min profiles per group:", min_profiles_per_group), category = "GPModeling")
  log_message("INFO", paste("Max depth:", max_depth, "cm"), category = "GPModeling")

  # Use Module 0 for initial data validation
  if (is.null(validation_config)) {
    validation_config <- get_default_configuration("validation")
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

  # Standardize property names using Module 0
  processed_data <- standardize_property_names(
    nrcs_combined_data,
    target_standard = "ssurgo"
  )

  # Apply unsuitable horizon detection from Module 0.
  # NOTE: the original source relied on magrittr's `.` placeholder here
  # (`%>% mutate(... = is_unsuitable(., ...))`), which substitutes the
  # piped-in data frame anywhere it appears in the RHS call - not just as
  # the first argument. The base R pipe (`|>`) has no equivalent behavior
  # for `.` nested inside another call, so converting to `|>` required
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

  # Use Module 0 safe coalesce for property standardization
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
#' Enhanced version that uses Module 0 validation utilities
#'
#' @param data Input NRCS data
#' @param min_profiles Minimum profiles per group
#' @param min_obs Minimum observations per group
#' @param target_groups Minimum number of adequate groups desired
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Selected grouping strategy name
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

  # Use Module 0 column validation
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
#' Enhanced validation using Module 0 utilities
#'
#' @param processed_data Processed NRCS data
#' @param properties Properties to validate
#' @param min_profiles Minimum profiles per group
#' @param min_observations Minimum observations per group
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results list
#' @export
validate_training_groups <- function(processed_data,
                                     properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
                                     min_profiles = 3,
                                     min_observations = 15,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== TRAINING GROUP VALIDATION ===", category = "GPModeling")

  # Use Module 0 validation framework
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

  # Property coverage assessment using Module 0 safe operations
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
# 2. GP MODEL BUILDING FUNCTIONS (Enhanced with Error Handling)
# ============================================================================

#' Build Stratified GP Models
#'
#' Enhanced version with proper Module 0 error handling and validation
#'
#' @param processed_nrcs_data Processed NRCS data from prepare_nrcs_training_data()
#' @param properties Character vector of properties to model
#' @param min_profiles_per_group Minimum profiles per group (default = 3)
#' @param min_observations_per_group Minimum observations per group (default = 15)
#' @param optimize_hyperparameters Whether to optimize GP hyperparameters (default = TRUE)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List of GP models organized by property and group
#' @export
build_stratified_gp_models <- function(processed_nrcs_data,
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

  return(gp_models)
}

#' Fit Individual GP Model
#'
#' Enhanced version with better error handling and diagnostics
#'
#' @param data Input data for the group
#' @param property Property name to model
#' @param optimize_hyperparameters Whether to optimize hyperparameters
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List containing GP model and diagnostics
#' @export
fit_individual_gp_model <- function(data, property, optimize_hyperparameters = TRUE,
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

    # Check for variation using Module 0 safe operations
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
      gp_model <- optimize_gp_hyperparameters(scaled_depths, values)
    } else {
      gp_model <- GPfit::GP_fit(X = as.matrix(scaled_depths), Y = values)
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
# 3. MULTIVARIATE CORRELATION PRESERVATION (From GP-depth-adjust.R)
# ============================================================================

#' Adjust Multiple Soil Properties While Preserving Correlations
#'
#' Core function from GP-depth-adjust.R that adjusts multiple simulated soil properties
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
#' @export
adjust_multivariate_depthwise_GP <- function(simulated_list, gp_models, depths,
                                             primary_property = NULL,
                                             verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== MULTIVARIATE GP DEPTH ADJUSTMENT ===", category = "GPModeling")

  # Input validation using Module 0 utilities
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

      # Apply adjustment using REFERENCE quantiles (not property-specific)
      adjusted_curr <- numeric(n_sims)

      for (j in 1:n_sims) {
        # Use the SAME quantile from reference property
        q <- reference_quantiles[j]

        # Find corresponding quantile value in current property's distribution
        quantile_value <- quantile(curr_values, probs = q, na.rm = TRUE)

        # Apply GP-based adjustment
        adjusted_curr[j] <- quantile_value + (prev_values[j] * gp_ratio - quantile_value)
      }

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
#' Helper function from GP-depth-adjust.R to validate correlation preservation
#'
#' @param original_list Original simulated property list
#' @param adjusted_list Adjusted property list
#' @param depths Depth vector
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#'
#' @export
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

# ============================================================================
# 4. ENHANCED SOIL PROPERTY SIMULATION (From GP-depth-adjust.R)
# ============================================================================

#' Select a Correlation Matrix and Property Set for Single-Cokey Simulation
#'
#' Chooses a single flat correlation matrix (and its associated property set)
#' to hand to [generate_monte_carlo_realizations()] for one cokey's data,
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
#'   (a matrix or `NULL`, letting [generate_monte_carlo_realizations()]
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
#' [generate_monte_carlo_realizations()] (as `result$simulation_data`) into a
#' long-format data frame with one row per horizon-realization combination,
#' matching the shape [apply_local_gp_adjustments()], [apply_nrcs_trend_adjustments()],
#' and [convert_to_property_matrices()] expect (a `simulation_number` column
#' plus one column per property, alongside the original horizon metadata).
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

#' Enhanced Soil Property Simulation with NRCS GP Models + Cholesky Correlations
#'
#' Complete implementation from GP-depth-adjust.R that combines triangular distribution
#' + Cholesky approach with NRCS GP models for realistic depth trends
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
#' @return Enhanced simulation results with realistic depth trends and preserved correlations
#' @export
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

  log_message("INFO", "=== ENHANCED SOIL PROPERTY SIMULATION ===", category = "GPModeling")
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
  # Choose a single correlation matrix to hand to generate_monte_carlo_realizations():
  # correlation_matrices/txt_correlation_matrices are documented as
  # "existing correlation matrices by genetic horizon" (genhz-keyed lists),
  # but generate_monte_carlo_realizations() takes one flat matrix per call -
  # use the matrix for cokey_data's most common genhz as a reasonable
  # single-cokey simplification (a cokey's horizons are usually dominated by
  # one or two genhz groups). Falls back to auto-detected properties and no
  # explicit correlation matrix (empirical estimation) if no usable
  # genhz-keyed or flat matrix is supplied.
  initial_simulations <- tryCatch({
    chosen <- select_simulation_correlation_matrix(cokey_data, correlation_matrices, txt_correlation_matrices)

    mc_result <- generate_monte_carlo_realizations(
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

  log_message("INFO", paste("Enhanced simulations completed with", nrow(enhanced_simulations), "rows"), category = "GPModeling")

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

  log_message("INFO", "=== ENHANCED SIMULATION COMPLETE ===", category = "GPModeling")

  return(enhanced_simulations)
}

# ============================================================================
# 5. SOIL MATCHING AND MODEL APPLICATION FUNCTIONS
# ============================================================================

#' Match Soils to GP Models
#'
#' Enhanced version with better error handling and fallback strategies
#'
#' @param simulated_cokeys Vector of cokeys from simulation data
#' @param nrcs_combined_data Original NRCS data used for GP training
#' @param gp_models Fitted GP models from build_stratified_gp_models()
#' @param property Property name for matching (default = "clay_pct")
#' @param matching_strategy Matching approach (default = "exact_cokey")
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data frame mapping cokeys to GP model groups
#' @export
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
#' Enhanced version with better error handling
#'
#' @param gp_model_info GP model information from fit_individual_gp_model()
#' @param new_depths Vector of depths for prediction
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Vector of predictions
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
# 6. MODEL VALIDATION AND DIAGNOSTICS (Enhanced)
# ============================================================================

#' Validate GP Models
#'
#' Enhanced validation with comprehensive diagnostics
#'
#' @param gp_models List of GP models from build_stratified_gp_models()
#' @param nrcs_data Original NRCS training data
#' @param validation_depths Depths for prediction testing
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results
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
# 7. SUPPORTING HELPER FUNCTIONS (Leveraging Module 0)
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
#' `perform_gp_cross_validation()` (validation-diagnostics.R) so the
#' fold-splitting logic is not duplicated across files.
#'
#' @param X Predictor matrix (one row per observation).
#' @param Y Response vector.
#' @param n_folds Number of cross-validation folds.
#' @param corr_candidates Named list of `corr` specs to pass to
#'   [GPfit::GP_fit()] (e.g. `list(type = "exponential", power = 1.95)`).
#'
#' @return List with `mean_rmse_by_candidate`, a named list of mean held-out
#'   RMSE per candidate (`NA` for a candidate that failed on every fold).
k_fold_gp_cv <- function(X, Y, n_folds, corr_candidates) {
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
        model <- GPfit::GP_fit(X = X[train_idx, , drop = FALSE], Y = Y[train_idx], corr = corr_spec)
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
#'
#' @return A `"GP"`-classed [GPfit::GP_fit()] model, with a `cv_results`
#'   attribute describing the cross-validation that selected it.
optimize_gp_hyperparameters <- function(X, Y, n_folds = 5) {
  tryCatch({
    X <- as.matrix(X)
    n <- length(Y)

    corr_candidates <- list(
      exponential_power1.95 = list(type = "exponential", power = 1.95),
      matern_nu1.5 = list(type = "matern", nu = 3 / 2),
      matern_nu2.5 = list(type = "matern", nu = 5 / 2)
    )

    effective_folds <- max(2, min(n_folds, n))

    # Cross-validation needs at least 2 training points per fold to be
    # meaningful; fall back to a single baseline model otherwise (matches
    # the previous stub's behavior for very small inputs).
    cv <- if (n >= 2 * effective_folds) {
      k_fold_gp_cv(X, Y, effective_folds, corr_candidates)
    } else {
      NULL
    }

    rmses <- if (!is.null(cv)) unlist(cv$mean_rmse_by_candidate) else numeric(0)

    if (length(rmses) > 0 && any(is.finite(rmses))) {
      best_name <- names(which.min(rmses))
      best_corr <- corr_candidates[[best_name]]

      final_model <- tryCatch({
        GPfit::GP_fit(X = X, Y = Y, corr = best_corr)
      }, error = function(e) {
        GPfit::GP_fit(X = X, Y = Y)
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
    baseline_model <- GPfit::GP_fit(X = X, Y = Y)
    attr(baseline_model, "cv_results") <- list(
      folds = effective_folds,
      corr_candidates = names(corr_candidates),
      best_candidate = NA_character_,
      mean_rmse_by_candidate = list()
    )
    return(baseline_model)

  }, error = function(e) {
    log_message("WARN", paste("Hyperparameter optimization failed, using default:", e$message), category = "GPModeling")
    return(GPfit::GP_fit(X = as.matrix(X), Y = Y))
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

# Trend assessment functions (using Module 0 safe operations)
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
#' (`multivariate-adjustment.R`), which already implements NRCS GP-model
#' depth-trend adjustment with correlation preservation. This wrapper's only
#' job is to infer the `properties` argument that function requires (as the
#' numeric, non-structural columns of `simulation_data`), since this
#' function's own signature does not accept one explicitly.
#'
#' @param simulation_data Simulated soil property data for one cokey.
#' @param nrcs_gp_models Fitted NRCS GP models (see [build_stratified_gp_models()]).
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
#' Delegates to the real [apply_local_gp_adjustments()] (`multivariate-adjustment.R`,
#' migrated from mod07), inferring the `properties` argument that function
#' requires (as the numeric, non-structural columns of `simulation_data`),
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
