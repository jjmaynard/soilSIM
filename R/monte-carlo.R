#' @title Monte Carlo Simulation Engine
#' @description Advanced Monte Carlo simulation for correlated soil property generation,
#'   including flexible per-property distribution families (via `distributions.R`'s
#'   percentile-triplet fitting) and compositional (ILR-based) texture handling.
#' @name monte_carlo_simulation
NULL

# ============================================================================
# 1. MAIN SIMULATION FUNCTIONS
# ============================================================================

#' Generate Monte Carlo Realizations of Soil Properties
#'
#' Master function for generating Monte Carlo realizations
#' for validation, error handling, logging, and data quality assessment.
#'
#' @param soil_data Data frame with soil property data (must have _l, _r, _h columns)
#' @param properties Character vector of properties to simulate
#' @param correlation_matrix Optional correlation matrix for properties
#' @param n_realizations Integer, number of Monte Carlo realizations (default = 1000)
#' @param simulation_config List of simulation configuration parameters
#' @param validate_inputs Logical, whether to validate inputs (default = TRUE)
#' @param parallel Logical, whether to use parallel processing (default = FALSE)
#' @param seed Integer, random seed for reproducibility
#' @param observed_data Optional named list keyed by property name (plus, for the
#'   texture group, `claytotal`/`sandtotal`/`silttotal` supplied together)
#'   providing a Bayesian-updating likelihood to fuse against each property's
#'   SSURGO-derived prior before simulating - see `fuse_observed_data_into_priors()`
#'   for the accepted shapes and behavior. `NULL` (default) skips fusion
#'   entirely, preserving today's prior-only behavior exactly.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List containing comprehensive simulation results
#'
#' @examples
#' \dontrun{
#' # Basic simulation
#' mc_results <- generate_monte_carlo_realizations(
#'   soil_data = soil_data,
#'   properties = c("sandtotal", "claytotal", "silttotal"),
#'   n_realizations = 1000
#' )
#'
#' # Advanced simulation with configuration
#' mc_results <- generate_monte_carlo_realizations(
#'   soil_data = soil_data,
#'   properties = c("sandtotal", "claytotal", "dbovendry"),
#'   correlation_matrix = cor_matrix,
#'   n_realizations = 5000,
#'   simulation_config = list(
#'     max_depth = 200,
#'     constraint_rules = "strict",
#'     distribution_type = "triangular",
#'     quality_threshold = 0.8
#'   ),
#'   parallel = TRUE,
#'   seed = 12345
#' )
#' }
#'
#' @export
generate_monte_carlo_realizations <- function(soil_data,
                                              properties,
                                              correlation_matrix = NULL,
                                              n_realizations = 1000,
                                              simulation_config = list(),
                                              validate_inputs = TRUE,
                                              parallel = FALSE,
                                              seed = NULL,
                                              observed_data = NULL,
                                              verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  start_time <- Sys.time()

  # Initialize logging if not already done
  if (is.null(getOption("soil_workflow_log_config"))) {
    setup_logging(log_level = "INFO")
  }

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    log_message("INFO", paste("Random seed set to:", seed), category = "MonteCarlo")
  }

  log_message("INFO", "=== MONTE CARLO SIMULATION STARTED ===", category = "MonteCarlo")

  # Merge with default configuration. normalize_monte_carlo_config() accepts
  # either a flat simulation_config (e.g. list(distribution_type="normal")) -
  # this function's own @examples have always shown that shape - or one
  # properly nested under monte_carlo; see that function's docs for why the
  # flat shape previously had NO effect at all.
  default_config <- get_monte_carlo_defaults()
  simulation_config <- normalize_monte_carlo_config(simulation_config, default_config)
  config <- merge_configurations(default_config, simulation_config)

  # Validate configuration
  config_validation <- validate_monte_carlo_config(config, n_realizations)
  if (!config_validation$valid) {
    handle_workflow_error(
      list(message = paste("Configuration validation failed:", paste(config_validation$errors, collapse = ", "))),
      context = "Configuration Validation",
      recovery_action = "stop"
    )
  }

  log_message("INFO", paste("Properties:", paste(properties, collapse = ", ")), category = "MonteCarlo")
  log_message("INFO", paste("Realizations:", n_realizations), category = "MonteCarlo")
  log_message("INFO", paste("Max depth:", config$monte_carlo$max_depth, "cm"), category = "MonteCarlo")

  # Step 1: Comprehensive input validation
  if (validate_inputs) {
    log_message("INFO", "Step 1: Validating input data", category = "MonteCarlo")

    input_validation <- validate_monte_carlo_inputs(
      soil_data = soil_data,
      properties = properties,
      correlation_matrix = correlation_matrix,
      n_realizations = n_realizations,
      config = config
    )

    if (!input_validation$valid) {
      handle_workflow_error(
        list(message = paste("Input validation failed:", paste(input_validation$errors, collapse = ", "))),
        context = "Input Validation",
        recovery_action = config$monte_carlo$error_recovery_action
      )
    }

    # Log validation warnings
    if (length(input_validation$warnings) > 0) {
      for (warning in input_validation$warnings) {
        log_message("WARN", warning, category = "MonteCarlo")
      }
    }
  }

  # Step 2: Property validation
  log_message("INFO", "Step 2: Validating soil properties", category = "MonteCarlo")

  property_validation <- validate_properties_with_synonyms(
    properties,
    property_lookup = "ssurgo",
    strict_mode = config$monte_carlo$strict_property_validation
  )

  if (!property_validation$valid) {
    handle_workflow_error(
      list(message = paste("Property validation failed:", paste(property_validation$errors, collapse = ", "))),
      context = "Property Validation",
      recovery_action = config$monte_carlo$error_recovery_action
    )
  }

  # Step 2.5: Resolve composition groups (e.g. texture: sand/clay/silt -> ilr1/ilr2)
  # `sim_properties` (possibly ILR-substituted) drives every step through the
  # simulation core; `properties` (the caller-facing, original vector) is used
  # everywhere else (validation, constraints, diagnostics, metadata).
  composition_plan <- resolve_composition_groups(properties, config)
  sim_properties <- composition_plan$sim_properties

  # Step 3: Data preparation
  log_message("INFO", "Step 3: Preparing simulation data", category = "MonteCarlo")

  simulation_data <- prepare_simulation_data(
    soil_data = soil_data,
    properties = sim_properties,
    config = config,
    composition_plan = composition_plan
  )

  if (nrow(simulation_data) == 0) {
    handle_workflow_error(
      list(message = "No suitable horizons found for simulation after filtering"),
      context = "Data Preparation",
      recovery_action = "stop"
    )
  }

  log_message("INFO", paste("Suitable horizons for simulation:", nrow(simulation_data)), category = "MonteCarlo")

  # Step 4: Parameter extraction and validation
  log_message("INFO", "Step 4: Extracting simulation parameters", category = "MonteCarlo")

  simulation_params <- prepare_simulation_parameters(
    simulation_data = simulation_data,
    properties = sim_properties,
    config = config,
    composition_plan = composition_plan
  )

  # Step 5: Correlation structure configuration
  log_message("INFO", "Step 5: Configuring correlation structure", category = "MonteCarlo")

  # An explicit user-supplied correlation_matrix is necessarily built against
  # the ORIGINAL property names (e.g. sandtotal/claytotal/silttotal) - there's
  # no principled way to remap a user-supplied 3-variable sub-block onto a
  # 2-variable ILR pair, so fall back to independent per-property simulation
  # for the texture group specifically in that combination (auto-estimation
  # has no such problem - it computes the matrix FROM sim_properties directly).
  active_group_with_user_matrix <- !is.null(correlation_matrix) &&
    any(vapply(composition_plan$groups, function(g) isTRUE(g$active), logical(1)))
  if (active_group_with_user_matrix) {
    log_message(
      "WARN",
      paste("A composition group is active and an explicit correlation_matrix was supplied;",
            "no principled remap exists from raw property names to ILR pseudo-properties -",
            "falling back to independent per-property simulation for the affected group(s)."),
      category = "MonteCarlo"
    )
    composition_plan <- list(sim_properties = properties, groups = lapply(composition_plan$groups, function(g) {
      g$active <- FALSE
      g
    }))
    sim_properties <- composition_plan$sim_properties
    simulation_params <- prepare_simulation_parameters(simulation_data, sim_properties, config, composition_plan)
  }

  # Step 4.5: Fuse observed field/lab data into SSURGO-derived priors (optional).
  # Deliberately placed AFTER the active_group_with_user_matrix recompute above
  # (not immediately after Step 4) - that recompute replaces simulation_params
  # wholesale, so fusing before it would have the fused posteriors silently
  # discarded and overwritten with fresh unfused priors.
  if (!is.null(observed_data)) {
    log_message("INFO", "Step 4.5: Fusing observed data into priors", category = "MonteCarlo")
    simulation_params <- fuse_observed_data_into_priors(
      simulation_params = simulation_params,
      simulation_data = simulation_data,
      sim_properties = sim_properties,
      properties = properties,
      observed_data = observed_data,
      composition_plan = composition_plan,
      config = config
    )
  }

  correlation_config <- configure_correlation_structure(
    simulation_params = simulation_params,
    properties = sim_properties,
    correlation_matrix = correlation_matrix,
    config = config,
    simulation_data = simulation_data
  )

  # Step 6: Distribution setup
  log_message("INFO", "Step 6: Setting up distributions", category = "MonteCarlo")

  distribution_setup <- setup_distributions(
    simulation_params = simulation_params,
    properties = sim_properties,
    config = config
  )

  # Step 7: Run Monte Carlo simulation
  log_message("INFO", "Step 7: Running Monte Carlo simulation", category = "MonteCarlo")

  track_progress(0, n_realizations, "Monte Carlo simulation")

  simulation_results <- tryCatch({
    if (parallel && n_realizations > config$monte_carlo$parallel_threshold) {
      run_parallel_simulation(
        simulation_params = simulation_params,
        distribution_setup = distribution_setup,
        correlation_config = correlation_config,
        n_realizations = n_realizations,
        properties = sim_properties,
        config = config
      )
    } else {
      run_sequential_simulation(
        simulation_params = simulation_params,
        distribution_setup = distribution_setup,
        correlation_config = correlation_config,
        n_realizations = n_realizations,
        properties = sim_properties,
        config = config
      )
    }
  }, error = function(e) {
    handle_workflow_error(e, "Monte Carlo Simulation", config$monte_carlo$error_recovery_action)
  })

  # Step 7.5: Collapse any simulated ILR pseudo-properties back to their raw
  # composition members. Runs ONCE on the fully-combined array (after
  # run_parallel_simulation()'s chunk-concatenation, not per-worker) since
  # ilr_inverse() is per-realization independent - correct and cheaper done
  # once here than n_cores times inside the parallel path.
  simulation_results <- restore_composition_properties(
    simulation_results, sim_properties, properties, composition_plan$groups
  )

  # Step 8: Apply constraints using enhanced constraint system
  log_message("INFO", "Step 8: Applying simulation constraints", category = "MonteCarlo")

  constrained_results <- apply_simulation_constraints(
    simulation_results = simulation_results,
    properties = properties,
    config = config,
    composition_plan = composition_plan
  )

  # Step 9: Comprehensive output validation
  log_message("INFO", "Step 9: Validating simulation output", category = "MonteCarlo")

  output_validation <- validate_simulation_output(
    simulation_results = constrained_results,
    properties = properties,
    config = config
  )

  # Step 10: Generate comprehensive diagnostics
  log_message("INFO", "Step 10: Generating simulation diagnostics", category = "MonteCarlo")

  diagnostics <- generate_simulation_diagnostics(
    original_data = simulation_data,
    simulation_results = constrained_results,
    properties = properties,
    correlation_config = correlation_config,
    config = config
  )

  # Step 11: Quality assessment
  log_message("INFO", "Step 11: Assessing simulation quality", category = "MonteCarlo")

  quality_assessment <- assess_simulation_quality(
    simulation_results = constrained_results,
    diagnostics = diagnostics,
    output_validation = output_validation,
    config = config
  )

  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "secs")

  # Compile comprehensive results
  final_results <- list(
    simulation_data = constrained_results,
    original_data = simulation_data,
    parameters = simulation_params,
    correlation_structure = correlation_config,
    distribution_setup = distribution_setup,
    validation = list(
      input_validation = if(validate_inputs) input_validation else NULL,
      property_validation = property_validation,
      output_validation = output_validation
    ),
    diagnostics = diagnostics,
    quality_assessment = quality_assessment,
    metadata = list(
      properties = properties,
      n_realizations = n_realizations,
      n_horizons = nrow(simulation_data),
      max_depth = config$monte_carlo$max_depth,
      processing_time = processing_time,
      parallel = parallel,
      seed = seed,
      timestamp = end_time,
      config_used = config,
      success_rate = output_validation$success_rate
    )
  )

  log_message("INFO", "=== MONTE CARLO SIMULATION COMPLETE ===", category = "MonteCarlo")
  log_message("INFO", paste("Processing time:", round(as.numeric(processing_time), 2), "seconds"), category = "MonteCarlo")
  log_message("INFO", paste("Success rate:", round(output_validation$success_rate * 100, 1), "%"), category = "MonteCarlo")
  log_message("INFO", paste("Overall quality score:", round(quality_assessment$overall_quality_score, 3)), category = "MonteCarlo")

  return(final_results)
}

#' Simulate Correlated Properties Using Enhanced Distribution Framework
#'
#' Core simulation function leveraging Module 0 statistical utilities
#'
#' @param simulation_params List containing simulation parameters for each horizon
#' @param correlation_matrix Correlation matrix for the properties
#' @param n_realizations Number of realizations to generate
#' @param properties Character vector of property names
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Array of simulated values, dimensioned horizons by properties by realizations
#'
#' @export
simulate_correlated_properties <- function(simulation_params,
                                                    correlation_matrix,
                                                    n_realizations,
                                                    properties,
                                                    config,
                                                    verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  n_horizons <- length(simulation_params)
  n_properties <- length(properties)

  log_message("DEBUG", paste("Simulating", n_horizons, "horizons,", n_properties, "properties,", n_realizations, "realizations"), category = "MonteCarlo")

  # Initialize result array
  simulated_values <- array(
    dim = c(n_horizons, n_properties, n_realizations),
    dimnames = list(
      horizon = 1:n_horizons,
      property = properties,
      realization = 1:n_realizations
    )
  )

  # Prepare correlation matrix
  if (is.null(correlation_matrix)) {
    correlation_matrix <- diag(n_properties)
    rownames(correlation_matrix) <- colnames(correlation_matrix) <- properties
  }

  # Use Module 0's safe correlation utilities
  correlation_matrix <- ensure_positive_definite_matrix(correlation_matrix)

  # Validate matrix
  matrix_validation <- validate_correlation_matrix(correlation_matrix, properties)
  if (!matrix_validation$valid) {
    log_message("WARN", paste("Correlation matrix issues:", matrix_validation$message), category = "MonteCarlo")
    # Fall back to identity matrix
    correlation_matrix <- diag(n_properties)
    rownames(correlation_matrix) <- colnames(correlation_matrix) <- properties
  }

  chol_matrix <- chol(correlation_matrix)

  # Generate correlated values for each horizon with progress tracking
  for (h in seq_len(n_horizons)) {
    if (h %% max(1, floor(n_horizons / 10)) == 0) {
      track_progress(h, n_horizons, "Simulating horizons")
    }

    horizon_params <- simulation_params[[h]]

    # Generate independent uniform random values
    uniform_matrix <- matrix(
      runif(n_properties * n_realizations),
      nrow = n_properties,
      ncol = n_realizations
    )

    # Transform to standard normal
    normal_matrix <- qnorm(uniform_matrix)

    # Apply correlation structure
    correlated_normal <- chol_matrix %*% normal_matrix

    # Transform back to uniform [0,1]
    correlated_uniform <- pnorm(correlated_normal)

    # Transform to appropriate distributions for each property
    for (p in seq_len(n_properties)) {
      prop_name <- properties[p]

      if (prop_name %in% names(horizon_params)) {
        params <- horizon_params[[prop_name]]

        # family travels WITH params (resolved once, per-property/per-horizon,
        # by extract_property_parameters()/prepare_simulation_parameters()) -
        # config$monte_carlo$distribution_type is only ever consulted as
        # THEIR fallback, not here.
        simulated_values[h, p, ] <- quantile_from_fit(
          correlated_uniform[p, ],
          family = params$family,
          fit = params$fit
        )
      } else {
        # Handle missing parameters
        simulated_values[h, p, ] <- NA
        log_message("DEBUG", paste("Missing parameters for", prop_name, "in horizon", h), category = "MonteCarlo")
      }
    }
  }

  # Final validation of simulated values
  validation_summary <- validate_simulated_array(simulated_values, properties, config)
  log_message("DEBUG", paste("Simulation validation - Finite values:",
                             round(validation_summary$finite_rate * 100, 1), "%"), category = "MonteCarlo")

  return(simulated_values)
}

#' Simulate Component Compositions (Enhanced)
#'
#' Enhanced component composition simulation
#'
#' @param component_data Data frame with component data
#' @param n_realizations Number of realizations
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return List with simulated component compositions and quality metrics
#'
#' @export
sim_component_compositions <- function(component_data,
                                                n_realizations = 1000,
                                                config = NULL,
                                                verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Use Module 0 configuration if not provided
  if (is.null(config)) {
    config <- get_monte_carlo_defaults()
  }

  log_message("INFO", paste("Simulating component compositions for", nrow(component_data), "components"), category = "MonteCarlo")

  # Validate component data
  component_validation <- validate_component_data(component_data)
  if (!component_validation$valid) {
    handle_workflow_error(
      list(message = paste("Component data validation failed:", component_validation$message)),
      context = "Component Validation",
      recovery_action = "stop"
    )
  }

  n_components <- nrow(component_data)

  # Generate distributions for each component
  component_realizations <- matrix(
    nrow = n_components,
    ncol = n_realizations
  )

  for (i in seq_len(n_components)) {
    # Handle missing values
    comp_data <- extract_component_parameters(component_data[i, ], config)

    # Generate values using configured distribution type
    component_realizations[i, ] <- generate_component_values(
      parameters = comp_data,
      n_realizations = n_realizations,
      distribution_type = config$monte_carlo$distribution_type
    )
  }

  # Normalize to ensure sum = 100% using enhanced normalization
  normalized_realizations <- normalize_component_realizations(
    component_realizations,
    target_sum = 100,
    method = config$monte_carlo$normalization_method
  )

  # Apply composition constraints if specified
  if (!is.null(config$monte_carlo$composition_constraints)) {
    normalized_realizations <- apply_composition_constraints(
      normalized_realizations,
      config$monte_carlo$composition_constraints
    )
  }

  # Quality assessment
  quality_metrics <- assess_component_quality(
    original_data = component_data,
    simulated_data = normalized_realizations,
    config = config
  )

  result <- list(
    realizations = t(normalized_realizations),
    component_data = component_data,
    n_realizations = n_realizations,
    n_components = n_components,
    quality_metrics = quality_metrics,
    validation = component_validation,
    constraints_applied = !is.null(config$monte_carlo$composition_constraints),
    metadata = list(
      distribution_type = config$monte_carlo$distribution_type,
      normalization_method = config$monte_carlo$normalization_method,
      timestamp = Sys.time()
    )
  )

  log_message("INFO", paste("Component simulation complete - Quality score:",
                            round(quality_metrics$overall_quality, 3)), category = "MonteCarlo")

  return(result)
}

# ============================================================================
# 2. ENHANCED DISTRIBUTION AND PARAMETER FUNCTIONS
# ============================================================================

#' Setup Distributions (Enhanced)
#'
#' Enhanced distribution setup
#'
#' @param simulation_params List of simulation parameters by horizon
#' @param properties Character vector of property names
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Enhanced distribution configuration
#'
#' @export
setup_distributions <- function(simulation_params, properties, config, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Setting up enhanced distribution framework", category = "MonteCarlo")

  distribution_config <- list(
    distributions = list(),
    summary_stats = list(),
    validation_results = list(),
    configuration = config$monte_carlo
  )

  n_horizons <- length(simulation_params)
  distribution_type <- config$monte_carlo$distribution_type

  # Process each horizon
  for (h in seq_len(n_horizons)) {
    horizon_params <- simulation_params[[h]]
    horizon_distributions <- list()

    for (prop in properties) {
      if (prop %in% names(horizon_params)) {
        params <- horizon_params[[prop]]

        # Enhanced parameter validation
        validation <- validate_distribution_parameters(
          parameters = params,
          distribution_type = distribution_type,
          property_name = prop
        )

        if (validation$valid) {
          horizon_distributions[[prop]] <- create_distribution_config(
            parameters = params,
            distribution_type = distribution_type,
            property_name = prop,
            validation = validation
          )
        } else {
          # Create fallback distribution
          horizon_distributions[[prop]] <- create_fallback_distribution(
            property_name = prop,
            representative_value = params$mode,
            distribution_type = distribution_type
          )
        }
      }
    }

    distribution_config$distributions[[h]] <- horizon_distributions
  }

  # Generate comprehensive summary  statistical utilities
  distribution_config$summary_stats <- summarize_distributions(
    distribution_config$distributions, properties, distribution_type
  )

  # Overall validation
  distribution_config$validation_results <- validate_distribution_setup(
    distribution_config$distributions, properties, config
  )

  log_message("DEBUG", paste("Distribution setup complete -",
                             distribution_config$validation_results$valid_distributions,
                             "valid out of",
                             distribution_config$validation_results$total_distributions),
              category = "MonteCarlo")

  return(distribution_config)
}

#' Prepare Simulation Parameters (Enhanced)
#'
#' Enhanced parameter preparation. Composition-group pseudo-properties (e.g.
#' `"ilr1"`/`"ilr2"`) are special-cased: BOTH are fit together, once per
#' horizon, via `distributions.R`'s `estimate_ilr_moments_mc()` from the
#' group's real member `_l/_r/_h` triplets - their joint covariance only
#' makes sense computed jointly, not by looping `extract_property_parameters()`
#' per pseudo name. The Monte Carlo covariance's ilr1<->ilr2 OFF-DIAGONAL is
#' intentionally discarded (only the marginal SDs are kept): the cross-horizon
#' empirical correlation-matrix step (`estimate_property_correlations()`) owns
#' ilr1<->ilr2 correlation instead, so the two sources of correlation
#' (per-observation uncertainty propagation vs. cross-sample empirical
#' correlation) are never double-counted.
#'
#' @param simulation_data Filtered simulation data
#' @param properties Properties to extract parameters for (may include
#'   composition-group pseudo-properties)
#' @param config Configuration settings
#' @param composition_plan Optional result of `resolve_composition_groups()`.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Enhanced simulation parameters
#'
#' @export
prepare_simulation_parameters <- function(simulation_data, properties, config, composition_plan = NULL,
                                           verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", paste("Extracting parameters for", length(properties), "properties from",
                             nrow(simulation_data), "horizons"), category = "MonteCarlo")

  n_horizons <- nrow(simulation_data)
  simulation_params <- vector("list", n_horizons)

  active_groups <- Filter(function(g) isTRUE(g$active), (composition_plan$groups %||% list()))
  pseudo_all <- unlist(lapply(active_groups, `[[`, "pseudo"))
  ordinary_properties <- setdiff(properties, pseudo_all)
  ilr_z <- qnorm((config$monte_carlo$lh_percentile %||% c(0.05, 0.95))[2])

  for (i in seq_len(n_horizons)) {
    horizon_data <- simulation_data[i, ]
    horizon_params <- list()

    for (prop in ordinary_properties) {
      tryCatch({
        # Extract parameters with enhanced error handling
        param_extraction <- extract_property_parameters(
          horizon_data = horizon_data,
          property_name = prop,
          config = config
        )

        if (param_extraction$valid) {
          horizon_params[[prop]] <- param_extraction$parameters
        }

      }, error = function(e) {
        log_message("DEBUG", paste("Error extracting", prop, "for horizon", i, ":", e$message),
                    category = "MonteCarlo")
      })
    }

    for (group in active_groups) {
      tryCatch({
        get_lrh <- function(m) c(
          l = horizon_data[[paste0(m, "_l")]],
          r = horizon_data[[paste0(m, "_r")]],
          h = horizon_data[[paste0(m, "_h")]]
        )
        lrh <- lapply(group$members, get_lrh)  # role order per group$members (default: sand, silt, clay)

        if (all(vapply(lrh, function(x) all(is.finite(x)), logical(1)))) {
          moments <- estimate_ilr_moments_mc(
            low_clay = lrh[[1]][["l"]], rep_clay = lrh[[1]][["r"]], high_clay = lrh[[1]][["h"]],
            low_sand = lrh[[2]][["l"]], rep_sand = lrh[[2]][["r"]], high_sand = lrh[[2]][["h"]],
            low_silt = lrh[[3]][["l"]], rep_silt = lrh[[3]][["r"]], high_silt = lrh[[3]][["h"]],
            z = ilr_z
          )
          horizon_params[[group$pseudo[1]]] <- list(
            family = "normal",
            fit = list(mean = unname(moments$mu[1]), sd = sqrt(unname(moments$Sigma[1, 1]))),
            source = "ilr_mc"
          )
          horizon_params[[group$pseudo[2]]] <- list(
            family = "normal",
            fit = list(mean = unname(moments$mu[2]), sd = sqrt(unname(moments$Sigma[2, 2]))),
            source = "ilr_mc"
          )
        }
      }, error = function(e) {
        log_message("DEBUG", paste("Error estimating ILR moments for horizon", i, ":", e$message),
                    category = "MonteCarlo")
      })
    }

    simulation_params[[i]] <- horizon_params
  }

  # Generate extraction summary
  extraction_summary <- summarize_parameter_extraction(
    simulation_params, properties, n_horizons
  )

  # Log extraction results
  for (prop in properties) {
    summary <- extraction_summary[[prop]]
    log_message("INFO", sprintf("%-15s: %3d complete, %3d estimated, %3d missing (%.1f%% coverage)",
                                prop, summary$complete, summary$estimated, summary$missing,
                                summary$coverage_rate * 100), category = "MonteCarlo")
  }

  return(simulation_params)
}

# ============================================================================
# 1.5 BAYESIAN-UPDATING WIRING (bridges bayesian-updating.R's standalone
# fusion primitives into the simulation pipeline; those primitives themselves
# are unchanged and remain independently usable/testable)
# ============================================================================

#' Fuse Observed Field/Lab Data into SSURGO-Derived Priors
#'
#' Orchestrates Bayesian updating against each property's per-horizon prior
#' fit produced by `prepare_simulation_parameters()`, using
#' `bayesian-updating.R`'s pure fusion primitives (`bayes_fuse()`,
#' `fuse_property()`, `fuse_texture_group_from_triplets()`). For every
#' property present in `observed_data`, every horizon's existing prior is
#' replaced by its fused posterior - each horizon keeps its OWN prior
#' (so posteriors still vary per horizon), but all incorporate the SAME
#' shared observed-data likelihood, since field/lab measurements aren't
#' typically tied to individual simulated horizons at that granularity.
#'
#' Two supported shapes per non-texture entry in `observed_data`:
#'   - a numeric vector of raw samples -> the fully general route
#'     (`fuse_property()`'s `bayesian_update()` grid-KDE path). The prior is
#'     first SAMPLED (via `quantile_from_fit()`, since that route needs both
#'     sides to already be raw vectors - it fits no parametric family to
#'     either side), and the resulting posterior sample is stored as an
#'     exact empirical quantile function (`family = "linear_cdf"`) rather
#'     than forced back into a possibly ill-fitting parametric family.
#'   - a family-native parameter list (`list(mean=, sd=)` for a normal or
#'     lognormal prior/likelihood; `list(shape1=, shape2=)` or
#'     `list(mean=, var=)` for a beta prior/likelihood) -> the closed-form
#'     route (`bayes_fuse()`), only supported for prior families
#'     `"normal"`/`"lognormal"`/`"beta"` (the only families with a conjugate
#'     route). Any other prior family (`"triangular"`/`"uniform"`/
#'     `"metalog"`/`"linear_cdf"`) SKIPS fusion for that (horizon, property)
#'     with a logged WARN, continuing to simulate from the unfused prior
#'     rather than erroring the whole pipeline.
#'
#' The texture composition group (if active) is handled jointly: supply all
#' three of `claytotal`/`sandtotal`/`silttotal` in `observed_data` (each
#' either a `c(low=,rep=,high=)` triplet, an unnamed length-3
#' `c(low,rep,high)` vector, or a raw sample vector reduced to one via its
#' own empirical quantiles at `config$monte_carlo$lh_percentile`) - partial
#' texture entries (1-2 of 3) skip texture fusion entirely with a WARN,
#' matching the same "some but not all composition members present"
#' fallback `resolve_composition_groups()` already uses. Fusion runs once
#' per horizon via `fuse_texture_group_from_triplets()` (that horizon's own
#' clay/sand/silt l/r/h as the prior side, the shared observed triplet as
#' the likelihood side), replacing
#' `simulation_params[[i]][["ilr1"]]`/`[["ilr2"]]`.
#'
#' Note: this only narrows each property's MARGINAL prior spread. It does
#' NOT change the cross-property correlation/dependency structure used by
#' the Cholesky copula in `simulate_correlated_properties()` -
#' `estimate_property_correlations()` (Step 5, when `auto_correlation=TRUE`)
#' estimates that from `simulation_data`'s raw representative values, which
#' fusion does not touch.
#'
#' @param simulation_params Per-horizon parameter list from
#'   `prepare_simulation_parameters()`.
#' @param simulation_data The filtered per-horizon data frame (same row
#'   order as `simulation_params`).
#' @param sim_properties Character vector `simulation_params` is keyed over
#'   (may include composition-group pseudo-properties).
#' @param properties The caller-facing, original property vector. Unused
#'   directly (fusion is keyed by `sim_properties`/`observed_data` names);
#'   kept for interface symmetry with sibling pipeline functions.
#' @param observed_data Named list keyed by property name (see shapes above).
#' @param composition_plan Result of `resolve_composition_groups()`.
#' @param config Simulation configuration.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return `simulation_params`, with fused entries replaced in place.
#' @export
fuse_observed_data_into_priors <- function(simulation_params, simulation_data, sim_properties,
                                            properties, observed_data, composition_plan, config,
                                            verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  n_horizons <- length(simulation_params)
  lh_probs <- config$monte_carlo$lh_percentile %||% c(0.05, 0.95)
  n_samples <- config$monte_carlo$bayesian_fusion_n_samples %||% 1000
  supported_closed_form <- c("normal", "lognormal", "beta")

  # --- texture (composition group) fusion, jointly via ILR -------------------
  active_groups <- Filter(function(g) isTRUE(g$active), (composition_plan$groups %||% list()))
  for (group in active_groups) {
    present <- intersect(group$members, names(observed_data))
    if (length(present) == 0) next
    if (length(present) < length(group$members)) {
      log_message(
        "WARN",
        paste0("observed_data supplies only ", length(present), "/", length(group$members),
               " texture members (", paste(present, collapse = ", "),
               ") - all three are required for joint fusion; skipping texture fusion."),
        category = "MonteCarlo"
      )
      next
    }

    lik_triplets <- setNames(
      lapply(group$members, function(m) as_lrh_triplet(observed_data[[m]], lh_probs)),
      c("clay", "sand", "silt")
    )
    z_shared <- qnorm(lh_probs[2])

    for (i in seq_len(n_horizons)) {
      horizon_data <- simulation_data[i, ]
      get_lrh <- function(m) c(
        l = horizon_data[[paste0(m, "_l")]],
        r = horizon_data[[paste0(m, "_r")]],
        h = horizon_data[[paste0(m, "_h")]]
      )
      lrh <- lapply(group$members, get_lrh)
      if (!all(vapply(lrh, function(x) all(is.finite(x)), logical(1)))) next

      prior_triplets <- setNames(
        lapply(lrh, function(x) c(x[["l"]], x[["r"]], x[["h"]])),
        c("clay", "sand", "silt")
      )

      fused <- tryCatch(
        fuse_texture_group_from_triplets(
          prior_triplets = prior_triplets, lik_triplets = lik_triplets,
          z_prior = z_shared, z_lik = z_shared, n_mc = 2000, n_samples = n_samples
        ),
        error = function(e) {
          log_message("WARN", paste("Texture fusion failed at horizon", i, ":", conditionMessage(e)),
                      category = "MonteCarlo")
          NULL
        }
      )
      if (is.null(fused)) next

      simulation_params[[i]][[group$pseudo[1]]] <- list(
        family = "normal",
        fit = list(mean = unname(fused$ilr_mu[1]), sd = sqrt(unname(fused$ilr_Sigma[1, 1]))),
        source = "ilr_mc_posterior"
      )
      simulation_params[[i]][[group$pseudo[2]]] <- list(
        family = "normal",
        fit = list(mean = unname(fused$ilr_mu[2]), sd = sqrt(unname(fused$ilr_Sigma[2, 2]))),
        source = "ilr_mc_posterior"
      )
    }
  }

  # --- ordinary (non-texture) property fusion --------------------------------
  # Only ACTIVE groups' members are excluded here - if a composition group is
  # inactive (e.g. only 2 of 3 texture members were requested in `properties`),
  # its members remain ordinary sim_properties and are still individually
  # fusable below.
  active_members <- unlist(lapply(active_groups, `[[`, "members"))
  ordinary_props <- setdiff(names(observed_data), active_members)
  ordinary_props <- intersect(ordinary_props, sim_properties)

  for (prop in ordinary_props) {
    likelihood <- observed_data[[prop]]
    is_vector_likelihood <- is.atomic(likelihood) && is.numeric(likelihood) && is.null(names(likelihood))

    for (i in seq_len(n_horizons)) {
      prior <- simulation_params[[i]][[prop]]
      if (is.null(prior)) next

      posterior <- tryCatch(
        fuse_one_property_prior(prior, likelihood, is_vector_likelihood, n_samples,
                                 supported_closed_form, prop, i),
        error = function(e) {
          log_message("WARN", paste("Fusion failed for", prop, "at horizon", i, ":", conditionMessage(e)),
                      category = "MonteCarlo")
          NULL
        }
      )
      if (!is.null(posterior)) {
        simulation_params[[i]][[prop]] <- posterior
      }
    }
  }

  simulation_params
}

#' Reduce an `observed_data` entry to an unnamed `c(low, rep, high)` triplet
#'
#' Accepts a `c(low=,rep=,high=)`/`list(low=,rep=,high=)` triplet, an
#' unnamed length-3 `c(low,rep,high)` vector, or a raw sample vector (reduced
#' via its own empirical quantiles at `lh_probs`).
#' @keywords internal
as_lrh_triplet <- function(x, lh_probs = c(0.05, 0.95)) {
  nm <- names(x)
  if (!is.null(nm) && all(c("low", "rep", "high") %in% nm)) {
    return(c(x[["low"]], x[["rep"]], x[["high"]]))
  }
  if (!is.null(nm) && all(c("l", "r", "h") %in% nm)) {
    return(c(x[["l"]], x[["r"]], x[["h"]]))
  }
  if (length(x) == 3) {
    return(c(x[[1]], x[[2]], x[[3]]))
  }
  # Raw sample vector - reduce to a triplet via its own empirical quantiles.
  q <- unname(quantile(x, probs = c(lh_probs[1], 0.5, lh_probs[2]), na.rm = TRUE))
  q
}

#' Fuse one property's per-horizon prior against a shared observed-data likelihood
#'
#' @param prior A `simulation_params[[i]][[prop]]` entry: `list(family=, fit=, source=)`.
#' @param likelihood The `observed_data[[prop]]` entry (raw vector or param list).
#' @param is_vector_likelihood Precomputed `TRUE`/`FALSE` for `likelihood`'s shape.
#' @param n_samples Sample size for the general route.
#' @param supported_closed_form Character vector of families with a conjugate route.
#' @param prop,horizon_index Used only for WARN message context.
#' @return A new `list(family=, fit=, source=)` posterior, or `NULL` to skip
#'   (leaving the prior unchanged).
#' @keywords internal
fuse_one_property_prior <- function(prior, likelihood, is_vector_likelihood, n_samples,
                                     supported_closed_form, prop, horizon_index) {
  family <- prior$family
  fit <- prior$fit

  if (is_vector_likelihood) {
    # General route: fuse_property() requires BOTH sides to be raw sample
    # vectors and fits no parametric family to either - sample the prior first.
    prior_samples <- quantile_from_fit(runif(n_samples), family = family, fit = fit)
    posterior_samples <- fuse_property(prior_samples, likelihood)
    sorted_vals <- sort(posterior_samples)
    probs <- (seq_along(sorted_vals) - 0.5) / length(sorted_vals)
    return(list(
      family = "linear_cdf",
      fit = list(probs = probs, values = sorted_vals),
      source = "bayesian_fusion_general"
    ))
  }

  # Closed-form route: only normal/lognormal/beta priors have a conjugate route.
  if (!(family %in% supported_closed_form)) {
    log_message(
      "WARN",
      paste0("Skipping Bayesian fusion for '", prop, "' at horizon ", horizon_index,
             ": prior family '", family, "' has no closed-form fusion route; ",
             "simulating from the unfused prior."),
      category = "MonteCarlo"
    )
    return(NULL)
  }

  if (family == "normal") {
    post <- bayes_update_normal_normal(fit$mean, fit$sd, likelihood$mean, likelihood$sd)
    return(list(family = "normal", fit = list(mean = post$mu, sd = post$sigma), source = "bayesian_fusion_normal"))
  }

  if (family == "lognormal") {
    # Prior fit is ALREADY log-space (fit_percentile_triplet()'s lognormal
    # branch fits fit_normal_triplet(log(l), log(r), log(h), ...)) - do NOT
    # call normal_to_lognormal_params() on it, that would double-convert.
    # The likelihood is expected in raw (natural, reported) space and is
    # converted to log-space to match.
    lik_log <- normal_to_lognormal_params(likelihood$mean, likelihood$sd)
    post <- bayes_update_normal_normal(fit$mean, fit$sd, lik_log$mu, lik_log$sigma)
    return(list(family = "lognormal", fit = list(mean = post$mu, sd = post$sigma), source = "bayesian_fusion_lognormal"))
  }

  if (family == "beta") {
    lower <- fit$lower
    upper <- fit$upper
    span <- upper - lower
    if (!is.null(likelihood$shape1)) {
      lik_alpha <- likelihood$shape1
      lik_beta <- likelihood$shape2
    } else {
      # Raw-scale mean/var - rescale onto the PRIOR's own [lower, upper]
      # support before fitting: fuse_beta() simply adds alpha/beta, which is
      # only mathematically valid if both sides share the same support.
      scaled_mean <- (likelihood$mean - lower) / span
      scaled_var <- likelihood$var / span^2
      lik_moments <- moments_to_beta(scaled_mean, scaled_var)
      lik_alpha <- lik_moments$alpha
      lik_beta <- lik_moments$beta
    }

    post <- fuse_beta(fit$shape1, fit$shape2, lik_alpha, lik_beta)
    if (!isTRUE(post$feasible)) {
      # Documented fallback: re-express both sides as Normal moments, fuse,
      # convert back - rather than accepting invalid (alpha, beta).
      prior_m <- beta_to_moments(fit$shape1, fit$shape2)
      lik_m <- beta_to_moments(lik_alpha, lik_beta)
      fallback <- bayes_update_normal_normal(prior_m$mean, sqrt(prior_m$var), lik_m$mean, sqrt(lik_m$var))
      fallback_beta <- moments_to_beta(fallback$mu, fallback$sigma^2)
      post <- list(alpha = fallback_beta$alpha, beta = fallback_beta$beta)
    }
    return(list(
      family = "beta",
      fit = list(shape1 = post$alpha, shape2 = post$beta, lower = lower, upper = upper),
      source = "bayesian_fusion_beta"
    ))
  }

  NULL
}

#' Configure Correlation Structure (Enhanced)
#'
#' Enhanced correlation configuration
#'
#' @param simulation_params Simulation parameters
#' @param properties Property names
#' @param correlation_matrix Optional correlation matrix
#' @param config Configuration settings
#' @param simulation_data Optional raw per-horizon data frame (passed through
#'   to `estimate_property_correlations()` so genhz-stratified estimation can
#'   be used when a `genhz` column is present).
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Enhanced correlation configuration
#'
#' @export
configure_correlation_structure <- function(simulation_params,
                                                     properties,
                                                     correlation_matrix = NULL,
                                                     config,
                                                     simulation_data = NULL,
                                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Configuring enhanced correlation structure", category = "MonteCarlo")

  correlation_config <- list(
    matrix = NULL,
    method = "none",
    properties = properties,
    validation = list(),
    statistics = list()
  )

  n_properties <- length(properties)

  # Determine correlation approach
  if (!is.null(correlation_matrix)) {
    # Validate user-provided matrix
    matrix_validation <- validate_correlation_matrix(correlation_matrix, properties)

    if (matrix_validation$valid) {
      correlation_config$matrix <- correlation_matrix
      correlation_config$method <- "user_provided"
      log_message("INFO", "Using user-provided correlation matrix", category = "MonteCarlo")
    } else {
      log_message("WARN", paste("Invalid correlation matrix:", matrix_validation$message), category = "MonteCarlo")
      correlation_config$matrix <- diag(n_properties)
      correlation_config$method <- "identity_fallback"
    }

  } else if (config$monte_carlo$auto_correlation) {
    # Estimate correlation from data
    log_message("INFO", "Estimating correlation structure from data", category = "MonteCarlo")

    estimated_correlation <- estimate_property_correlations(
      simulation_params = simulation_params,
      properties = properties,
      config = config,
      simulation_data = simulation_data
    )

    if (estimated_correlation$valid) {
      correlation_config$matrix <- estimated_correlation$correlation_matrix
      correlation_config$method <- "estimated_from_data"
      correlation_config$estimation_info <- estimated_correlation$info
    } else {
      correlation_config$matrix <- diag(n_properties)
      correlation_config$method <- "identity_estimation_failed"
    }

  } else {
    # Use identity matrix (independent simulation)
    correlation_config$matrix <- diag(n_properties)
    correlation_config$method <- "identity_default"
    log_message("INFO", "Using independent simulation (identity matrix)", category = "MonteCarlo")
  }

  # Set property names
  rownames(correlation_config$matrix) <- properties
  colnames(correlation_config$matrix) <- properties

  # Ensure positive definiteness
  correlation_config$matrix <- ensure_positive_definite_matrix(correlation_config$matrix)

  # Calculate matrix statistics
  correlation_config$statistics <- calculate_matrix_statistics(correlation_config$matrix)

  # Final validation
  correlation_config$validation <- validate_correlation_matrix(
    correlation_config$matrix, properties
  )

  log_message("INFO", paste("Correlation method:", correlation_config$method), category = "MonteCarlo")
  log_message("DEBUG", paste("Condition number:",
                             round(correlation_config$statistics$condition_number, 2)), category = "MonteCarlo")

  return(correlation_config)
}

# ============================================================================
# 3. ENHANCED CONSTRAINT AND VALIDATION FUNCTIONS
# ============================================================================

#' Apply Simulation Constraints (Enhanced)
#'
#' Enhanced constraint application
#'
#' @param simulation_results Array of simulation results
#' @param properties Property names
#' @param config Configuration settings
#' @param composition_plan Optional result of `resolve_composition_groups()`;
#'   when a group is active, its sum-to-100 constraint is already exact by
#'   construction (via `restore_composition_properties()`'s ILR inverse), so
#'   `get_sum_constraints()` skips it rather than re-rescaling.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Constrained simulation results with quality metrics
#'
#' @export
apply_simulation_constraints <- function(simulation_results, properties, config, composition_plan = NULL,
                                          verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Applying enhanced simulation constraints", category = "MonteCarlo")

  constraint_rules <- get_constraint_rules(properties, config, composition_plan)
  constrained_results <- simulation_results
  constraint_summary <- list()

  # Apply range constraints  range validation
  if (!is.null(constraint_rules$range_constraints)) {
    log_message("DEBUG", "Applying range constraints", category = "MonteCarlo")

    range_result <- apply_range_constraints_batch(
      constrained_results, constraint_rules$range_constraints, properties
    )
    constrained_results <- range_result$data
    constraint_summary$range_applied <- range_result$summary
  }

  # Apply sum constraints (e.g., texture properties)
  if (!is.null(constraint_rules$sum_constraints)) {
    log_message("DEBUG", "Applying sum constraints", category = "MonteCarlo")

    sum_result <- apply_sum_constraints(
      constrained_results, constraint_rules$sum_constraints, properties
    )
    constrained_results <- sum_result$data
    constraint_summary$sum_applied <- sum_result$summary
  }

  # Apply relationship constraints
  if (!is.null(constraint_rules$relationship_constraints)) {
    log_message("DEBUG", "Applying relationship constraints", category = "MonteCarlo")

    relationship_result <- apply_relationship_constraints(
      constrained_results, constraint_rules$relationship_constraints, properties
    )
    constrained_results <- relationship_result$data
    constraint_summary$relationship_applied <- relationship_result$summary
  }

  # Apply physical plausibility constraints
  if (!is.null(constraint_rules$physical_constraints)) {
    log_message("DEBUG", "Applying physical constraints", category = "MonteCarlo")

    physical_result <- apply_physical_constraints(
      constrained_results, constraint_rules$physical_constraints, properties
    )
    constrained_results <- physical_result$data
    constraint_summary$physical_applied <- physical_result$summary
  }

  # Summary logging
  applied_constraints <- names(constraint_summary)
  log_message("INFO", paste("Applied constraints:", paste(applied_constraints, collapse = ", ")),
              category = "MonteCarlo")

  # Add constraint metadata to results
  attr(constrained_results, "constraint_summary") <- constraint_summary

  return(constrained_results)
}

#' Validate Monte Carlo Inputs (Enhanced)
#'
#' Comprehensive input validation
#'
#' @param soil_data Input soil data
#' @param properties Properties to simulate
#' @param correlation_matrix Optional correlation matrix
#' @param n_realizations Number of realizations
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results
#'
#' @export
validate_monte_carlo_inputs <- function(soil_data, properties, correlation_matrix,
                                        n_realizations, config,
                                        verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Validating Monte Carlo inputs", category = "MonteCarlo")

  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    data_validation = list(),
    property_validation = list(),
    parameter_validation = list()
  )

  # Basic parameter validation
  param_specs <- list(
    n_realizations = list(type = "numeric", required = TRUE, range = c(1, 100000))
  )

  param_validation <- validate_parameters(
    list(n_realizations = n_realizations),
    param_specs,
    strict_mode = TRUE
  )

  if (!param_validation$valid) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(validation_results$errors, param_validation$errors)
  }

  validation_results$parameter_validation <- param_validation

  # Data quality validation
  required_columns <- c("cokey")
  numeric_columns <- paste0(rep(properties, each = 3), c("_l", "_r", "_h"))

  # FIXED: Safe quality thresholds access
  quality_thresholds <- tryCatch({
    config$monte_carlo$quality_thresholds %||%
      config$quality_thresholds %||%
      list(completeness_threshold = 0.8, missing_data_threshold = 0.2)
  }, error = function(e) {
    list(completeness_threshold = 0.8, missing_data_threshold = 0.2)
  })

  data_validation <- validate_data_quality(
    soil_data,
    required_columns = required_columns,
    numeric_columns = intersect(numeric_columns, names(soil_data)),
    quality_thresholds = quality_thresholds
  )

  validation_results$data_validation <- data_validation

  # FIXED: Safe quality score access
  data_quality_score <- tryCatch({
    data_validation$overall_quality$score %||%
      data_validation$quality_score %||%
      0.8  # Default
  }, error = function(e) {
    0.8  # Fallback
  })

  min_quality_score <- tryCatch({
    config$monte_carlo$min_quality_score %||% 0.6
  }, error = function(e) {
    0.6
  })

  if (data_quality_score < min_quality_score) {
    validation_results$warnings <- c(
      validation_results$warnings,
      paste("Low data quality score:", round(data_quality_score, 3))
    )
  }

  # Property validation
  property_validation <- validate_properties_with_synonyms(
    properties,
    property_lookup = "ssurgo",
    strict_mode = FALSE
  )

  validation_results$property_validation <- property_validation

  if (!property_validation$valid) {
    strict_validation <- tryCatch({
      config$monte_carlo$strict_property_validation %||% FALSE
    }, error = function(e) {
      FALSE
    })

    if (strict_validation) {
      validation_results$valid <- FALSE
      validation_results$errors <- c(validation_results$errors, property_validation$errors)
    } else {
      validation_results$warnings <- c(validation_results$warnings, property_validation$warnings)
    }
  }

  # Correlation matrix validation if provided
  if (!is.null(correlation_matrix)) {
    correlation_validation <- validate_correlation_matrix(correlation_matrix, properties)

    if (!correlation_validation$valid) {
      validation_results$warnings <- c(
        validation_results$warnings,
        paste("Correlation matrix issue:", correlation_validation$message)
      )
    }
  }

  # FIXED: Check data sufficiency using the fixed function
  suitable_data_check <- check_data_sufficiency(soil_data, properties, config)
  if (!suitable_data_check$sufficient) {
    validation_results$warnings <- c(
      validation_results$warnings,
      paste("Limited data available for simulation:", suitable_data_check$message)
    )
  }

  log_message("DEBUG", paste("Input validation complete - Valid:", validation_results$valid,
                             "Errors:", length(validation_results$errors),
                             "Warnings:", length(validation_results$warnings)),
              category = "MonteCarlo")

  return(validation_results)
}

#' Validate Simulation Output (Enhanced)
#'
#' Enhanced output validation
#'
#' @param simulation_results Simulation results array
#' @param properties Property names
#' @param config Configuration settings
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Enhanced validation results
#'
#' @export
validate_simulation_output <- function(simulation_results, properties, config,
                                        verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("DEBUG", "Validating simulation output", category = "MonteCarlo")

  validation_results <- list(
    valid = TRUE,
    success_rate = 0,
    property_validation = list(),
    quality_metrics = list(),
    issues = character(0)
  )

  # Dimension validation
  dims <- dim(simulation_results)
  if (length(dims) != 3 || dims[2] != length(properties)) {
    validation_results$valid <- FALSE
    validation_results$issues <- c(validation_results$issues, "Invalid result dimensions")
    return(validation_results)
  }

  n_horizons <- dims[1]
  n_properties <- dims[2]
  n_realizations <- dims[3]
  total_values <- n_horizons * n_realizations

  # Property-specific validation
  for (i in seq_along(properties)) {
    prop <- properties[i]
    prop_values <- simulation_results[, i, ]

    # Use Module 0 utilities for statistical analysis
    prop_stats <- analyze_property_simulation_quality(prop_values, prop, config)
    validation_results$property_validation[[prop]] <- prop_stats

    # Check success rate
    if (prop_stats$success_rate < config$monte_carlo$min_success_rate) {
      validation_results$issues <- c(
        validation_results$issues,
        paste("Low success rate for", prop, ":", round(prop_stats$success_rate * 100, 1), "%")
      )
    }

    # Check for extreme values  outlier detection
    outliers <- detect_outliers(
      as.vector(prop_values),
      method = "iqr",
      threshold = config$monte_carlo$outlier_threshold
    )

    # mean(outliers, na.rm=TRUE) is NaN (not FALSE) when prop_values is
    # entirely non-finite (e.g. a property missing from every horizon) -
    # `if (NaN > x)` errors with "missing value where TRUE/FALSE needed"
    # rather than just skipping the check, which previously crashed the
    # whole pipeline on that edge case instead of degrading gracefully.
    outlier_rate <- mean(outliers, na.rm = TRUE)
    if (is.finite(outlier_rate) && outlier_rate > config$monte_carlo$max_outlier_rate) {
      validation_results$issues <- c(
        validation_results$issues,
        paste("High outlier rate for", prop, ":", round(outlier_rate * 100, 1), "%")
      )
    }
  }

  # Overall quality metrics
  validation_results$quality_metrics <- calculate_simulation_quality_metrics(
    simulation_results, properties, config
  )

  # Overall success rate
  total_finite <- sum(sapply(validation_results$property_validation, function(x) x$finite_values))
  validation_results$success_rate <- total_finite / (total_values * length(properties))

  # Final validation decision
  if (validation_results$success_rate < config$monte_carlo$min_success_rate) {
    validation_results$valid <- FALSE
  }

  log_message("INFO", paste("Output validation:", ifelse(validation_results$valid, "PASSED", "FAILED")),
              category = "MonteCarlo")
  log_message("INFO", paste("Overall success rate:", round(validation_results$success_rate * 100, 1), "%"),
              category = "MonteCarlo")

  if (length(validation_results$issues) > 0) {
    log_message("WARN", paste("Validation issues found:", length(validation_results$issues)),
                category = "MonteCarlo")
  }

  return(validation_results)
}

# ============================================================================
# 4. CONFIGURATION AND HELPER FUNCTIONS
# ============================================================================

#' Get Monte Carlo Default Configuration
#'
#' Returns default configuration optimized for Monte Carlo simulation
#'
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Default Monte Carlo configuration
#'
#' @export
get_monte_carlo_defaults <- function(verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  base_config <- tryCatch({
    get_default_configuration("full")
  }, error = function(e) {
    # Fallback if Module 0 function not available
    list()
  })

  # Complete Monte Carlo specific defaults
  monte_carlo_config <- list(
    # Basic simulation parameters
    distribution_type = "triangular",   # fallback family when a property has no property_distributions entry
    # What SSURGO's _l/_h are assumed to represent, as probabilities - a
    # documented, overridable assumption (SSURGO's low/high are professional-
    # judgment bounds, not confirmed exact percentiles) rather than a silently
    # hardcoded one.
    lh_percentile = c(0.05, 0.95),
    # Per-property distribution-family overrides, e.g.
    # list(claytotal = list(family = "beta", bounds = c(0, 100))) - resolved
    # by extract_property_parameters() ahead of the global distribution_type
    # fallback above.
    property_distributions = list(),
    # Properties that must be simulated JOINTLY (not independently) because
    # they're compositionally constrained (e.g. sand+silt+clay=100). See
    # resolve_composition_groups()/restore_composition_properties()
    # (distributions.R). `members` sets which real property occupies
    # `ilr_forward()`/`ilr_inverse()`'s position 1/2/3 (their `clay`/`sand`/
    # `silt` parameter names are positional-role placeholders, not an
    # identity requirement - see the header comment above those functions).
    # This default, (sandtotal, silttotal, claytotal), matches the
    # sequential binary partition `compositions::ilr()` used to build the
    # optional KSSL reference correlation matrices (`kssl-reference-correlations.R`)
    # - sand vs {silt,clay}, then silt vs clay - so that fallback feature's
    # ilr1/ilr2 correlations are directly reusable. This is a pure internal
    # reparameterization: simulated clay/sand/silt output is statistically
    # invariant to which role order is used (see test-distributions.R's
    # role-order invariance test).
    composition_groups = list(
      texture = list(members = c("sandtotal", "silttotal", "claytotal"), pseudo = c("ilr1", "ilr2"))
    ),
    max_depth = 250,
    auto_correlation = FALSE,
    # When auto_correlation estimation can't produce a matrix (too little
    # data), what to fall back to: "identity" (independent simulation,
    # today's exact behavior) or "kssl_global" (an opt-in fallback to the
    # KSSL reference correlation matrices - see kssl-reference-correlations.R
    # and estimate_property_correlations()). Defaults to "identity" - zero
    # behavior change unless a caller explicitly opts in.
    correlation_fallback = "identity",
    parallel_threshold = 1000,

    # Quality control
    min_quality_score = 0.6,
    min_success_rate = 0.8,
    outlier_threshold = 2.5,
    max_outlier_rate = 0.1,

    # FIXED: Add missing configuration parameters
    min_data_points = 5,           # Minimum data points per property
    minimum_observations = 10,      # Minimum observations for analysis

    # Error handling
    error_recovery_action = "warn",
    strict_property_validation = FALSE,
    strict_unsuitability = FALSE,

    # Constraints
    apply_range_constraints = TRUE,
    apply_sum_constraints = TRUE,
    apply_relationship_constraints = TRUE,
    apply_physical_constraints = FALSE,

    # Component simulation
    normalization_method = "proportional",
    composition_constraints = NULL,

    # Advanced options
    use_robust_methods = TRUE,
    preserve_correlations = TRUE,

    # Quality thresholds with safe defaults
    quality_thresholds = list(
      completeness_threshold = 0.8,
      missing_data_threshold = 0.2,
      outlier_threshold = 0.05,
      correlation_threshold = 0.05
    )
  )

  # Merge with base configuration safely
  if (length(base_config) > 0) {
    result_config <- merge_configurations_safe(base_config, list(monte_carlo = monte_carlo_config))
  } else {
    result_config <- list(monte_carlo = monte_carlo_config)
  }

  return(result_config)
}

#' Normalize a Monte Carlo `simulation_config` to be Properly Nested
#'
#' `generate_monte_carlo_realizations()`'s own `@examples` have always shown
#' a FLAT `simulation_config` (e.g. `list(distribution_type = "normal", max_depth = 200)`),
#' but `get_monte_carlo_defaults()` nests every Monte Carlo setting under
#' `$monte_carlo`, and `merge_configurations()` merges strictly by matching
#' key path - so a flat user config silently landed as a new, unused
#' top-level `config$distribution_type` and had NO effect at all (confirmed
#' by testing: `distribution_type = "normal"` passed flat left
#' `config$monte_carlo$distribution_type` at its default `"triangular"`,
#' every time). This accepts either shape: already-nested
#' (`list(monte_carlo = list(...))`) configs pass through unchanged; a flat
#' config whose names match known `monte_carlo` settings gets wrapped.
#'
#' @param simulation_config User-supplied config (flat or nested).
#' @param default_config Result of `get_monte_carlo_defaults()`, used only to
#'   recognize which flat key names belong under `monte_carlo`.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return `simulation_config`, wrapped under `monte_carlo` if it was flat.
#' @export
normalize_monte_carlo_config <- function(simulation_config, default_config,
                                          verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(simulation_config) || length(simulation_config) == 0) {
    return(simulation_config)
  }
  if ("monte_carlo" %in% names(simulation_config)) {
    return(simulation_config)
  }
  monte_carlo_keys <- names(default_config$monte_carlo)
  if (any(names(simulation_config) %in% monte_carlo_keys)) {
    return(list(monte_carlo = simulation_config))
  }
  simulation_config
}

#' Safe Configuration Merging
#'
#' Safe configuration merging that handles missing base configurations
#'
#' @param base_config Default/base configuration list.
#' @param new_config User-supplied override configuration list.
#' @return The merged configuration list.
merge_configurations_safe <- function(base_config, new_config) {
  tryCatch({
    if (exists("merge_configurations")) {
      return(merge_configurations(base_config, new_config))
    } else {
      # Simple merge fallback
      result <- base_config
      for (name in names(new_config)) {
        if (name %in% names(result) && is.list(result[[name]]) && is.list(new_config[[name]])) {
          # Merge nested lists
          for (sub_name in names(new_config[[name]])) {
            result[[name]][[sub_name]] <- new_config[[name]][[sub_name]]
          }
        } else {
          result[[name]] <- new_config[[name]]
        }
      }
      return(result)
    }
  }, error = function(e) {
    return(new_config)  # Return new config if merge fails
  })
}

#' Validate Monte Carlo Configuration
#'
#' Validate Monte Carlo configuration
#'
#' @param config Configuration to validate
#' @param n_realizations Number of realizations (for validation)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Validation results
#'
#' @export
validate_monte_carlo_config <- function(config, n_realizations, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Define parameter specifications
  param_specs <- list(
    distribution_type = list(type = "character", required = TRUE,
                             choices = c("triangular", "normal", "uniform", "beta",
                                         "lognormal", "metalog", "linear_cdf", "auto")),
    max_depth = list(type = "numeric", required = TRUE, range = c(10, 1000)),
    min_quality_score = list(type = "numeric", required = TRUE, range = c(0, 1)),
    min_success_rate = list(type = "numeric", required = TRUE, range = c(0, 1)),
    outlier_threshold = list(type = "numeric", required = TRUE, range = c(1, 5)),
    max_outlier_rate = list(type = "numeric", required = TRUE, range = c(0, 1)),
    error_recovery_action = list(type = "character", required = TRUE,
                                 choices = c("stop", "warn", "continue", "retry")),
    normalization_method = list(type = "character", required = TRUE,
                                choices = c("proportional", "additive", "multiplicative")),
    correlation_fallback = list(type = "character", required = TRUE,
                                choices = c("identity", "kssl_global"))
  )

  # Extract monte carlo config if nested
  if ("monte_carlo" %in% names(config)) {
    config_to_validate <- config$monte_carlo
  } else {
    config_to_validate <- config
  }

  # Use Module 0's validate_parameters
  validation <- validate_parameters(config_to_validate, param_specs, strict_mode = FALSE)

  # Additional Monte Carlo specific validations
  if (n_realizations > 50000 && (!is.null(config_to_validate$parallel) && !config_to_validate$parallel)) {
    validation$warnings <- c(validation$warnings,
                             "Large number of realizations without parallel processing may be slow")
  }

  return(validation)
}

# ============================================================================
# SUPPORTING HELPER FUNCTIONS (LEVERAGING Module 0)
# ============================================================================

#' Prepare Simulation Data (Enhanced)
#'
#' Enhanced data preparation
#'
#' @param soil_data Input soil data
#' @param properties Properties to prepare (may include composition-group
#'   pseudo-properties, e.g. "ilr1"/"ilr2" - see `composition_plan`)
#' @param config Configuration settings
#' @param composition_plan Optional result of `resolve_composition_groups()`;
#'   when supplied, data-availability checks are resolved back to the real
#'   underlying SSURGO columns (e.g. `sandtotal_r`/`claytotal_r`/`silttotal_r`)
#'   for any pseudo-property in `properties`, since a literal `"ilr1_r"`
#'   column never exists.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Prepared simulation data
#'
#' @export
prepare_simulation_data <- function(soil_data, properties, config, composition_plan = NULL,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Use Module 0's unsuitable horizon detection
  if (!"unsuitable_horizon" %in% names(soil_data)) {
    soil_data$unsuitable_horizon <- is_unsuitable(
      soil_data,
      strict_mode = config$monte_carlo$strict_unsuitability
    )
  }

  # Apply filters
  suitable_mask <- !soil_data$unsuitable_horizon

  # Apply depth filter if available
  if ("hzdepb_r" %in% names(soil_data)) {
    depth_mask <- soil_data$hzdepb_r <= config$monte_carlo$max_depth
    final_mask <- suitable_mask & depth_mask
  } else {
    final_mask <- suitable_mask
  }

  # Filter for minimum data availability, resolving any composition-group
  # pseudo-properties (e.g. "ilr1"/"ilr2") back to their real SSURGO columns first
  real_properties <- resolve_real_properties(properties, composition_plan)
  data_availability_mask <- check_property_data_availability(soil_data, real_properties, config)
  final_mask <- final_mask & data_availability_mask

  simulation_data <- soil_data[final_mask, ]

  log_message("DEBUG", paste("Data preparation: filtered from", nrow(soil_data), "to",
                             nrow(simulation_data), "horizons"), category = "MonteCarlo")

  return(simulation_data)
}

#' Resolve composition-group pseudo-properties to their real backing columns
#'
#' @param properties Character vector, possibly including pseudo-properties
#'   (e.g. "ilr1"/"ilr2") from an active composition group.
#' @param composition_plan Optional result of `resolve_composition_groups()`.
#' @return Character vector of real property names (pseudo-properties expanded
#'   to their group's `members`; everything else passed through unchanged).
#' @keywords internal
resolve_real_properties <- function(properties, composition_plan = NULL) {
  if (is.null(composition_plan)) {
    return(properties)
  }

  active_groups <- Filter(function(g) isTRUE(g$active), composition_plan$groups)
  real <- character(0)
  for (p in properties) {
    matched_group <- Find(function(g) p %in% g$pseudo, active_groups)
    if (!is.null(matched_group)) {
      real <- c(real, matched_group$members)
    } else {
      real <- c(real, p)
    }
  }
  unique(real)
}

# qtriangular()/transform_to_distribution()/ensure_positive_definite_matrix()/
# validate_correlation_matrix() previously lived here - they now live in
# distributions.R (as quantile_triangular(), quantile_from_fit(),
# ensure_positive_definite_matrix(), validate_correlation_matrix()) so
# statistics.R and bayesian-updating.R can share them too. No behavior
# change for the matrix utilities; transform_to_distribution() is replaced
# by direct quantile_from_fit() calls in simulate_correlated_properties().

# Sequential/parallel simulation dispatch (both real - despite the stale
# section name below, neither is a placeholder; both delegate to the real
# simulate_correlated_properties()).
run_sequential_simulation <- function(simulation_params, distribution_setup, correlation_config, n_realizations, properties, config) {
  return(simulate_correlated_properties(simulation_params, correlation_config$matrix, n_realizations, properties, config))
}

#' Run Monte Carlo simulation across multiple worker processes
#'
#' simulate_correlated_properties() generates all n_realizations for every
#' horizon in one call, with no state carried across horizons or
#' realizations, so it can safely be called once per worker with a smaller
#' n_realizations chunk and the resulting arrays concatenated along the
#' realization dimension. Follows the same Windows-PSOCK-vs-mclapply
#' pattern as mod07_multivariate_adjustment.R's process_cokeys_parallel(),
#' with a sequential fallback on any error.
#'
#' @param simulation_params List of per-horizon parameter lists.
#' @param distribution_setup Unused; kept for interface compatibility with `run_sequential_simulation()`.
#' @param correlation_config Result of `configure_correlation_structure()`.
#' @param n_realizations Number of realizations to generate.
#' @param properties Character vector of property names.
#' @param config Simulation configuration.
#' @return Same shape as simulate_correlated_properties(): an array
#'   dimensioned horizon by property by realization.
run_parallel_simulation <- function(simulation_params, distribution_setup, correlation_config, n_realizations, properties, config) {
  n_cores <- config$monte_carlo$n_cores %||% max(1, parallel::detectCores() - 1)
  n_cores <- min(n_cores, n_realizations)

  if (n_cores <= 1) {
    return(run_sequential_simulation(simulation_params, distribution_setup, correlation_config, n_realizations, properties, config))
  }

  # Split n_realizations into n_cores near-equal contiguous chunks.
  base_size <- n_realizations %/% n_cores
  remainder <- n_realizations %% n_cores
  chunk_sizes <- rep(base_size, n_cores) + c(rep(1, remainder), rep(0, n_cores - remainder))
  chunk_sizes <- chunk_sizes[chunk_sizes > 0]

  correlation_matrix <- correlation_config$matrix

  # future_seed = TRUE: simulate_correlated_properties() draws random values (rnorm() etc.), so
  # each chunk needs future's parallel-safe RNG streams - a genuine improvement over the prior
  # mclapply()/parLapply() path, which had no controlled cross-worker RNG stream at all.
  chunk_results <- run_parallel_lapply(
    chunk_sizes,
    function(chunk_n) simulate_correlated_properties(simulation_params, correlation_matrix, chunk_n, properties, config),
    n_cores = length(chunk_sizes),
    future_seed = TRUE,
    op_name = "Monte Carlo simulation",
    sequential_fallback = function() {
      run_sequential_simulation(simulation_params, distribution_setup, correlation_config, n_realizations, properties, config)
    }
  )

  # run_parallel_lapply()'s sequential_fallback already returns a full-shaped
  # [horizon, property, n_realizations] array (not a list of per-chunk arrays) when the parallel
  # path fails - detect and pass that straight through rather than trying to re-chunk-combine it.
  if (!is.list(chunk_results)) {
    return(chunk_results)
  }

  # Combine per-worker [horizon, property, chunk_n] arrays into one
  # [horizon, property, n_realizations] array by filling contiguous
  # realization-index slices (avoids adding an abind dependency).
  n_horizons <- length(simulation_params)
  n_properties <- length(properties)
  combined <- array(
    dim = c(n_horizons, n_properties, n_realizations),
    dimnames = list(horizon = 1:n_horizons, property = properties, realization = 1:n_realizations)
  )
  offset <- 0
  for (chunk in chunk_results) {
    chunk_n <- dim(chunk)[3]
    combined[, , (offset + 1):(offset + chunk_n)] <- chunk
    offset <- offset + chunk_n
  }

  combined
}

# Component composition helper functions
#
# Real prototype source (code/sim-functions.R:sim_component_comp()) confirms
# component_data's actual shape: comppct_l/comppct_r/comppct_h (component
# percent low/representative/high), with missing comppct_l/comppct_h filled
# from comppct_r -/+ 2. The functions below implement that shape for real
# (previously untagged placeholders returning hardcoded/pass-through values).

#' Validate Component Composition Data
#'
#' @param component_data Data frame expected to have a `comppct_r` column.
#' @return List with `valid`, `message`.
validate_component_data <- function(component_data) {
  if (!is.data.frame(component_data) || nrow(component_data) == 0) {
    return(list(valid = FALSE, message = "component_data must be a non-empty data frame"))
  }
  if (!"comppct_r" %in% names(component_data)) {
    return(list(valid = FALSE, message = "component_data must have a comppct_r column"))
  }
  if (all(is.na(component_data$comppct_r))) {
    return(list(valid = FALSE, message = "comppct_r is entirely missing"))
  }
  list(valid = TRUE, message = "Component data validation passed")
}

#' Extract Component Distribution Parameters
#'
#' Extracts `min`/`mode`/`max` for one component row from its
#' `comppct_l`/`comppct_r`/`comppct_h` triplet, filling missing
#' `comppct_l`/`comppct_h` from `comppct_r -/+ 2` (matching the legacy
#' `sim_component_comp()` fallback).
#'
#' @param component_row One-row data frame/list for a single component.
#' @param config Unused (kept for interface compatibility).
#' @return `list(min=, mode=, max=)`.
extract_component_parameters <- function(component_row, config) {
  comppct_r <- component_row$comppct_r

  if (is.null(comppct_r) || length(comppct_r) == 0 || is.na(comppct_r)) {
    return(list(min = 0, mode = 10, max = 20))
  }

  comppct_l <- component_row$comppct_l
  if (is.null(comppct_l) || length(comppct_l) == 0 || is.na(comppct_l)) {
    comppct_l <- comppct_r - 2
  }

  comppct_h <- component_row$comppct_h
  if (is.null(comppct_h) || length(comppct_h) == 0 || is.na(comppct_h)) {
    comppct_h <- comppct_r + 2
  }

  list(min = max(0, comppct_l), mode = comppct_r, max = comppct_h)
}

generate_component_values <- function(parameters, n_realizations, distribution_type) {
  return(runif(n_realizations, parameters$min, parameters$max))  # Placeholder
}

normalize_component_realizations <- function(component_realizations, target_sum, method) {
  return(apply(component_realizations, 2, function(x) x / sum(x) * target_sum))
}

#' Apply Composition Constraints
#'
#' Clamps `normalized_realizations` (an `[n_components, n_realizations]`
#' matrix) to `constraints$min`/`constraints$max` when supplied.
#'
#' @param normalized_realizations Matrix of normalized component realizations.
#' @param constraints List, optionally with `min`/`max`, or `NULL`.
#' @return Constrained matrix (same dimensions as input).
apply_composition_constraints <- function(normalized_realizations, constraints) {
  if (is.null(constraints) || (is.null(constraints$min) && is.null(constraints$max))) {
    return(normalized_realizations)
  }

  lo <- constraints$min %||% -Inf
  hi <- constraints$max %||% Inf

  constrained <- pmin(pmax(normalized_realizations, lo), hi)
  dim(constrained) <- dim(normalized_realizations)
  dimnames(constrained) <- dimnames(normalized_realizations)
  constrained
}

#' Assess Component Composition Quality
#'
#' Compares each component's mean simulated value against its original
#' `comppct_r`, as a relative-error-based quality score.
#'
#' @param original_data Original component data with a `comppct_r` column.
#' @param simulated_data `[n_components, n_realizations]` matrix of
#'   simulated/normalized component values.
#' @param config Unused (kept for interface compatibility).
#' @return List with `overall_quality`, `mean_relative_error`.
assess_component_quality <- function(original_data, simulated_data, config) {
  if (is.null(original_data) || !("comppct_r" %in% names(original_data)) ||
      is.null(simulated_data) || nrow(original_data) != nrow(simulated_data)) {
    return(list(overall_quality = NA_real_, mean_relative_error = NA_real_))
  }

  sim_means <- rowMeans(simulated_data, na.rm = TRUE)
  orig_values <- original_data$comppct_r
  denom <- pmax(abs(orig_values), 1e-8)
  relative_errors <- abs(sim_means - orig_values) / denom

  list(
    overall_quality = max(0, 1 - mean(relative_errors, na.rm = TRUE)),
    mean_relative_error = mean(relative_errors, na.rm = TRUE)
  )
}

#' Extract Simulation Parameters for One Horizon/Property (Enhanced)
#'
#' Resolves the property's distribution family - from
#' `config$monte_carlo$property_distributions[[property_name]]$family`, else
#' `config$monte_carlo$distribution_type` (default `"triangular"`) - and fits
#' family-appropriate parameters from its SSURGO `_l/_r/_h` triplet via
#' `distributions.R`'s `fit_percentile_triplet()`. Previously this only ever
#' returned `list(min=, mode=, max=)` regardless of the requested
#' `distribution_type`, so `"normal"`/`"beta"` silently produced `NA` values
#' downstream (their shape parameters were never actually computed).
#'
#' @param horizon_data One-row data frame/list for a single horizon.
#' @param property_name Property name (without `_l`/`_r`/`_h` suffix).
#' @param config Simulation configuration.
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return `list(valid=, parameters=list(family=, fit=, source=))`.
#' @export
extract_property_parameters <- function(horizon_data, property_name, config,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  l_col <- paste0(property_name, "_l")
  r_col <- paste0(property_name, "_r")
  h_col <- paste0(property_name, "_h")

  l_val <- horizon_data[[l_col]]
  r_val <- horizon_data[[r_col]]
  h_val <- horizon_data[[h_col]]

  registry_entry <- config$monte_carlo$property_distributions[[property_name]]
  family <- registry_entry$family %||% config$monte_carlo$distribution_type %||% "triangular"
  lh_probs <- registry_entry$lh_percentile %||% config$monte_carlo$lh_percentile %||% c(0.05, 0.95)
  bounds <- registry_entry$bounds

  if (all(is.finite(c(l_val, r_val, h_val)))) {
    fitted <- fit_percentile_triplet(l_val, r_val, h_val, family, lh_probs = lh_probs, bounds = bounds)
    return(list(
      valid = TRUE,
      parameters = list(family = fitted$family, fit = fitted$fit, source = "complete")
    ))
  } else if (is.finite(r_val)) {
    # Insufficient data to fit the requested family - fall back to a simple
    # +/-20% triangular estimate around the representative value, but still
    # tag family="triangular" explicitly so quantile_from_fit() dispatches
    # correctly regardless of the global/registry setting.
    range_est <- r_val * 0.2
    return(list(
      valid = TRUE,
      parameters = list(
        family = "triangular",
        fit = list(min = max(0, r_val - range_est), mode = r_val, max = r_val + range_est),
        source = "estimated"
      )
    ))
  } else {
    return(list(valid = FALSE, parameters = NULL))
  }
}

summarize_parameter_extraction <- function(simulation_params, properties, n_horizons) {
  summary <- list()
  for (prop in properties) {
    complete <- sum(sapply(simulation_params, function(x) prop %in% names(x) && x[[prop]]$source == "complete"))
    # "estimated" also counts composition-group pseudo-properties fit via
    # estimate_ilr_moments_mc() (source == "ilr_mc") - otherwise they'd be
    # misreported as 100% missing despite always being present when their
    # group is active.
    estimated <- sum(sapply(simulation_params, function(x) prop %in% names(x) && x[[prop]]$source %in% c("estimated", "ilr_mc")))
    missing <- n_horizons - complete - estimated

    summary[[prop]] <- list(
      complete = complete,
      estimated = estimated,
      missing = missing,
      coverage_rate = (complete + estimated) / n_horizons
    )
  }
  return(summary)
}

#' Estimate a Correlation Matrix from Simulation Parameters (Real Implementation)
#'
#' Builds one representative-value-per-horizon data frame - the per-horizon
#' fitted mean/mode for each property (for composition-group pseudo-properties
#' like ilr1/ilr2 there's no raw `_r` column, so the fitted mean IS the
#' representative value) - plus `genhz` from `simulation_data` if present,
#' then calls `distributions.R`'s `estimate_correlation_matrix_robust()`
#' (`Hmisc::rcorr()` + PD-repair, genhz-stratified when possible). Previously
#' this was a placeholder always returning `valid=FALSE`, so
#' `auto_correlation=TRUE` silently fell back to the identity matrix.
#'
#' When `config$monte_carlo$correlation_fallback == "kssl_global"` (opt-in;
#' defaults to `"identity"`, today's exact behavior), a group that fails
#' empirical estimation falls back to that group's own KSSL reference
#' correlation matrix (`kssl-reference-correlations.R`) instead of being
#' dropped, and the final identity fallback becomes a KSSL-pooled matrix
#' instead. `genhz` is auto-derived from `simulation_data$hzname` via
#' `classify_genhz()` when this is requested and `simulation_data$genhz`
#' isn't already present - scoped narrowly behind the opt-in flag so a
#' caller who has `hzname` but never requests `"kssl_global"` sees no
#' behavior change at all.
#'
#' @param simulation_params List of per-horizon parameter lists.
#' @param properties Character vector of property names.
#' @param config Simulation configuration.
#' @param simulation_data Optional raw per-horizon data frame (same row order
#'   as `simulation_params`); its `genhz` column, if present, enables
#'   genhz-stratified estimation (an explicit `genhz` column always takes
#'   precedence over auto-derivation from `hzname`).
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's
#'   log level so \code{INFO}-level progress messages print for the duration
#'   of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return `list(valid=, correlation_matrix=, info=list(method=, n_obs=))`.
#' @export
estimate_property_correlations <- function(simulation_params, properties, config, simulation_data = NULL,
                                            verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  representative_value <- function(params) {
    if (is.null(params)) return(NA_real_)
    fit <- params$fit
    switch(params$family %||% "triangular",
      triangular = fit$mode,
      uniform = mean(c(fit$min, fit$max)),
      normal = fit$mean,
      lognormal = exp(fit$mean),
      beta = fit$lower + (fit$shape1 / (fit$shape1 + fit$shape2)) * (fit$upper - fit$lower),
      metalog = fit$fallback_values[match(0.5, fit$fallback_probs)],
      linear_cdf = fit$values[match(0.5, fit$probs)],
      NA_real_
    )
  }

  rep_matrix <- do.call(rbind, lapply(simulation_params, function(horizon_params) {
    vapply(properties, function(p) representative_value(horizon_params[[p]]) %||% NA_real_, numeric(1))
  }))
  rep_df <- as.data.frame(rep_matrix)
  names(rep_df) <- properties

  n_properties <- length(properties)
  global_fallback <- diag(n_properties)
  rownames(global_fallback) <- colnames(global_fallback) <- properties
  global_fallback_method <- "global_fallback"

  use_kssl_fallback <- isTRUE(config$monte_carlo$correlation_fallback == "kssl_global")

  group_var <- NULL
  if (!is.null(simulation_data) && "genhz" %in% names(simulation_data) && nrow(simulation_data) == nrow(rep_df)) {
    rep_df$genhz <- simulation_data$genhz
    group_var <- "genhz"
  } else if (use_kssl_fallback && !is.null(simulation_data) &&
             "hzname" %in% names(simulation_data) && nrow(simulation_data) == nrow(rep_df)) {
    rep_df$genhz <- classify_genhz(simulation_data$hzname)
    group_var <- "genhz"
  }

  group_fallback_matrices <- NULL
  if (use_kssl_fallback) {
    kssl_pooled <- build_kssl_fallback_matrix(properties, genhz = NULL)
    if (!is.null(kssl_pooled)) {
      global_fallback <- kssl_pooled
      global_fallback_method <- "kssl_global_fallback"
    }
    if (!is.null(group_var)) {
      present_genhz <- unique(rep_df$genhz[!is.na(rep_df$genhz)])
      group_fallback_matrices <- stats::setNames(
        lapply(present_genhz, function(g) build_kssl_fallback_matrix(properties, genhz = g)),
        as.character(present_genhz)
      )
      group_fallback_matrices <- Filter(Negate(is.null), group_fallback_matrices)
    }
  }

  result <- estimate_correlation_matrix_robust(
    rep_df, group_var = group_var, min_group_n = 5, global_fallback = global_fallback,
    group_fallback_matrices = group_fallback_matrices, global_fallback_method = global_fallback_method
  )
  # estimate_correlation_matrix_robust() orders columns however rep_df's
  # numeric columns happened to be found - reorder/rename to guarantee the
  # matrix matches `properties`' exact order for downstream Cholesky use.
  ordered_matrix <- result$matrix[properties, properties, drop = FALSE]

  list(
    valid = TRUE,
    correlation_matrix = ordered_matrix,
    info = list(method = result$method, n_obs = result$n_obs)
  )
}

calculate_matrix_statistics <- function(matrix) {
  eigenvals <- eigen(matrix, only.values = TRUE)$values
  return(list(
    eigenvalues = eigenvals,
    condition_number = max(eigenvals) / min(eigenvals[eigenvals > 1e-12]),
    determinant = det(matrix)
  ))
}

get_constraint_rules <- function(properties, config, composition_plan = NULL) {
  return(list(
    range_constraints = get_range_constraints(properties),
    sum_constraints = get_sum_constraints(properties, composition_plan),
    relationship_constraints = get_relationship_constraints(properties),
    physical_constraints = if(config$monte_carlo$apply_physical_constraints) get_physical_constraints(properties) else NULL
  ))
}

get_range_constraints <- function(properties) {
  constraints <- list()
  for (prop in properties) {
    if (prop %in% c("sandtotal", "claytotal", "silttotal")) {
      constraints[[prop]] <- c(0, 100)
    } else if (prop == "dbovendry") {
      constraints[[prop]] <- c(0.3, 3.0)
    } else if (prop %in% c("ph1to1h2o", "phh2o")) {
      constraints[[prop]] <- c(2.5, 11.0)
    }
  }
  return(constraints)
}

#' Determine sum-to-100 constraints (e.g. texture) to apply post-simulation
#'
#' @param properties Character vector of property names.
#' @param composition_plan Optional result of `resolve_composition_groups()`;
#'   when the texture group is active, sum-to-100 is already exact by
#'   construction (via `restore_composition_properties()`'s ILR inverse), so
#'   no proportional-rescale entry is returned for it - re-rescaling would
#'   distort the correlation structure the ILR path exists to preserve.
get_sum_constraints <- function(properties, composition_plan = NULL) {
  if (isTRUE(composition_plan$groups$texture$active)) {
    return(NULL)
  }

  texture_props <- intersect(properties, c("sandtotal", "claytotal", "silttotal"))
  if (length(texture_props) >= 2) {
    return(list(texture = list(properties = texture_props, target_sum = 100, tolerance = 5)))
  }
  return(NULL)
}

get_relationship_constraints <- function(properties) {
  water_props <- intersect(properties, c("wthirdbar", "wfifteenbar"))
  if (length(water_props) == 2) {
    return(list(water_retention = list(lower_property = "wfifteenbar", upper_property = "wthirdbar")))
  }
  return(NULL)
}

get_physical_constraints <- function(properties) {
  # Real, soil-science-grounded plausible-range bounds (same source values as
  # data-infilling.R's apply_basic_range_limits()). Only
  # includes properties NOT already covered by get_range_constraints() above,
  # to avoid duplicating the same clamp under two different constraint
  # categories.
  physical_limits <- list(
    partdensity = c(2.0, 3.2),
    cec7 = c(0, 200),
    cec82 = c(0, 200),
    ecec = c(0, 200),
    om = c(0, 100),
    oc = c(0, 60),
    rfv = c(0, 95),
    fragvol = c(0, 95),
    wthirdbar = c(0, 70),
    wfifteenbar = c(0, 60),
    awc = c(0, 0.5)
  )

  constraints <- list()
  for (prop in intersect(properties, names(physical_limits))) {
    constraints[[prop]] <- physical_limits[[prop]]
  }
  return(constraints)
}

# Additional helper functions for enhanced functionality
apply_range_constraints_batch <- function(results, constraints, properties) {
  adjustments <- 0
  for (i in seq_along(properties)) {
    prop <- properties[i]
    if (prop %in% names(constraints)) {
      range_vals <- constraints[[prop]]
      before <- results[, i, ]
      results[, i, ] <- pmax(results[, i, ], range_vals[1])
      results[, i, ] <- pmin(results[, i, ], range_vals[2])
      adjustments <- adjustments + sum(before != results[, i, ], na.rm = TRUE)
    }
  }
  return(list(data = results, summary = list(adjustments = adjustments)))
}

apply_sum_constraints <- function(results, constraints, properties) {
  adjustments <- 0
  for (constraint_name in names(constraints)) {
    constraint <- constraints[[constraint_name]]
    prop_indices <- match(constraint$properties, properties)
    prop_indices <- prop_indices[!is.na(prop_indices)]

    if (length(prop_indices) >= 2) {
      for (h in seq_len(dim(results)[1])) {
        for (r in seq_len(dim(results)[3])) {
          current_sum <- sum(results[h, prop_indices, r])
          if (current_sum > 0 && abs(current_sum - constraint$target_sum) > 0.01) {
            results[h, prop_indices, r] <- results[h, prop_indices, r] * constraint$target_sum / current_sum
            adjustments <- adjustments + 1
          }
        }
      }
    }
  }
  return(list(data = results, summary = list(adjustments = adjustments)))
}

apply_relationship_constraints <- function(results, constraints, properties) {
  adjustments <- 0
  for (constraint_name in names(constraints)) {
    constraint <- constraints[[constraint_name]]
    lower_idx <- match(constraint$lower_property, properties)
    upper_idx <- match(constraint$upper_property, properties)

    if (!is.na(lower_idx) && !is.na(upper_idx)) {
      lower_vals <- results[, lower_idx, ]
      upper_vals <- results[, upper_idx, ]
      violated <- !is.na(lower_vals) & !is.na(upper_vals) & (lower_vals > upper_vals)

      if (any(violated)) {
        # Clip the lower-bound property down to the upper-bound property's
        # value wherever the ordering is violated (e.g. wfifteenbar <= wthirdbar).
        lower_vals[violated] <- upper_vals[violated]
        results[, lower_idx, ] <- lower_vals
        adjustments <- adjustments + sum(violated)
      }
    }
  }
  return(list(data = results, summary = list(adjustments = adjustments)))
}

apply_physical_constraints <- function(results, constraints, properties) {
  # Same clamp logic as apply_range_constraints_batch(), applied against the
  # physical-bounds table from get_physical_constraints() instead of the
  # general range-constraints table.
  adjustments <- 0
  for (i in seq_along(properties)) {
    prop <- properties[i]
    if (prop %in% names(constraints)) {
      range_vals <- constraints[[prop]]
      before <- results[, i, ]
      results[, i, ] <- pmax(results[, i, ], range_vals[1])
      results[, i, ] <- pmin(results[, i, ], range_vals[2])
      adjustments <- adjustments + sum(before != results[, i, ], na.rm = TRUE)
    }
  }
  return(list(data = results, summary = list(adjustments = adjustments)))
}

#' Check Data Sufficiency (FIXED VERSION)
#'
#' Fixed version with safe configuration access and robust error handling
#'
#' @param soil_data Input soil data
#' @param properties Properties to check
#' @param config Configuration object
#' @return Data sufficiency results
#'
check_data_sufficiency <- function(soil_data, properties, config) {

  # FIXED: Safe configuration access with defaults
  min_data_points <- tryCatch({
    config$monte_carlo$min_data_points %||%
      config$min_data_points %||%
      5  # Safe default
  }, error = function(e) {
    5  # Fallback default
  })

  sufficient_properties <- 0
  property_details <- list()

  for (prop in properties) {
    r_col <- paste0(prop, "_r")

    if (r_col %in% names(soil_data)) {
      # FIXED: Robust data availability check
      available_data <- tryCatch({
        sum(!is.na(soil_data[[r_col]]) & is.finite(soil_data[[r_col]]))
      }, error = function(e) {
        0  # Return 0 if calculation fails
      })

      property_details[[prop]] <- list(
        column_exists = TRUE,
        available_data = available_data,
        sufficient = available_data >= min_data_points
      )

      if (available_data >= min_data_points) {
        sufficient_properties <- sufficient_properties + 1
      }
    } else {
      property_details[[prop]] <- list(
        column_exists = FALSE,
        available_data = 0,
        sufficient = FALSE
      )
    }
  }

  # Overall sufficiency assessment
  total_properties <- length(properties)
  required_sufficient <- max(1, floor(total_properties * 0.5))  # At least 50% of properties
  sufficient <- sufficient_properties >= required_sufficient

  return(list(
    sufficient = sufficient,
    message = paste(sufficient_properties, "of", total_properties, "properties have sufficient data"),
    details = list(
      sufficient_properties = sufficient_properties,
      total_properties = total_properties,
      min_data_points = min_data_points,
      required_sufficient = required_sufficient,
      property_details = property_details
    )
  ))
}

check_property_data_availability <- function(soil_data, properties, config) {
  # Check if each row has data for at least one property
  has_data <- rep(FALSE, nrow(soil_data))

  for (i in seq_len(nrow(soil_data))) {
    for (prop in properties) {
      r_col <- paste0(prop, "_r")
      if (r_col %in% names(soil_data) && !is.na(soil_data[[r_col]][i])) {
        has_data[i] <- TRUE
        break
      }
    }
  }

  return(has_data)
}

analyze_property_simulation_quality <- function(prop_values, prop_name, config) {
  finite_values <- sum(is.finite(prop_values))
  total_values <- length(prop_values)

  # range(x, na.rm=TRUE) on an entirely non-finite vector (e.g. a property
  # missing from every horizon) warns "no non-missing arguments to min/max"
  # and returns c(Inf, -Inf) - report NA instead, which is the honest value
  # and avoids that spurious warning.
  value_range <- if (finite_values == 0) c(NA_real_, NA_real_) else range(prop_values, na.rm = TRUE)

  return(list(
    finite_values = finite_values,
    total_values = total_values,
    success_rate = finite_values / total_values,
    mean_value = mean(prop_values, na.rm = TRUE),
    sd_value = sd(prop_values, na.rm = TRUE),
    range = value_range
  ))
}

#' Calculate Simulation Quality Metrics
#'
#' Computes real per-property finite-value and outlier rates from the
#' `[horizon, property, realization]` simulation array, and an overall
#' quality score combining them.
#'
#' @param simulation_results `[horizon, property, realization]` array.
#' @param properties Character vector of property names (column order of
#'   `simulation_results`'s second dimension).
#' @param config Simulation configuration (uses `config$monte_carlo$outlier_threshold`).
#' @return List with `per_property` (named list of `finite_rate`/`outlier_rate`),
#'   `overall_finite_rate`, `overall_outlier_rate`, `overall_quality`.
calculate_simulation_quality_metrics <- function(simulation_results, properties, config) {
  per_property <- stats::setNames(vector("list", length(properties)), properties)

  for (i in seq_along(properties)) {
    values <- simulation_results[, i, ]
    finite_rate <- mean(is.finite(values))

    outliers <- tryCatch(
      detect_outliers(as.vector(values), method = "iqr", threshold = config$monte_carlo$outlier_threshold),
      error = function(e) rep(FALSE, length(values))
    )
    outlier_rate <- mean(outliers, na.rm = TRUE)

    per_property[[properties[i]]] <- list(
      finite_rate = finite_rate,
      outlier_rate = if (is.finite(outlier_rate)) outlier_rate else NA_real_
    )
  }

  overall_finite_rate <- mean(vapply(per_property, function(x) x$finite_rate, numeric(1)))
  overall_outlier_rate <- mean(vapply(per_property, function(x) x$outlier_rate %||% 0, numeric(1)), na.rm = TRUE)
  overall_quality <- max(0, min(1, overall_finite_rate - overall_outlier_rate))

  list(
    per_property = per_property,
    overall_finite_rate = overall_finite_rate,
    overall_outlier_rate = overall_outlier_rate,
    overall_quality = overall_quality
  )
}

#' Generate Simulation Diagnostics
#'
#' Computes real per-property simulated vs. original mean/SD comparisons.
#'
#' @param original_data Original SSURGO-derived data.
#' @param simulation_results Constrained `[horizon, property, realization]` array.
#' @param properties Character vector of property names.
#' @param correlation_config Correlation configuration used for simulation.
#' @param config Unused (kept for interface compatibility).
#' @return List with `per_property` (named list of `sim_mean`/`sim_sd`/
#'   `orig_mean`/`orig_sd`), `correlation_method`, `n_properties`, `timestamp`.
generate_simulation_diagnostics <- function(original_data, simulation_results, properties, correlation_config, config) {
  per_property <- stats::setNames(vector("list", length(properties)), properties)

  for (i in seq_along(properties)) {
    prop <- properties[i]
    sim_values <- as.vector(simulation_results[, i, ])
    sim_values <- sim_values[is.finite(sim_values)]

    orig_values <- if (prop %in% names(original_data)) {
      original_data[[prop]][is.finite(original_data[[prop]])]
    } else {
      numeric(0)
    }

    per_property[[prop]] <- list(
      sim_mean = if (length(sim_values) > 0) mean(sim_values) else NA_real_,
      sim_sd = if (length(sim_values) > 1) stats::sd(sim_values) else NA_real_,
      orig_mean = if (length(orig_values) > 0) mean(orig_values) else NA_real_,
      orig_sd = if (length(orig_values) > 1) stats::sd(orig_values) else NA_real_
    )
  }

  list(
    per_property = per_property,
    correlation_method = correlation_config$method %||% NA_character_,
    n_properties = length(properties),
    timestamp = Sys.time()
  )
}

#' Assess Simulation Quality
#'
#' Derives real constraint-satisfaction and statistical-consistency scores
#' from `diagnostics$per_property` (simulated vs. original mean/SD), and
#' combines them with `output_validation$success_rate` into an overall
#' quality score.
#'
#' @param simulation_results Unused (kept for interface compatibility).
#' @param diagnostics Result of [generate_simulation_diagnostics()].
#' @param output_validation Result of `validate_simulation_output()`, with
#'   a `success_rate` field.
#' @param config Unused (kept for interface compatibility).
#' @return List with `overall_quality_score`, `component_scores`
#'   (`data_quality`/`constraint_satisfaction`/`statistical_consistency`),
#'   `recommendations`.
assess_simulation_quality <- function(simulation_results, diagnostics, output_validation, config) {
  data_quality <- output_validation$success_rate %||% NA_real_

  per_property <- diagnostics$per_property %||% list()

  constraint_scores <- vapply(per_property, function(p) {
    if (is.null(p$orig_mean) || !is.finite(p$orig_mean) || !is.finite(p$sim_mean)) return(NA_real_)
    denom <- max(abs(p$orig_mean), 1e-8)
    max(0, 1 - abs(p$sim_mean - p$orig_mean) / denom)
  }, numeric(1))
  constraint_satisfaction <- if (any(is.finite(constraint_scores))) mean(constraint_scores, na.rm = TRUE) else NA_real_

  consistency_scores <- vapply(per_property, function(p) {
    if (is.null(p$orig_sd) || !is.finite(p$orig_sd) || !is.finite(p$sim_sd) || p$orig_sd <= 0) return(NA_real_)
    max(0, 1 - abs(p$sim_sd - p$orig_sd) / p$orig_sd)
  }, numeric(1))
  statistical_consistency <- if (any(is.finite(consistency_scores))) mean(consistency_scores, na.rm = TRUE) else NA_real_

  component_scores <- list(
    data_quality = data_quality,
    constraint_satisfaction = constraint_satisfaction,
    statistical_consistency = statistical_consistency
  )

  valid_scores <- unlist(component_scores)
  valid_scores <- valid_scores[is.finite(valid_scores)]
  overall_quality_score <- if (length(valid_scores) > 0) mean(valid_scores) else NA_real_

  recommendations <- if (is.finite(overall_quality_score) && overall_quality_score >= 0.8) {
    "Simulation completed successfully"
  } else {
    "Simulation quality below target - review property_validation/diagnostics for details"
  }

  list(
    overall_quality_score = overall_quality_score,
    component_scores = component_scores,
    recommendations = recommendations
  )
}

#' Validate a Horizon's Fitted Distribution Parameters (Real Implementation)
#'
#' Thin wrapper around `distributions.R`'s `validate_fit_parameters()`.
#' Previously this was a stub always returning `valid=TRUE` regardless of
#' input, so nothing ever caught a malformed/undefined parameter set.
#'
#' @param parameters A params object from `extract_property_parameters()`
#'   (`list(family=, fit=, source=)`).
#' @param distribution_type Fallback family if `parameters$family` is absent.
#' @param property_name Unused; kept for interface compatibility with callers.
#' @return `list(valid=, message=)`.
validate_distribution_parameters <- function(parameters, distribution_type, property_name) {
  family <- parameters$family %||% distribution_type
  fit <- parameters$fit %||% parameters
  validate_fit_parameters(family, fit)
}

create_distribution_config <- function(parameters, distribution_type, property_name, validation) {
  return(list(
    parameters = parameters,
    family = parameters$family %||% distribution_type,
    property_name = property_name,
    valid = validation$valid
  ))
}

create_fallback_distribution <- function(property_name, representative_value, distribution_type) {
  range_est <- representative_value * 0.2
  return(list(
    parameters = list(
      family = "triangular",
      fit = list(min = max(0, representative_value - range_est),
                 mode = representative_value,
                 max = representative_value + range_est)
    ),
    family = "triangular",
    property_name = property_name,
    valid = FALSE,
    source = "fallback"
  ))
}

#' Summarize Distribution Setup
#'
#' Real summary of `distributions` (a per-horizon list of per-property
#' distribution configs from `create_distribution_config()`/
#' `create_fallback_distribution()`, both of which carry a `valid` and a
#' `family` field): counts and per-family breakdown.
#'
#' @param distributions Per-horizon list of per-property distribution configs.
#' @param properties Unused (kept for interface compatibility).
#' @param distribution_type Fallback family label for entries with no `family`.
#' @return List with `total_distributions`, `valid_distributions`,
#'   `distribution_type`, `family_counts`.
summarize_distributions <- function(distributions, properties, distribution_type) {
  all_entries <- unlist(distributions, recursive = FALSE)
  total <- length(all_entries)
  valid <- sum(vapply(all_entries, function(x) isTRUE(x$valid), logical(1)))
  families <- table(vapply(all_entries, function(x) x$family %||% distribution_type, character(1)))

  list(
    total_distributions = total,
    valid_distributions = valid,
    distribution_type = distribution_type,
    family_counts = as.list(families)
  )
}

#' Validate Distribution Setup
#'
#' Real validity-rate computation over `distributions` (see
#' [summarize_distributions()] for the expected shape).
#'
#' @param distributions Per-horizon list of per-property distribution configs.
#' @param properties Unused (kept for interface compatibility).
#' @param config Unused (kept for interface compatibility).
#' @return List with `valid_distributions`, `total_distributions`, `validity_rate`.
validate_distribution_setup <- function(distributions, properties, config) {
  all_entries <- unlist(distributions, recursive = FALSE)
  total <- length(all_entries)
  valid <- sum(vapply(all_entries, function(x) isTRUE(x$valid), logical(1)))

  list(
    valid_distributions = valid,
    total_distributions = total,
    validity_rate = if (total > 0) valid / total else NA_real_
  )
}

validate_simulated_array <- function(simulated_values, properties, config) {
  finite_count <- sum(is.finite(simulated_values))
  total_count <- length(simulated_values)

  return(list(
    finite_rate = finite_count / total_count,
    total_values = total_count,
    finite_values = finite_count
  ))
}
