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

  # PERF: process_single_cokey() previously re-filtered the FULL multi-cokey simulation_data
  # (dplyr::filter(cokey == !!cokey)) once per cokey - O(rows x cokeys) instead of O(rows). Split
  # once here instead (mirrors the split()-based fix already applied to
  # maybe_adjust_soil_data_depth_trend()/run_stage1_fusion_group() elsewhere in this package), and
  # pass each cokey's own pre-split subset down - see PERFORMANCE_IMPROVEMENT_PLAN.md Tier 1.
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
#' @param config Optional Monte Carlo config (as from \code{get_monte_carlo_defaults()}) whose
#'   \code{monte_carlo$vertical_correlation_method} selects between \code{"joint_copula"}
#'   (default as of \code{VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md} Phase 13 - dispatches to
#'   \code{preserve_correlation_structure_joint()}, drawing depth correlation and property
#'   correlation simultaneously) and \code{"gp_quantile_retrofit"} (the original algorithm,
#'   still fully supported as an explicit opt-out - dispatches to
#'   \code{preserve_correlation_structure()}). \code{NULL} (default) resolves to
#'   \code{"joint_copula"}, matching \code{get_monte_carlo_defaults()}'s own default - set
#'   \code{config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit"} explicitly to
#'   opt back into the original behavior. Under \code{"joint_copula"},
#'   \code{config$monte_carlo$vertical_correlation_gating} (default \code{FALSE}) separately
#'   controls whether \code{bound_sd}-based discontinuity gating (Phase 1c/1d) is applied - kept
#'   independent of the core method choice since its numeric defaults are not yet empirically
#'   calibrated (Phase 8).
#' @param gp_models Optional named list of fitted GP models (as `fit_local_gp_model_single()`
#'   returns), keyed by property - passed through to \code{preserve_correlation_structure_joint()}
#'   when \code{vertical_correlation_method = "joint_copula"}, so its depth kernel can reuse each
#'   property's already-fitted, already-cross-validated length-scale (see
#'   \code{extract_depth_length_scale()}) instead of falling back to a full-depth-range default.
#'   Ignored entirely under \code{"gp_quantile_retrofit"}.
#' @return Adjusted simulation data
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

  # Which vertical-correlation method to use. DEFAULT AS OF Phase 13: "joint_copula" (flipped from
  # "gp_quantile_retrofit" - see get_monte_carlo_defaults()'s own extended comment for the full
  # decision trail). A NULL/missing config, or a config that simply doesn't set this key, both
  # resolve to this same default - kept in sync with get_monte_carlo_defaults()'s own default so
  # "no config passed" means the same thing everywhere in this package, whether or not a caller
  # goes through get_monte_carlo_defaults() first. Explicitly set
  # config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit" to opt back into the
  # original algorithm.
  vertical_correlation_method <- config$monte_carlo$vertical_correlation_method %||% "joint_copula"

  # Discontinuity gating (build_depth_correlation_kernel()'s boundary_distinctness suppression,
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 1c/1d) is a SEPARATE opt-in from the core
  # joint_copula method itself (Phase 8) - bound_sd is attached unconditionally upstream
  # (attach_osd_boundary_distinctness() in simulate_ssurgo_mapunit_draws()), so without this flag
  # there would be no way to use joint_copula WITHOUT gating whenever OSD lookup succeeds. Its
  # numeric defaults (distinctness_range/min_gate_weight) are not yet empirically calibrated
  # against real KSSL/SSURGO lag correlations (see the decision-points section of
  # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md) - defaults to FALSE so the two decisions (core
  # method vs. gating strength) can be made independently.
  vertical_correlation_gating <- isTRUE(config$monte_carlo$vertical_correlation_gating)

  # Apply multivariate adjustment with enhanced error handling
  adjusted_matrices <- tryCatch({
    if (preserve_correlations && length(property_matrices) >= 2) {
      if (identical(vertical_correlation_method, "joint_copula")) {
        # bound_sd (OSD boundary distinctness - VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase
        # 1b) is broadcast identically across every simulation realization at a given depth
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

  # Merge with original data using Module 8 safe operations
  result_data <- merge_adjusted_data(cokey_data, adjusted_data, available_properties)

  return(result_data)
}

# ============================================================================
# 1b. VERTICAL-CORRELATION REDESIGN: JOINT DEPTH x PROPERTY COPULA (Phase 2)
# ============================================================================
#
# See VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md. `preserve_correlation_structure()` below
# retrofits vertical correlation onto already-independent-across-depth draws via a sequential
# gp_ratio nudge keyed off ONE "primary property"'s rank - the two functions here instead draw
# depth correlation (R_depth, from build_depth_correlation_kernel()) and property correlation
# (R_prop, the same flat matrix simulate_correlated_triangular() already uses) SIMULTANEOUSLY from
# a single Kronecker-separable joint distribution, so both are satisfied by construction rather
# than approximated by a retrofit. Phase 3's preserve_correlation_structure_joint() wires these two
# functions into a drop-in alternative with the same signature/contract as
# preserve_correlation_structure() itself.

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
#' @export
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
#' @export
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
      # intended "no variation (or no data) - keep original values" contract. Found via
      # VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 10's full-AOI benchmark against real
      # (messy) Salinas Valley SSURGO data - synthetic test fixtures never happened to include an
      # all-NA property column.
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

#' Preserve Correlation Structure via a Joint Depth x Property Copula (Phase 3)
#'
#' Drop-in alternative to `preserve_correlation_structure()` - same required parameters, in the
#' same order, so existing call sites work unchanged - that replaces its sequential
#' `gp_ratio`/single-"primary-property" retrofit with the joint Kronecker-copula sampler from
#' `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md` Phase 2 (`sample_joint_depth_property_copula()` +
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
#'   applied to every property with a usable model, averaged across them (reusing the
#'   already-fitted, already-cross-validated GP fits per `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md`
#'   Phase 1, rather than a new estimation step). When `NULL`/empty (e.g. this function is called
#'   standalone, without the fitted models available), falls back to a length-scale spanning the
#'   full depth range (`diff(range(depths))`) - a conservative "moderate smooth correlation across
#'   the whole profile" default that keeps this function usable without requiring GP models.
#' @param boundary_distinctness Optional per-depth `bound_sd` vector, passed through to
#'   `build_depth_correlation_kernel()` for discontinuity gating (Phase 1c/1d). `NULL` (default)
#'   skips gating.
#' @param kernel `"exponential"` (default) or `"matern"` - passed through to
#'   `build_depth_correlation_kernel()`.
#'
#' @return List of adjusted property matrices - same shape/contract as
#'   `preserve_correlation_structure()`'s return value. Degrades gracefully to the original,
#'   unadjusted `property_matrices` (with a warning) on insufficient dimensions or a sampling/
#'   mapping failure, matching `preserve_correlation_structure()`'s own graceful-failure contract.
#' @export
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
  # values - the same source adjust_multivariate_depthwise_GP()'s own verification step already
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
#' @param config Optional Monte Carlo config, passed through to `apply_gp_depth_trends()` -
#'   `config$monte_carlo$vertical_correlation_method` (default `"joint_copula"` as of
#'   `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md` Phase 13; set to `"gp_quantile_retrofit"` to opt
#'   back into the original algorithm) reaches the NRCS/regional GP path the same way it already
#'   reaches the local-GP path (Phase 6/11). `NULL` (default) resolves to `"joint_copula"`,
#'   matching `get_monte_carlo_defaults()`'s own default.
#' @return Adjusted simulation data
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

  # Extract the actual fitted NRCS GP model objects (not just their predictions) for the joint
  # method's depth-kernel length-scale reuse (extract_depth_length_scale(), Phase 1) - same
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

  # Apply GP depth trends using enhanced function
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
#' @param gp_control Passed through to `fit_local_gp_models()`/`fit_local_gp_model_single()`'s
#'   `gp_control` - see `fit_local_gp_model_single()`'s docs for why the default is much smaller
#'   than `GPfit::GP_fit()`'s own default.
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
                                       gp_control = c(20, 10, 2),
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
  # config$monte_carlo$vertical_correlation_method = "joint_copula"
  # (VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 4/6) is actually reachable from this
  # function's own already-fitted local_gp_models, not just by calling apply_gp_depth_trends()
  # directly.
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
#' Enhanced version with Module 8 validation and configuration management.
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
#' @export
fit_local_gp_models <- function(cokey_data, properties, config = NULL, gp_control = c(20, 10, 2),
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
#' Enhanced version with Module 8 error handling.
#'
#' @param cokey_data Simulation data
#' @param local_predictions Local GP predictions
#' @param unique_depths Depth vector
#' @param preserve_correlations Whether to preserve correlations
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so
#'   \code{INFO}-level progress messages print for the duration of this call (default
#'   \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @param config Optional config, passed straight through to `apply_gp_depth_trends()` - lets
#'   `config$monte_carlo$vertical_correlation_method` (default `"joint_copula"` as of
#'   `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md` Phase 13; set to `"gp_quantile_retrofit"` to opt
#'   back into the original algorithm) reach this call site. `NULL` (default) resolves to
#'   `"joint_copula"`, matching `get_monte_carlo_defaults()`'s own default.
#' @param gp_models Optional named list of fitted local GP models (as `apply_local_gp_adjustments()`
#'   already has in scope via `fit_local_gp_models()`), passed straight through to
#'   `apply_gp_depth_trends()` so the joint-copula depth kernel can reuse their fitted
#'   length-scales. Ignored under `"gp_quantile_retrofit"`.
#' @return Adjusted simulation data
#' @export
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
      # Vectorized lookup instead of a per-cell dplyr::filter()/pull() over the whole
      # data frame (was O(depths * sims * nrow(simulation_data)), with per-call dplyr/rlang
      # NSE overhead on top - profiling on a real AOI showed this single loop accounting for
      # ~60% of simulate_ssurgo_mapunit_draws()'s total runtime). match() against a combined
      # depth/simulation_number key finds the FIRST matching row for each (depth, sim_num)
      # cell, preserving the original's value[1] semantics exactly; unmatched or NA values
      # both correctly collapse to NA. Order of `target_key` matches matrix()'s column-major
      # fill (depth varies fastest, matching row index; sim_number varies slowest, matching
      # column index).
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

  metadata_cols <- intersect(
    c("cokey", "compname", "mukey", "hzdept_r", "hzdepb_r", "simulation_number", "unique_id"),
    names(original_data)
  )

  # Vectorized "first matching row" lookup instead of a per-(depth, sim_num) cell
  # dplyr::filter()/slice(1) over the whole original_data data frame - the same anti-pattern
  # already fixed in convert_to_property_matrices() (was O(depths * sims * nrow(original_data))
  # with per-call dplyr/rlang NSE overhead; profiling on a real AOI showed this the
  # next-largest remaining bottleneck at ~48% of simulate_ssurgo_mapunit_draws()'s total
  # runtime after that first fix). Grid order (depth outer/slower, sim_num inner/faster)
  # matches the original nested loop's row order exactly.
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

# Enhanced sequential processing with Module 8 progress tracking
process_cokeys_sequential <- function(cokey_groups, unique_cokeys, properties,
                                      gp_models, cokey_mapping, use_nrcs_gp, use_local_gp,
                                      preserve_correlations, config) {

  log_message("INFO", "Running integration sequentially", category = "MultivarAdjust")

  results <- list()

  for (i in seq_along(unique_cokeys)) {
    cokey <- unique_cokeys[i]

    # Progress tracking using Module 8
    track_progress(i, length(unique_cokeys), "Processing cokeys", update_frequency = 10)

    result <- process_single_cokey(cokey_groups[[as.character(cokey)]], cokey, properties, gp_models,
                                            cokey_mapping, use_nrcs_gp, use_local_gp,
                                            preserve_correlations, config)

    results[[i]] <- result
  }

  return(results)
}

# Enhanced single cokey processing with Module 8 error handling
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

      # config threaded through (VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md Phase 11) - previously
      # omitted entirely, so "joint_copula" was unreachable via the NRCS/regional GP path
      # regardless of what a caller's config said, even though the local-GP branch just below
      # already received it (Phase 6).
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
  # Vectorized: quantile() already accepts a vector of probs and computes every requested
  # quantile from a SINGLE sort of curr_values. The original per-replicate loop called
  # quantile(curr_values, probs = q, ...) once per j - curr_values never changes across
  # iterations, so this was n_sims separate full sorts of the same data instead of one.
  # Profiling on a real AOI showed this loop alone accounting for ~44% of the whole SSURGO
  # simulation pipeline's total runtime. quantile()'s failure modes (e.g. curr_values all-NA)
  # depend on curr_values as a whole, not on which individual prob was requested, so a single
  # tryCatch around the vectorized call is behavior-equivalent to the original's per-element
  # fallback - curr_values are matrix rows of length n_sims at every real call site, so
  # returning curr_values unadjusted on failure matches the original's per-j
  # "fall back to curr_values[j]" exactly.
  tryCatch({
    quantile_values <- stats::quantile(curr_values, probs = reference_quantiles, na.rm = TRUE, names = FALSE)
    quantile_values + (prev_values * gp_ratio - quantile_values)
  }, error = function(e) {
    curr_values
  })
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

#' @section Performance:
#' Previously did a per-row `which()` full-table scan of `result_data` for every row of
#' `adjusted_data` (an O(n_adjusted x n_result) join) - `Rprof()` profiling (10,000-row synthetic
#' benchmark) confirmed this as a real per-cokey hot path cost
#' (PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4). Replaced with a single vectorized key match. A row
#' only updates `result_data` when its (`hzdept_r`, `simulation_number`) key matches **exactly
#' one** `result_data` row (the original's `length(match_idx) == 1` contract, silently preserved
#' - zero or multiple matches are skipped, not an error) and its own value for that property is
#' non-`NA`. When multiple `adjusted_data` rows share the same key, R's vectorized `[<-`
#' assignment applies them in order and the last one wins - verified to match the original
#' sequential loop's last-write-wins behavior exactly (confirmed empirically:
#' `x[c(2,2,3)] <- c(10,20,30)` yields `x[2] == 20`, not `10`).
merge_adjusted_data <- function(cokey_data, adjusted_data, available_properties) {
  result_data <- cokey_data

  key_result <- paste(result_data$hzdept_r, result_data$simulation_number, sep = "\a")
  key_adj <- paste(adjusted_data$hzdept_r, adjusted_data$simulation_number, sep = "\a")

  # A key must be unique WITHIN result_data (exactly one candidate row) to match the original
  # "length(match_idx) == 1" contract - match() alone only finds the first occurrence and can't
  # tell duplicates from a genuine single match.
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

#' @section Performance:
#' Previously scaled texture properties one row at a time via `data[i, texture_props]`
#' data.frame row-slicing - `Rprof()`-free benchmarking alone made this obvious (31.61s for
#' 50,000 synthetic rows despite a trivial per-row body, PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4),
#' matching the same anti-pattern already fixed in `related_property_estimation()`'s texture
#' branch (195x there). Replaced with a `rowSums()`-based vectorization operating on the whole
#' texture-column matrix at once. `texture_sum > 0 & !is.na(texture_sum)` preserves the
#' original's exact `&&`-based NA handling (R's `&`/`&&` both resolve `NA & FALSE` to `FALSE`,
#' so an `NA` sum still correctly skips scaling for that row either way).
apply_cross_property_constraints <- function(data, properties) {
  # Enhanced texture sum constraint using Module 8 validation
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
