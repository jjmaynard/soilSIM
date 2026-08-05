#' @title Van Genuchten / ROSETTA-Based Available Water Storage (AWS) Modeling
#'
#' @description Closed-form van Genuchten water-retention evaluation, Monte
#'   Carlo simulation of available water holding capacity (AWHC) from
#'   `soilDB::ROSETTA()`-derived pedotransfer parameters, and depth-sliced
#'   available-water-storage summaries for a set of soil components. Ported
#'   from `code_ref/brdf/property_simulation.R` (identical to the superseded
#'   `code/sim-functions.R` copy). Self-contained: unlike
#'   `R/depth-simulation.R`'s functions, none of these three depend on any
#'   other not-yet-ported helper from that file (e.g. `sim_component_comp()`).
#' @name aws_simulation
NULL

#' Evaluate the van Genuchten Water Retention Curve
#'
#' Closed-form van Genuchten (1980) volumetric water content at a given
#' matric potential.
#'
#' @param h Matric potential (cmH2O, typically negative).
#' @param alpha,n Van Genuchten shape parameters.
#' @param theta_r,theta_s Residual and saturated volumetric water content.
#'
#' @return Numeric vector of volumetric water content, `theta(h)`.
#'
#' @examples
#' van_genuchten(h = -100, alpha = 0.02, n = 1.3, theta_r = 0.05, theta_s = 0.45)
#' @export
van_genuchten <- function(h, alpha, n, theta_r, theta_s) {
  m <- 1 - (1 / n)
  as.numeric(theta_r + (theta_s - theta_r) / ((1 + (abs(alpha * h))^n)^m))
}


#' Simulate Available Water Holding Capacity from van Genuchten Parameters
#'
#' Monte Carlo-simulates field capacity (FC) and permanent wilting point
#' (PWP) water content - and their difference, available water holding
#' capacity (AWHC) - for each row (soil layer/component) of van Genuchten
#' parameters, by drawing `n_simulations` random samples per parameter from
#' `stats::rnorm(mean, sd)` and evaluating `van_genuchten()` at the FC/PWP
#' matric potentials.
#'
#' @section Ported-as-is behavioral quirks:
#' Two aspects of this function are preserved exactly as in the legacy source
#' rather than "fixed," since both are plausibly intentional design choices
#' (not unconditional-crash bugs) and changing either would alter the actual
#' simulated values every caller receives:
#' - `set.seed(123)` is called **inside** the per-row loop, so every row
#'   re-seeds and draws from an identically-seeded random stream rather than
#'   an evolving one across rows.
#' - Sampled `alpha`/`n` values are back-transformed via `10^(...)`, i.e. the
#'   input `data$alpha`/`data$npar` (and their SDs) are treated as already
#'   being in log10 space. This is a standard technique for keeping van
#'   Genuchten shape parameters positive after adding Gaussian noise, but
#'   whether `soilDB::ROSETTA()`'s actual reported values are meant to be
#'   interpreted this way is not verified by this port.
#'
#' @param data A data frame with one row per soil layer/component, and
#'   columns `alpha`, `sd_alpha`, `npar`, `sd_npar`, `theta_r`, `sd_theta_r`,
#'   `theta_s`, `sd_theta_s`, `layerID` - exactly the shape produced by
#'   `soilDB::ROSETTA(..., include.sd = TRUE)` plus a caller-added `layerID`
#'   (see `calculate_aws_df()`).
#' @param n_simulations Number of Monte Carlo draws per row (default 100).
#'
#' @return A named list (one element per row, keyed `paste(layerID, i,
#'   sep = "_")`), each a data frame of `n_simulations` draws with columns
#'   `alpha`, `n`, `theta_r`, `theta_s`, `sim_num`, `theta_fc`, `theta_pwp`,
#'   `AWHC`. Rows with any missing van Genuchten parameter are skipped
#'   (absent from the result).
#'
#' @export
simulate_vg_aws <- function(data, n_simulations = 100) {

  # Initialize result list to store the simulations for each component
  simulation_results <- list()

  # Iterate over each row in the data
  for (i in seq_len(nrow(data))) {

    # Extract van Genuchten parameters and their uncertainties for the current component
    alpha_mean    <- data$alpha[i]
    alpha_sd      <- data$sd_alpha[i]
    n_mean        <- data$npar[i]
    n_sd          <- data$sd_npar[i]
    theta_r_mean  <- data$theta_r[i]
    theta_r_sd    <- data$sd_theta_r[i]
    theta_s_mean  <- data$theta_s[i]
    theta_s_sd    <- data$sd_theta_s[i]

    # If any of the parameters are NA, skip this row
    if (is.na(alpha_mean) | is.na(n_mean) | is.na(theta_r_mean) | is.na(theta_s_mean)) {
      next
    }

    # Generate random samples for the van Genuchten parameters based on the provided means and standard deviations
    set.seed(123)  # Set seed for reproducibility
    alpha_samples   <- stats::rnorm(n_simulations, mean = alpha_mean, sd = alpha_sd)
    n_samples       <- stats::rnorm(n_simulations, mean = n_mean, sd = n_sd)
    theta_r_samples <- stats::rnorm(n_simulations, mean = theta_r_mean, sd = theta_r_sd)
    theta_s_samples <- stats::rnorm(n_simulations, mean = theta_s_mean, sd = theta_s_sd)

    # Convert matric potentials from kPa to cmH2O (1 kPa = 10.19716 cmH2O)
    h_fc  <- -33 * 10.19716    # Field capacity (FC) matric potential in cmH2O
    h_pwp <- -1500 * 10.19716  # Permanent wilting point (PWP) matric potential in cmH2O

    # Create a dataframe to hold the sampled parameters
    results <- data.frame(
      alpha   = 10^(alpha_samples),   # Convert alpha back from log scale
      n       = 10^(n_samples),       # Convert npar back from log scale
      theta_r = theta_r_samples,
      theta_s = theta_s_samples,
      sim_num = seq_len(n_simulations)
    )

    # Calculate water retention at FC and PWP for each set of sampled parameters using the van Genuchten function
    results <- results |>
      dplyr::rowwise() |>
      dplyr::mutate(
        theta_fc  = van_genuchten(h_fc, alpha, n, theta_r, theta_s),   # Water content at Field Capacity
        theta_pwp = van_genuchten(h_pwp, alpha, n, theta_r, theta_s)    # Water content at Permanent Wilting Point
      ) |>
      dplyr::ungroup()

    # Calculate Available Water Holding Capacity (AWHC) as the difference between FC and PWP
    results$AWHC <- results$theta_fc - results$theta_pwp

    # Store the results in the list with a key combining component name and layerID
    simulation_results[[paste(data$layerID[i], i, sep = "_")]] <- results
  }

  # Return the list of simulation results for each component
  return(simulation_results)
}


#' Minimal na.rm-mean `slab.fun` for `aqp::slab()`
#'
#' Replacement for the no-longer-available `aqp::mean_na()` (removed/renamed
#' in current `aqp` versions) - preserves its old single-value-per-slab
#' contract, which `aqp::slab()`'s current default `slab.fun`
#' (`slab_function(method = "numeric")`) does not (it returns quantile
#' columns instead of a single `value` column).
#'
#' @param values Numeric vector of observations within one depth slab.
#' @param ... Ignored (`aqp::slab()` passes additional arguments positionally).
#' @return A single numeric value, `mean(values, na.rm = TRUE)`.
#' @noRd
.aws_slab_mean <- function(values, ...) mean(values, na.rm = TRUE)

#' Calculate Available Water Storage by Depth Interval
#'
#' Runs `soilDB::ROSETTA()` on soil texture/bulk-density/water-retention
#' inputs to derive van Genuchten pedotransfer parameters, Monte
#' Carlo-simulates available water holding capacity per horizon via
#' `simulate_vg_aws()`, then depth-slices and summarizes mean AWHC per
#' component over standard depth intervals via `aqp::slab()`.
#'
#' @section Known limitation:
#' This function requires **live network access** - `soilDB::ROSETTA()`
#' POSTs to `https://www.handbook60.org/api/v1/rosetta/<version>` (via
#' `httr::POST()`) rather than computing pedotransfer parameters locally.
#' Calling this without internet connectivity, or if `handbook60.org` is
#' unreachable, will fail.
#'
#' @param sim_data_df A data frame with one row per horizon, with columns
#'   `sand_total`, `silt_total`, `clay_total`, `bulk_density_third_bar`,
#'   `water_retention_third_bar`, `water_retention_15_bar` (ROSETTA's
#'   expected input variable names - note these differ from the rest of
#'   `soilSIM`'s `sandtotal_r`/`claytotal_r`/etc. SSURGO-derived naming
#'   convention, since they're ROSETTA's own API contract), plus `compname`,
#'   `hzdept_r`, `hzdepb_r`, `cokey`.
#'
#' @return A long-format data frame with one row per `cokey` per depth slab
#'   actually spanned by that component's horizons (slab boundaries per
#'   `slab.structure = c(0, 5, 15, 30, 60, 100)`), with columns `cokey`,
#'   `top`, `bottom`, `AWHC` (mean available water holding capacity over that
#'   slab). The `tidyr::pivot_wider()` step widens on `variable` (always
#'   `"AWHC"` here, since only one property is slabbed), so it does not
#'   collapse depth slabs into columns - it exists to match `aqp::slab()`'s
#'   long output shape to a plain `value`-column contract.
#'
#' @export
calculate_aws_df <- function(sim_data_df) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("calculate_aws_df() requires the 'httr' package (used internally by soilDB::ROSETTA()) to be installed.")
  }

  # Run Rosetta using soil properties
  variables <- c(
    "sand_total", "silt_total", "clay_total",
    "bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar"
  )
  rosetta_data <- soilDB::ROSETTA(sim_data_df, vars = variables, v = "3", include.sd = TRUE)

  # Create layerID by combining compname and hzdept_r
  sim_data_df$layerID <- paste(sim_data_df$compname, sim_data_df$hzdept_r, sep = "_")
  rosetta_data$layerID <- sim_data_df$layerID

  # Simulate available water storage using van Genuchten parameters
  aws <- simulate_vg_aws(rosetta_data)

  # Use Map to add new columns to each data frame in the list
  aws <- Map(function(df, top_value, bottom_value, compname_value, cokey_value) {
    df$top <- top_value
    df$bottom <- bottom_value
    df$compname <- compname_value
    df$cokey <- cokey_value
    return(df)
  }, aws, sim_data_df$hzdept_r, sim_data_df$hzdepb_r, sim_data_df$compname, sim_data_df$cokey)

  # Combine the list of data frames into one
  sim_aws_df <- dplyr::bind_rows(aws)

  # Assign depth structure to the data frame using aqp. Note: the Map() step
  # above attaches the depth boundaries as `top`/`bottom` (not `hzdept_r`/
  # `hzdepb_r`) - the formula below must reference the columns that actually
  # exist on sim_aws_df, not sim_data_df's original column names.
  aqp::depths(sim_aws_df) <- cokey ~ top + bottom

  # Compute slab summaries over specified depth intervals (e.g., 0, 5, 15, 30, 60, 100 cm).
  # `aqp::mean_na` no longer exists in current aqp versions (this legacy code
  # predates its removal/rename) - `.aws_slab_mean()` below is a minimal
  # local replacement with the same contract (a single na.rm-mean per slab,
  # yielding a `value` column), since aqp's current default `slab.fun`
  # (`slab_function(method = "numeric")`) instead returns quantile columns
  # with no `value` column, which would break the pivot_wider() step below.
  prof.slab <- aqp::slab(sim_aws_df,
                         fm = cokey ~ AWHC,
                         slab.structure = c(0, 5, 15, 30, 60, 100),
                         slab.fun = .aws_slab_mean)

  # Reshape the slab output to wide format
  prof.slab.1 <- prof.slab |>
    dplyr::select(-contributing_fraction) |>
    tidyr::pivot_wider(names_from = variable, values_from = value)

  # Clean up missing values if necessary (this line may be adjusted based on desired behavior)
  prof.slab.1[is.na(prof.slab.1)] <- NA

  return(prof.slab.1)
}
