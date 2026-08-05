#' @title Soil Profile Depth Simulation
#'
#' @description SSURGO-based soil-profile horizon-depth simulation: top-down/bottom-up
#'   depth simulation, thickness-variability estimation via repeated simulation,
#'   OSD-distinctness-derived boundary perturbation, and
#'   `aqp::SoilProfileCollection`-level orchestration (single profile, whole
#'   collection, parallelized collection, or by-mukey). Ported from
#'   `code_ref/brdf/depth_simulation.R`, the canonical, generalized version of
#'   this logic (superseding earlier copies kept for historical reference in
#'   `code/scratch/`, not part of this port).
#'
#'   This is the first place `aqp`/`SoilProfileCollection` objects are used in
#'   `soilSIM`; the random triangular-distribution sampler these functions rely
#'   on, `tri_dist()`, lives in `R/distributions.R` (it was originally defined
#'   in the *companion* legacy file, `code_ref/brdf/property_simulation.R`, not
#'   in `depth_simulation.R` itself). That companion file's correlated-property
#'   and van Genuchten/AWS modeling functions remain out of scope for this port.
#' @name depth_simulation
NULL

#' Query SSURGO Database for Soil Data by Mukey
#'
#' Queries the NRCS SSURGO database (`mapunit`/`component`/`chorizon`/`chfrags`)
#' for horizon morphology, texture, bulk density, water-retention, and
#' rock-fragment-volume data for the given map unit key(s).
#'
#' This deliberately does not reuse `execute_ssurgo_query_working()` (which
#' also joins `mapunit`/`component`/`chorizon`/`chfrags`): that function only
#' ever selects `chf.fragsize_r`, never `chf.fragvol_l/r/h`, so its rock
#' fragment aggregation path (`aggregate_rock_fragment_volume_working()`)
#' currently produces no data. This function's entire purpose is summing
#' `fragvol_l/r/h` per horizon, so delegating would silently drop rock
#' fragment volume - a real regression, not a style difference.
#'
#' @param mukeys A vector of string or numeric map unit keys (mukeys) representing the map units to retrieve data for.
#'
#' @return A data frame with soil horizon data for the given mukeys, including aggregated rock fragment data.
#'
#' @export
get_aws_data_by_mukey <- function(mukeys) {
  # Step 1: Format mukey(s) for SQL
  formatted_mukey <- paste0("(", paste0("'", mukeys, "'", collapse = ","), ")")

  # Step 2: Construct the SQL query
  query <- sprintf("
    SELECT
      -- contextual data
      co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
      -- horizon morphology
      ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
      -- soil texture parameters
      sandtotal_l, sandtotal_r, sandtotal_h,
      silttotal_l, silttotal_r, silttotal_h,
      claytotal_l, claytotal_r, claytotal_h,
      -- bulk density and water retention parameters
      dbovendry_l, dbovendry_r, dbovendry_h,
      wthirdbar_l, wthirdbar_r, wthirdbar_h,
      wfifteenbar_l, wfifteenbar_r, wfifteenbar_h,
      -- rock fragment volume from chfrags table
      chf.fragvol_l AS rfv_l, chf.fragvol_r AS rfv_r, chf.fragvol_h AS rfv_h
    FROM mapunit AS mu
      INNER JOIN component AS co ON mu.mukey = co.mukey
      INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
      LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey
    WHERE mu.mukey IN %s
    ORDER BY co.cokey, ch.hzdept_r ASC;", formatted_mukey)

  # Step 3: Submit the query and fetch results
  result <- soilDB::SDA_query(query)

  # Step 4: Group by `chkey` and sum rock fragment values (rfv_l, rfv_r, rfv_h)
  rfv_sum <- result |>
    dplyr::group_by(chkey) |>
    dplyr::summarise(
      rfv_l = sum(rfv_l, na.rm = TRUE),
      rfv_r = sum(rfv_r, na.rm = TRUE),
      rfv_h = sum(rfv_h, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Step 5: Merge the summed values back to the main dataset and remove duplicates
  result <- result |>
    dplyr::select(-dplyr::all_of(c("rfv_l", "rfv_r", "rfv_h"))) |>
    dplyr::left_join(rfv_sum, by = "chkey") |>
    dplyr::distinct()

  # Step 6: Return the final result
  return(result)
}


#--------------------------------------------------------------------------------------------------------
# Functions to run soil profile depth simulations

#' Query OSD Data and Convert Horizon Distinctness to Offset
#'
#' Retrieves Official Series Description (OSD) data for the soil series
#' (component names) present in `horizon_data`, extracts horizon distinctness,
#' and converts the distinctness codes into boundary-offset values via
#' `aqp::hzDistinctnessCodeToOffset()`.
#'
#' @param horizon_data Data frame containing SSURGO horizon data, with at
#'   least `compname` and `hzname` columns.
#'
#' @return A data frame containing the following columns:
#'         - `id`: Unique identifier for the soil profile (the component name).
#'         - `hzname`: The name of the soil horizon.
#'         - `distinctness`: The distinctness code for the horizon (e.g., A, C, G, D).
#'         - `genhz`: The generalized horizon designation.
#'         - `bound_sd`: The calculated offset value from the distinctness code.
#'
#' @examples
#' \dontrun{
#'   horizon_data <- data.frame(
#'     compname = c("amador", "pentz", "pardee", "auburn", "loafercreek", "millvilla"),
#'     hzname = c("A", "Bt", "R", "A", "Bw", "C")
#'   )
#'   result <- query_osd_distinctness(horizon_data)
#'   head(result)
#' }
#' @export
query_osd_distinctness <- function(horizon_data) {
  # Subset SSURGO horizon data to keep only compname and hzname
  ssurgo_horizon_data <- horizon_data |>
    dplyr::select(compname, hzname)

  # Ensure 'genhz' column exists in ssurgo_horizon_data
  if (!"genhz" %in% colnames(ssurgo_horizon_data)) {
    ssurgo_horizon_data$genhz <- aqp::generalizeHz(
      ssurgo_horizon_data$hzname,
      new = c('O', 'A', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^A', '^B', '^C', '^Cr', '^R'),
      ordered = TRUE
    )
  }

  # Rename 'compname' to 'id'
  colnames(ssurgo_horizon_data)[colnames(ssurgo_horizon_data) == "compname"] <- "id"

  # Fetch OSD data for the list of soil component names, reusing cached per-series results
  # where available (see fetch_osd_horizons_cached()'s Known quirk note re: id casing) -
  # avoids redundant live soilDB::fetchOSD() calls when the same series appears in multiple
  # profiles/components within or across sessions (e.g. simulate_profile_depths_by_mukey()
  # calls this once per component, and components frequently share a series name).
  osd_horizon_data <- fetch_osd_horizons_cached(unique(horizon_data$compname))

  # Generate generalized horizon codes for the OSD data
  osd_horizon_data$genhz <- aqp::generalizeHz(
    osd_horizon_data$hzname,
    new = c('O', 'A', 'B', 'C', 'Cr', 'R'),
    pattern = c('O', '^\\d*A', '^\\d*B', '^\\d*C', '^\\d*Cr', '^\\d*R'),
    ordered = TRUE
  )

  # Identify missing genhz values in ssurgo_horizon_data that are not in osd_horizon_data
  missing_genhz <- setdiff(ssurgo_horizon_data$genhz, osd_horizon_data$genhz)

  # If there are missing genhz values, append them with NA for distinctness
  if (length(missing_genhz) > 0) {
    missing_rows <- ssurgo_horizon_data[ssurgo_horizon_data$genhz %in% missing_genhz, c("id", "hzname", "genhz")]
    missing_rows$distinctness <- NA
    osd_horizon_data <- base::rbind(osd_horizon_data, missing_rows)
    message("Added missing genhz values: ", paste(missing_genhz, collapse = ", "))
  } else {
    message("No missing genhz values found.")
  }

  # Infill missing distinctness values using a custom function
  osd_horizon_data <- infill_missing_distinctness(osd_horizon_data)

  # Convert the distinctness codes to offset values
  osd_horizon_data$bound_sd <- aqp::hzDistinctnessCodeToOffset(osd_horizon_data$distinctness)

  # Rename the columns for clarity
  names(osd_horizon_data) <- c("id", "hzname", "distinctness", "genhz", "bound_sd")

  return(osd_horizon_data)
}

#' Fetch OSD horizon distinctness data for a set of soil series names, with per-series caching
#'
#' `soilDB::fetchOSD()` is a live network call. Multiple profiles/components in the same
#' mukey/AOI frequently share the same series name (e.g. `simulate_profile_depths_by_mukey()`
#' calls `query_osd_distinctness()` once per component), so fetching unconditionally on every
#' call re-requests identical data redundantly. This wrapper caches each series' raw
#' `id`/`hzname`/`distinctness` rows on disk - the same `cache_get()`/`cache_set()` mechanism
#' `R/raster-cache.R` already uses for AOI-keyed raster fetches, just with a series-name-keyed
#' cache key instead - fetching only the series not already cached.
#'
#' @section Known quirk preserved (not a bug introduced here): `soilDB::fetchOSD()` returns its
#'   `id` column upper-cased regardless of the requested case (e.g. requesting `"amador"`
#'   returns rows with `id == "AMADOR"`) - confirmed via a live call before writing this
#'   function. Splitting a combined multi-series fetch into per-series cache entries is done via
#'   case-insensitive matching on `id` to account for this. This casing difference was already
#'   present in `soilDB::fetchOSD()`'s output before this caching layer existed, and downstream
#'   code (`simulate_and_perturb_soil_profiles()`) was already unaffected by it, since it joins
#'   `query_osd_distinctness()`'s output by `genhz` only, never by `id`.
#'
#' @param compnames Character vector of soil series names, as would be passed directly to
#'   `soilDB::fetchOSD()`.
#' @return A data frame with columns `id`, `hzname`, `distinctness` - one row per horizon, across
#'   all requested series (any already-cached plus any newly fetched).
#' @keywords internal
fetch_osd_horizons_cached <- function(compnames) {
  compnames <- unique(compnames)
  keys <- paste0("osd_series_", vapply(toupper(compnames), digest::digest, character(1)))
  cached <- stats::setNames(lapply(keys, cache_get), compnames)

  missing <- compnames[vapply(cached, is.null, logical(1))]
  if (length(missing) > 0) {
    s <- soilDB::fetchOSD(missing)

    if (!all(c("distinctness", "top", "bottom", "hzname", "id") %in% names(s@horizons))) {
      stop("Necessary fields (distinctness, id, top, bottom, hzname) are not available in the OSD data.")
    }

    fetched <- s@horizons[, c("id", "hzname", "distinctness")]
    for (cn in missing) {
      rows <- fetched[toupper(fetched$id) == toupper(cn), , drop = FALSE]
      key <- paste0("osd_series_", digest::digest(toupper(cn)))
      cache_set(key, "osd_series", rows)
      cached[[cn]] <- rows
    }
  }

  do.call(rbind, cached)
}


#' Infill Missing Distinctness Values for Horizons
#'
#' Fills in missing `distinctness` values for common horizon names (O, A, B, C, R, Cr)
#' based on typical boundary characteristics observed in the field.
#'
#' ## Default Boundary Distinctness Assignments:
#' - **O Horizons (Organic Layers):** diffuse
#' - **A Horizons (Topsoil):** clear
#' - **B Horizons (Subsoil):** gradual
#' - **C Horizons (Parent Material):** gradual
#' - **Cr Horizons (Weathered Bedrock):** gradual
#' - **R Horizons (Bedrock):** abrupt
#'
#' @param horizon_data A data frame containing horizon data. It must include `hzname` and `distinctness` columns.
#'
#' @return A data frame with missing `distinctness` values infilled based on the horizon name or generalized horizon group.
#'
#' @examples
#' \dontrun{
#'   df <- data.frame(
#'     hzname = c("A", "B", "C", "R", "O", "Cr"),
#'     distinctness = c(NA, "gradual", NA, NA, "diffuse", NA),
#'     stringsAsFactors = FALSE
#'   )
#'   df <- infill_missing_distinctness(df)
#'   print(df)
#' }
#' @export
infill_missing_distinctness <- function(horizon_data) {
  # List of default boundary distinctness values for common horizons
  default_distinctness <- list(
    "O"  = "diffuse",  # Organic horizon
    "A"  = "clear",    # Topsoil
    "E"  = "clear",    # Subsoil
    "B"  = "gradual",  # Subsoil
    "C"  = "gradual",  # Parent material
    "Cr" = "gradual",  # Weathered bedrock
    "R"  = "abrupt"    # Bedrock
  )

  # Assign a generalized horizon group (genhz) if it's not already present
  if (!"genhz" %in% colnames(horizon_data)) {
    horizon_data$genhz <- aqp::generalizeHz(
      horizon_data$hzname,
      new = c('O', 'A', 'E', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^\\d*A', '^\\d*E', '^\\d*B', '^\\d*C', '^\\d*Cr', '^\\d*R'),
      ordered = TRUE
    )
  }

  # Function to assign default distinctness based on horizon name or generalized horizon group
  fill_default_distinctness <- function(hzname, genhz, distinctness) {
    if (!is.na(distinctness)) {
      return(distinctness)  # If distinctness exists, return it
    }
    # Use default distinctness based on genhz first, then hzname if available
    if (!is.na(genhz) && genhz %in% names(default_distinctness)) {
      return(default_distinctness[[genhz]])
    } else if (!is.na(hzname) && hzname %in% names(default_distinctness)) {
      return(default_distinctness[[hzname]])
    } else {
      return(NA)  # Return NA if no match is found
    }
  }

  # Apply the default distinctness infilling using mapply
  horizon_data$distinctness <- mapply(
    fill_default_distinctness,
    horizon_data$hzname,
    horizon_data$genhz,
    horizon_data$distinctness
  )

  return(horizon_data)
}


#' Infill Missing Depth Variability
#'
#' Infills missing depth variability values in a horizon data frame by
#' replacing missing top and bottom depth bounds with the representative
#' depth value plus or minus 2, ensuring that the resulting values are not
#' negative.
#'
#' @param horizon_data A data frame containing horizon depth data. It must include the following columns:
#'   - `hzdept_r`: Representative top depth.
#'   - `hzdept_l`: Top low value (may be missing).
#'   - `hzdept_h`: Top high value (may be missing).
#'   - `hzdepb_r`: Representative bottom depth.
#'   - `hzdepb_l`: Bottom low value (may be missing).
#'   - `hzdepb_h`: Bottom high value (may be missing).
#'
#' @return A data frame with missing depth variability values filled in. Missing top low values are replaced with
#'   `hzdept_r - 2` and missing top high values with `hzdept_r + 2`. Similarly, missing bottom low values are replaced
#'   with `hzdepb_r - 2` and missing bottom high values with `hzdepb_r + 2`. All values are ensured to be non-negative.
#'
#' @examples
#' data <- data.frame(
#'   hzdept_r = c(10, 20),
#'   hzdept_l = c(NA, 18),
#'   hzdept_h = c(NA, NA),
#'   hzdepb_r = c(30, 40),
#'   hzdepb_l = c(NA, NA),
#'   hzdepb_h = c(NA, 42)
#' )
#' infill_missing_depth_variability(data)
#' @export
infill_missing_depth_variability <- function(horizon_data) {
  # Infill missing top low values (hzdept_l) with hzdept_r - 2, ensuring non-negative results
  horizon_data$hzdept_l <- ifelse(
    is.na(horizon_data$hzdept_l),
    pmax(horizon_data$hzdept_r - 2, 0),
    horizon_data$hzdept_l
  )

  # Infill missing top high values (hzdept_h) with hzdept_r + 2
  horizon_data$hzdept_h <- ifelse(
    is.na(horizon_data$hzdept_h),
    pmax(horizon_data$hzdept_r + 2, 0),
    horizon_data$hzdept_h
  )

  # Infill missing bottom low values (hzdepb_l) with hzdepb_r - 2
  horizon_data$hzdepb_l <- ifelse(
    is.na(horizon_data$hzdepb_l),
    pmax(horizon_data$hzdepb_r - 2, 0),
    horizon_data$hzdepb_l
  )

  # Infill missing bottom high values (hzdepb_h) with hzdepb_r + 2
  horizon_data$hzdepb_h <- ifelse(
    is.na(horizon_data$hzdepb_h),
    pmax(horizon_data$hzdepb_r + 2, 0),
    horizon_data$hzdepb_h
  )

  return(horizon_data)
}


#' Top-down Simulation of Soil Profile
#'
#' Simulates a soil profile by calculating the top and bottom depths for each
#' horizon in a top-down manner. The first horizon always starts at 0 cm, and
#' the top of each subsequent horizon is set equal to the bottom depth of the
#' previous horizon. The bottom depth of each horizon is simulated using a
#' triangular distribution via `tri_dist()`, defined by a lower bound
#' (`hzdepb_l`), a representative value (`hzdepb_r`), and an upper bound
#' (`hzdepb_h`).
#'
#' @param horizon_data A data frame containing horizon data. It should include the following columns:
#'   - `hzname`: Name of the horizon.
#'   - `hzdepb_l`: Lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: Representative bottom depth of the horizon.
#'   - `hzdepb_h`: Upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @examples
#' horizon_data <- data.frame(
#'   hzname = c("A", "B", "C"),
#'   hzdepb_l = c(20, 40, 60),
#'   hzdepb_r = c(25, 45, 65),
#'   hzdepb_h = c(30, 50, 70)
#' )
#' profile <- simulate_soil_profile_top_down(horizon_data)
#' print(profile)
#' @export
simulate_soil_profile_top_down <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Set the top of the first horizon to 0
  result$top[1] <- 0

  for (i in 1:n) {
    if (i > 1) {
      # Set the top of the current horizon to be the bottom of the previous horizon
      result$top[i] <- result$bottom[i - 1]
    }

    # Determine bounds for bottom depth simulation
    bottom_low  <- max(horizon_data$hzdepb_l[i], result$top[i])
    bottom_high <- max(bottom_low, horizon_data$hzdepb_h[i])
    bottom_mode <- min(max(horizon_data$hzdepb_r[i], bottom_low), bottom_high)

    # Generate a single random sample using the tri_dist function
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # Ensure the bottom depth is not less than the top depth
    if (bottom_depth < result$top[i]) {
      bottom_depth <- result$top[i]
    }
    result$bottom[i] <- bottom_depth
  }

  return(result)
}


#' Bottom-up Simulation of Soil Profile
#'
#' Simulates a soil profile from the bottom upward by calculating the top and
#' bottom depths for each horizon. The simulation starts from the bottom
#' horizon and proceeds upward, setting the bottom depth of each horizon to
#' the top depth of the horizon immediately below. For each horizon, the
#' bottom and top depths are simulated using a triangular distribution via
#' `tri_dist()`. The function ensures that there are no gaps between horizons
#' and that the top of the uppermost horizon is set to 0.
#'
#' @param horizon_data A data frame containing horizon data. The data frame must include the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `hzdept_l`: The lower bound for the top depth of the horizon.
#'   - `hzdept_r`: The representative top depth of the horizon.
#'   - `hzdept_h`: The upper bound for the top depth of the horizon.
#'   - `hzdepb_l`: The lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: The representative bottom depth of the horizon.
#'   - `hzdepb_h`: The upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @examples
#' horizon_data <- data.frame(
#'   hzname = c("A", "B", "C"),
#'   hzdept_l = c(0, 20, 40),
#'   hzdept_r = c(0, 25, 45),
#'   hzdept_h = c(0, 30, 50),
#'   hzdepb_l = c(20, 40, 60),
#'   hzdepb_r = c(25, 45, 65),
#'   hzdepb_h = c(30, 50, 70)
#' )
#' profile <- simulate_soil_profile_bottom_up(horizon_data)
#' print(profile)
#' @export
simulate_soil_profile_bottom_up <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Process horizons from bottom to top
  for (i in n:1) {
    # Simulate bottom depth using tri_dist
    bottom_low  <- horizon_data$hzdepb_l[i]
    bottom_high <- horizon_data$hzdepb_h[i]
    bottom_mode <- horizon_data$hzdepb_r[i]

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    if (i == n) {
      # For the bottom horizon, simulate top depth
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    } else {
      # For other horizons, set bottom depth to the top depth of the horizon below
      bottom_depth <- result$top[i + 1]

      # Simulate top depth using tri_dist
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    }

    # Assign the simulated depths to the result
    result$top[i]    <- top_depth
    result$bottom[i] <- bottom_depth
  }

  # Ensure the top depth of the uppermost horizon is 0
  result$top[1] <- 0

  # Adjust the bottom depth of the first horizon if necessary
  if (result$bottom[1] < result$top[1]) {
    result$bottom[1] <- result$top[1]
  }

  # Ensure no gaps between horizons
  for (i in 1:(n - 1)) {
    if (result$bottom[i] != result$top[i + 1]) {
      result$bottom[i] <- result$top[i + 1]
    }
  }

  return(result)
}


#' Simulate Soil Profile Thickness
#'
#' Simulates the thickness of soil horizons for a given soil profile using
#' both top-down and bottom-up simulation methods. It performs multiple
#' simulations and calculates the standard deviation of horizon thickness
#' (i.e., the difference between representative bottom and top depths)
#' across the simulations.
#'
#' @param horizon_data A dataframe containing horizon data. It must include the following columns:
#'        - `hzname`: The name of each horizon (e.g., A, B, C).
#'        - `hzdept_r`: The representative top depth of each horizon.
#'        - `hzdepb_r`: The representative bottom depth of each horizon.
#'        - `hzdept_l`: The low estimate for the top depth of each horizon.
#'        - `hzdept_h`: The high estimate for the top depth of each horizon.
#'        - `hzdepb_l`: The low estimate for the bottom depth of each horizon.
#'        - `hzdepb_h`: The high estimate for the bottom depth of each horizon.
#' @param n_simulations Integer, number of simulations to perform (default is 500).
#'
#' @return A dataframe containing the following columns:
#'         - `hzname`: The horizon name.
#'         - `top`: The representative top depth of the horizon (from `hzdept_r`).
#'         - `bottom`: The representative bottom depth of the horizon (from `hzdepb_r`).
#'         - `thickness_sd`: The standard deviation of the horizon thickness (bottom - top) across the simulations.
#'
#' @examples
#' horizon_data <- data.frame(
#'   hzname = c("A", "B", "C"),
#'   hzdept_r = c(0, 20, 35),
#'   hzdepb_r = c(20, 35, 50),
#'   hzdept_l = c(0, 15, 30),
#'   hzdept_h = c(0, 25, 40),
#'   hzdepb_l = c(15, 30, 45),
#'   hzdepb_h = c(25, 40, 60)
#' )
#' set.seed(123)
#' summarized_results <- simulate_soil_profile_thickness(horizon_data, n_simulations = 500)
#' print(summarized_results)
#'
#' @export
simulate_soil_profile_thickness <- function(horizon_data, n_simulations = 500) {
  # Step 1: Infill missing depth variability values
  horizon_data <- infill_missing_depth_variability(horizon_data)

  results_list <- list()

  # Loop over the number of simulations
  for (sim in 1:n_simulations) {
    # Randomly select the simulation direction: top-down or bottom-up
    if (stats::runif(1) < 0.5) {
      # Simulate profile using the top-down method
      result <- simulate_soil_profile_top_down(horizon_data)
      method <- "Top-Down"
    } else {
      # Simulate profile using the bottom-up method
      result <- simulate_soil_profile_bottom_up(horizon_data)
      method <- "Bottom-Up"
    }

    # Calculate thickness for each horizon as the difference between bottom and top depths
    result$Thickness <- result$bottom - result$top

    # Record the simulation method and simulation number
    result$Method <- method
    result$Simulation <- sim

    # Store the current simulation result
    results_list[[sim]] <- result
  }

  # Combine all simulation results into a single dataframe
  all_results <- do.call(rbind, results_list)

  # Reorder columns for clarity: Simulation, Method, hzname, top, bottom, Thickness
  all_results <- all_results[, c("Simulation", "Method", "hzname", "top", "bottom", "Thickness")]

  # Group results by horizon name and calculate the standard deviation of thickness for each horizon
  horizon_sd <- all_results |>
    dplyr::group_by(hzname) |>
    dplyr::summarise(thickness_sd = round(stats::sd(Thickness), 2)) |>
    dplyr::ungroup()

  # Merge the representative top and bottom depths from horizon_data
  summarized_results <- merge(
    horizon_data[, c("hzname", "hzdept_r", "hzdepb_r")],
    horizon_sd,
    by = "hzname"
  )

  # Rename columns for consistency: hzdept_r -> top, hzdepb_r -> bottom
  summarized_results <- summarized_results |>
    dplyr::rename(
      top = hzdept_r,
      bottom = hzdepb_r
    )

  # Return the summarized results
  return(summarized_results)
}


#' Simulate and Perturb Soil Profiles
#'
#' Takes a `SoilProfileCollection` object and perturbs its horizons to
#' simulate variability in soil profile thickness and boundary depths. The
#' simulation is performed in two steps: first, horizon thickness is
#' perturbed based on simulated thickness variability (using a top-down or
#' bottom-up simulation method), and then the boundary depths are perturbed
#' using distinctness-derived offsets. If the profile contains only one
#' horizon, the function simply replicates the profile.
#'
#' @section Known limitation:
#' This function requires `soil_profile`'s horizons to already carry a
#' `sim_comppct` column (it derives the number of simulations to run from
#' `unique(horizons(soil_profile)$sim_comppct)`). `sim_component_comp()`
#' (`R/property-simulation.R`) produces this column, but at component
#' (`cokey`) grain, not horizon grain - callers must
#' `dplyr::left_join(horizon_data, sim_component_comp(component_data), by = "cokey")`
#' before calling this function; it will still error with a missing-column
#' condition if that join hasn't been done first.
#'
#' @param soil_profile A SoilProfileCollection object containing soil profile data.
#'
#' @return A SoilProfileCollection object with perturbed horizon depths.
#'
#' @export
simulate_and_perturb_soil_profiles <- function(soil_profile) {
  # Step 1: Extract horizons from the soil profile and select relevant columns
  soil_profile@horizons <- aqp::horizons(soil_profile) |>
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct,
                  hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, sim_comppct)

  horizon_data <- aqp::horizons(soil_profile)
  horizon_data$genhz <- aqp::generalizeHz(
    horizon_data$hzname,
    new = c('O','A','B','C','Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C','^\\d*Cr','^\\d*R'),
    ordered = TRUE
  )
  n_simulations <- unique(horizon_data$sim_comppct)  # Set number of simulations based on sim_comppct

  # Check if there's only one horizon (e.g., R horizon only)
  if (nrow(horizon_data) == 1) {
    message("Only one horizon in the profile. Skipping perturbation for this profile.")

    # Return the original soil profile replicated n_simulations times (no perturbation needed)
    replicated_profiles <- list()
    for (i in 1:n_simulations) {
      # Create a copy of the soil profile and assign a unique profile ID
      new_profile <- soil_profile
      aqp::profile_id(new_profile) <- paste0(aqp::profile_id(new_profile), "_sim_", i)
      replicated_profiles[[i]] <- new_profile
    }
    simulated_profiles <- aqp::combine(replicated_profiles)
    return(simulated_profiles)
  }

  # Step 1 (continued): Simulate soil profile thickness
  simulated_thickness <- simulate_soil_profile_thickness(horizon_data, n_simulations)
  simulated_thickness$genhz <- aqp::generalizeHz(
    simulated_thickness$hzname,
    new = c('O','A','B','C','Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C','^\\d*Cr','^\\d*R'),
    ordered = TRUE
  )

  # Step 2: Query OSD distinctness and get bound_sd values
  distinctness_data <- query_osd_distinctness(horizon_data)
  bound_lut <- distinctness_data |>
    dplyr::group_by(id, genhz) |>
    dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE)) |>
    dplyr::ungroup()

  # Step 3: Merge the simulated thickness data and distinctness data
  combined_data <- merge(simulated_thickness, bound_lut, by = "genhz") |>
    dplyr::arrange(top)

  # Step 4: Add thickness standard deviation and bound_sd to the soil profile horizons
  horizon_data <- horizon_data |>
    dplyr::left_join(combined_data |> dplyr::select(hzname, thickness_sd, bound_sd) |> dplyr::distinct(), by = "hzname")
  soil_profile@horizons <- horizon_data |>
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct, hzdept_r, hzdepb_r, thickness_sd, bound_sd)

  # Step 3: First perturbation (horizon thickness)
  perturbed_profiles_thickness <- aqp::perturb(
    soil_profile,
    n = n_simulations,
    thickness.attr = "thickness_sd"
  )

  # Step 4: Second perturbation (boundary distinctness)
  list_of_profiles <- list()
  for (i in seq_along(perturbed_profiles_thickness)) {
    single_profile <- perturbed_profiles_thickness[i]
    perturbed_profile <- aqp::perturb(
      single_profile,
      n = 1,
      boundary.attr = "bound_sd",
      id = single_profile@site$id
    )
    list_of_profiles[[i]] <- perturbed_profile
  }

  simulated_profiles <- aqp::combine(list_of_profiles)

  # Step 5: Identify out-of-range horizons and adjust depths
  adjust_out_of_range_profiles <- function(simulated_profiles, horizon_data) {
    sim_horizons <- aqp::horizons(simulated_profiles) |>
      dplyr::select(id, mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

    merged_data <- merge(
      sim_horizons,
      horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
      by = "hzname",
      all.x = TRUE
    )

    merged_data <- merged_data |>
      dplyr::rowwise() |>
      dplyr::mutate(
        top = dplyr::if_else(is.na(hzdept_l), top, pmax(top, hzdept_l)),
        bottom = dplyr::if_else(is.na(hzdepb_h), bottom, pmin(bottom, hzdepb_h))
      ) |>
      dplyr::ungroup()

    simulated_profiles_hzdata <- aqp::horizons(simulated_profiles) |>
      dplyr::left_join(merged_data |> dplyr::select(hzname, id, top, bottom), by = c("id", "hzname"))
    simulated_profiles@horizons <- simulated_profiles_hzdata

    return(simulated_profiles)
  }

  simulated_profiles <- adjust_out_of_range_profiles(simulated_profiles, horizon_data)

  # Step 6: Return the adjusted simulated profiles
  return(simulated_profiles)
}


#' Run Soil Profile Depth Simulations for a SoilProfileCollection
#'
#' Simulates soil profile depths for each profile in a `SoilProfileCollection`
#' by applying `simulate_and_perturb_soil_profiles()` to each individual
#' profile, then combines the results into a single `SoilProfileCollection`.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value for setting the random seed (default is 123) for reproducible simulations.
#'
#' @return A SoilProfileCollection object containing the simulated soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_collection <- simulate_profile_depths_by_collection(soil_collection, seed = 123)
#' }
#'
#' @export
simulate_profile_depths_by_collection <- function(soil_collection, seed = 123) {

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Loop over each profile in the SoilProfileCollection
  for (i in seq_along(soil_collection)) {
    # Extract the current soil profile
    soil_profile <- soil_collection[i, ]

    # Simulate and perturb the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Run Soil Profile Depth Simulations in Parallel for a SoilProfileCollection
#'
#' Simulates soil profile depths for each profile in a `SoilProfileCollection`
#' in parallel using the `future`/`future.apply` packages, then combines the
#' results into a single `SoilProfileCollection`.
#'
#' @section Note on parallel workers:
#' `future::multisession` workers are fresh R processes. `future`/`globals`
#' auto-detect that `simulate_single_profile()` calls a `soilSIM`-namespaced
#' function and attach the package in each worker automatically - this only
#' works when `soilSIM` is actually installed and attached in the calling
#' session, not merely `devtools::load_all()`'d.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value to set the random seed for reproducibility (default is 123).
#' @param n_cores Integer, number of cores to use for parallel processing (default is 6).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_profiles <- simulate_profile_depths_by_collection_parallel(
#'     soil_collection, seed = 123, n_cores = 6
#'   )
#'   print(simulated_profiles)
#' }
#'
#' @export
simulate_profile_depths_by_collection_parallel <- function(soil_collection, seed = 123, n_cores = 6) {

  # Set up the parallel backend using the future package
  future::plan(future::multisession, workers = n_cores)

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Define a function to simulate a single soil profile
  simulate_single_profile <- function(i) {
    soil_profile <- soil_collection[i, ]  # Select the ith profile
    cat("Processing profile ID:", soil_profile$id, "\n")
    simulate_and_perturb_soil_profiles(soil_profile)
  }

  # Run the simulation in parallel using future_lapply
  result <- tryCatch({
    all_simulated_profiles <- future.apply::future_lapply(seq_along(soil_collection), simulate_single_profile)

    # Revert to sequential processing
    future::plan(future::sequential)

    # Combine all simulated profiles into a single SoilProfileCollection
    combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

    return(combined_simulated_profiles)
  }, error = function(e) {
    cat("Error in processing profiles:\n")
    cat("Error message:", e$message, "\n")
    # Optionally return NULL to indicate failure
    return(NULL)
  })

  return(result)
}


#' Simulate Profile Depths by Mukey
#'
#' Queries soil data for a given map unit key (mukey) using the SSURGO
#' database, sets up the corresponding soil profile data, and then simulates
#' and perturbs the soil profiles for that mukey via
#' `simulate_and_perturb_soil_profiles()`.
#'
#' @param mukey A string or numeric value representing the map unit key to query.
#' @param n_simulations Integer, the number of triangular draws per component used to derive
#'   `sim_comppct` via `sim_component_comp()` (default is 100).
#' @param seed An integer to set the random seed for reproducibility (default is 123).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_profiles <- simulate_profile_depths_by_mukey("123456", n_simulations = 100, seed = 123)
#'   print(simulated_profiles)
#' }
#'
#' @export
simulate_profile_depths_by_mukey <- function(mukey, n_simulations = 100, seed = 123) {
  # Step 1: Query the data based on the mukey
  mu_data <- get_aws_data_by_mukey(mukeys = mukey)

  if (is.null(mu_data)) {
    stop("No data found for the provided mukey.")
  }

  # Step 1.5: Derive sim_comppct (sim_component_comp() produces it at cokey grain) and join
  # it onto every horizon row by cokey - simulate_and_perturb_soil_profiles() requires it at
  # horizon grain. Mirrors the already-verified pattern in R/ssurgo-simulation.R's
  # simulate_ssurgo_mapunit_draws() (closes the gap this function used to document as a
  # Known limitation: callers previously had to perform this join themselves).
  component_data <- sim_component_comp(mu_data, n_simulations = n_simulations)
  mu_data <- dplyr::left_join(
    mu_data, component_data[, c("cokey", "sim_comppct")],
    by = "cokey"
  )

  # Step 2: Set up the soil profile data
  set.seed(seed)  # For reproducibility
  # `id` (not `compname`) is the profile-ID column convention the rest of
  # this file assumes (query_osd_distinctness() renames compname -> id
  # itself, and adjust_out_of_range_profiles()/evaluate_simulated_depths()
  # both hardcode a literal `id` column) - derive it before calling
  # aqp::depths<- so downstream horizon selects don't fail with
  # "object 'id' not found".
  mu_data$id <- mu_data$compname
  aqp::depths(mu_data) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(mu_data) <- "hzname"

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 3: Loop over each profile in mu_data
  for (i in seq_along(mu_data)) {
    # Select the current profile
    soil_profile <- mu_data[i, ]

    # Step 4: Simulate and perturb the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Step 5: Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Step 6: Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Evaluate Simulated Profile Depths
#'
#' Evaluates whether the simulated soil profile depths fall outside the
#' specified range. Extracts the simulated top and bottom depths from a
#' `SoilProfileCollection`, merges these with the original horizon data, and
#' flags horizons where the simulated top is less than the lower bound or
#' the simulated bottom exceeds the upper bound.
#'
#' @param simulated_profiles A SoilProfileCollection object containing simulated soil profiles.
#' @param horizon_data A data frame with original horizon data, including columns:
#'        - `hzname`: Horizon name.
#'        - `hzdept_l`: Lower bound for the top depth.
#'        - `hzdepb_h`: Upper bound for the bottom depth.
#'
#' @return A data frame containing rows (horizons) where the simulated depths are out of range.
#'
#' @examples
#' \dontrun{
#'   out_of_range <- evaluate_simulated_depths(simulated_profiles, horizon_data)
#'   head(out_of_range)
#' }
#'
#' @export
evaluate_simulated_depths <- function(simulated_profiles, horizon_data) {
  # Extract horizons from the simulated profiles
  sim_horizons <- aqp::horizons(simulated_profiles) |>
    dplyr::select(mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

  # Merge simulated horizons with original horizon data for comparison
  merged_data <- merge(
    sim_horizons,
    horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
    by = "hzname",
    all.x = TRUE
  )

  # Check if simulated top and bottom depths fall outside the specified range
  merged_data$out_of_range <- with(merged_data, (top < hzdept_l) | (bottom > hzdepb_h))

  # Return rows where simulated depths are out of range
  out_of_range_data <- merged_data[merged_data$out_of_range == TRUE, ]

  return(out_of_range_data)
}
