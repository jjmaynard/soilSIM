#' @title SSURGO Data Acquisition & Management
#' @description Functions for downloading, processing, and validating SSURGO soil survey data
#' @name ssurgo_data
NULL

# ==============================================================================
# CORE DOWNLOAD FUNCTIONS
# ==============================================================================

#' Download SSURGO Tabular Data with Comprehensive Processing
#'
#' Enhanced version of the proven download_ssurgo_tabular function with additional
#' caching, validation, and reporting capabilities while maintaining compatibility
#' with existing soil property simulation workflows.
#'
#' @param aoi_wkt Character. Well-Known Text (WKT) representation of the area of interest
#' @param properties Character vector. Soil properties to download.
#'   Default: c("w3b", "w15b", "db", "cec", "rfv", "clay", "ph", "sand", "silt", "soc")
#' @param include_restrictions Logical. Whether to include horizon restriction data for
#'   unsuitable horizon detection (default: TRUE)
#' @param cache_dir Character. Directory for caching downloaded data (default: NULL for no caching)
#' @param force_download Logical. If TRUE, bypass cache and download fresh data (default: FALSE)
#' @param validate_data Logical. Perform comprehensive data validation (default: TRUE)
#' @param verbose Logical. Provide detailed progress messages (default: FALSE)
#'
#' @return List containing:
#'   \itemize{
#'     \item ssurgo_data: Combined horizon and component data with restriction flags
#'     \item mu: Spatial map unit data (from soilDB::mukey.wcs)
#'     \item metadata: Download and processing metadata
#'     \item validation_results: Data validation results (if validate_data = TRUE)
#'     \item cache_info: Cache usage information (if caching enabled)
#'     \item components_missing_horizons: Data frame (`mukey`/`cokey`/`compname`/`comppct_l/r/h`),
#'       one row per component that has a real `comppct` but zero `chorizon` rows in SDA AND no
#'       AOI sibling (same `compname`) to recover a profile from - empty if none. Always empty for
#'       a cache hit (see Component recovery section below).
#'     \item components_recovered: Data frame (`mukey`/`cokey`/`compname`/`n_sibling_cokeys_used`),
#'       one row per component whose horizon profile WAS synthesized from AOI siblings - empty if
#'       none. Always empty for a cache hit (see Component recovery section below).
#'   }
#'
#' @section Component recovery (on by default):
#' `execute_ssurgo_query_working()`'s `INNER JOIN chorizon` silently drops any real component with
#' zero `chorizon` rows in SDA - a common, verified SSURGO data-completeness gap, especially for
#' minor components. This function recovers such components automatically: for each one, it
#' searches this same AOI's own fetched data for other cokeys sharing the missing component's
#' `compname` that DO have full horizon data, and averages their profiles (aligned by
#' `classify_genhz()` group) to synthesize a representative one - see
#' `synthesize_component_horizons_from_siblings()`'s docs for the exact averaging rule. Synthesized
#' rows are tagged `component_synthesized = TRUE` and `infill_method = "component_aoi_average"` in
#' `ssurgo_data` for traceability. A component with no AOI sibling to recover from is left out of
#' `ssurgo_data` (same as before this feature existed) but reported in
#' `components_missing_horizons` instead of vanishing with zero trace.
#'
#' @section Cache invalidation:
#' A cache entry (`cache_dir`) written before this feature shipped, or before it was extended,
#' predates component recovery entirely - `ssurgo_data` in that cached entry won't contain any
#' synthesized rows, and a cache hit always returns empty `components_missing_horizons`/
#' `components_recovered` placeholders regardless of what a fresh download would find (only counts,
#' not the full data frames, are persisted via `metadata`). Clear the cache to pick up recovered
#' components.
#'
#' @export
download_ssurgo_tabular <- function(aoi_wkt,
                                    properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "ph1to1h2o",
                                                   "cec7", "om", "wthirdbar", "wfifteenbar"),
                                    include_restrictions = TRUE,
                                    cache_dir = NULL,
                                    force_download = FALSE,
                                    validate_data = TRUE,
                                    verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  start_time <- Sys.time()

  if (verbose) {
    log_message("INFO", "SSURGO Data Download Started", category = "Download")
  }

  # Step 1: Input validation using validation function
  if (verbose) log_message("INFO", "Validating inputs", category = "Download")

  # Use the FIXED validation function instead of the problematic one
  #
  # NOTE: the original source passed config_file/log_validation here, but
  # validate_download_inputs_ssurgo() (unlike its
  # validate_download_inputs_ssurgo_with_config() sibling) does not accept
  # those parameters - R CMD check's static "unused arguments" analysis
  # confirmed this call could never actually succeed. Since this code path
  # is network-dependent and untested (no prior coverage could have caught
  # it at runtime), dropping the two invalid arguments is the minimal fix
  # that makes the call match the function it actually calls.
  input_validation <- validate_download_inputs_ssurgo(
    aoi_wkt = aoi_wkt,
    properties = properties,
    include_restrictions = include_restrictions
  )

  if (!input_validation$valid) {
    stop("Input validation failed: ", paste(input_validation$errors, collapse = "; "))
  }

  if (verbose && length(input_validation$warnings) > 0) {
    for (warning_msg in input_validation$warnings) {
      log_message("WARN", warning_msg, category = "Download")
    }
  }

  # Step 2: Check cache if enabled
  cache_result <- NULL
  if (!is.null(cache_dir) && !force_download) {
    if (verbose) log_message("INFO", "Checking cache", category = "Download")

    cache_result <- check_ssurgo_cache(
      aoi_wkt = aoi_wkt,
      properties = properties,
      include_restrictions = include_restrictions,
      cache_dir = cache_dir,
      verbose = verbose
    )

    if (!is.null(cache_result) && cache_result$cache_hit) {
      if (verbose) log_message("INFO", "Cache hit! Loading data from cache", category = "Download")

      # Validate cached data if requested
      validation_results <- NULL
      if (validate_data) {
        validation_results <- validate_data_quality(
          cache_result$data,
          required_columns = c("cokey", "hzdept_r", "hzdepb_r")
        )
      }

      # NOTE: components_missing_horizons/components_recovered are NOT persisted across a cache
      # round-trip (only their counts, via cache_result$metadata's
      # components_missing_horizons_count/components_recovered_count fields, if the cache entry
      # was written by a post-recovery-feature download) - a cache hit returns empty placeholders
      # here rather than extending the on-disk cache schema. A cache entry written BEFORE this
      # feature shipped also predates component recovery entirely (ssurgo_data itself won't
      # contain any synthesized rows) - clear the cache to pick up recovered components.
      return(list(
        ssurgo_data = cache_result$data,
        mu = cache_result$mu,
        metadata = cache_result$metadata,
        validation_results = validation_results,
        cache_info = list(
          cache_hit = TRUE,
          cache_file = cache_result$cache_file,
          cache_timestamp = cache_result$timestamp
        ),
        components_missing_horizons = data.frame(mukey = character(0), cokey = character(0),
          compname = character(0), comppct_l = numeric(0), comppct_r = numeric(0),
          comppct_h = numeric(0), stringsAsFactors = FALSE),
        components_recovered = data.frame(mukey = character(0), cokey = character(0),
          compname = character(0), n_sibling_cokeys_used = integer(0), stringsAsFactors = FALSE)
      ))
    }
  }

  # Step 3: Setup property configuration
  if (verbose) log_message("INFO", "Setting up property configuration", category = "Download")

  ssurgo_lookup <- create_ssurgo_property_lookup_working()

  # Validate requested properties
  valid_properties <- input_validation$validation_metadata$property_stats$valid_property_names
  if (length(valid_properties) == 0) {
    stop("No valid properties specified after validation")
  }

  if (verbose) {
    log_message("INFO", paste("Using", length(valid_properties), "valid properties:",
                              paste(valid_properties, collapse = ", ")), category = "Download")
  }

  # Step 4: Process AOI and get map unit keys
  if (verbose) log_message("INFO", "Processing area of interest", category = "Download")

  spatial_result <- process_aoi_and_get_mukeys_working(aoi_wkt, verbose = verbose)
  aoi <- spatial_result$aoi
  mu <- spatial_result$mu
  mukey_list <- spatial_result$mukey_list

  if (length(mukey_list) == 0) {
    stop("No map unit keys found for the specified AOI")
  }

  if (verbose) {
    log_message("INFO", paste("Found", length(mukey_list), "unique map unit keys"), category = "Download")
  }

  # Step 5: Execute SSURGO query
  if (verbose) log_message("INFO", "Executing SSURGO database query", category = "Download")

  ssurgo_data <- execute_ssurgo_query_working(
    mukey_list = mukey_list,
    properties = properties,
    ssurgo_lookup = ssurgo_lookup,
    include_restrictions = include_restrictions,
    verbose = verbose
  )

  # Step 6: Process rock fragment volume if requested
  if ("rfv" %in% properties && all(c("rfv_l", "rfv_r", "rfv_h") %in% names(ssurgo_data))) {
    if (verbose) log_message("INFO", "Aggregating rock fragment volume", category = "Download")
    ssurgo_data <- aggregate_rock_fragment_volume_working(ssurgo_data, verbose = verbose)
  }

  # Step 7: Add restriction and unsuitable horizon indicators using Module 8
  if (include_restrictions && nrow(ssurgo_data) > 0) {
    if (verbose) log_message("INFO", "Processing restriction and horizon suitability data", category = "Download")
    ssurgo_data <- add_restriction_indicators_working(ssurgo_data, verbose = verbose)
  } else if (!include_restrictions) {
    if (verbose) log_message("INFO", "Adding basic unsuitable horizon detection", category = "Download")
    # Use Module 8 function for unsuitable horizon detection
    ssurgo_data$unsuitable_horizon <- is_unsuitable(ssurgo_data, hzname_col = "hzname")
  }

  # Step 7.5: Recover components with a real comppct but zero chorizon rows - a real, verified
  # SSURGO data-completeness gap (execute_ssurgo_query_working()'s INNER JOIN chorizon silently
  # drops such components with zero trace). On by default: search this same AOI's ssurgo_data for
  # other cokeys sharing the missing component's compname and average their profiles (see
  # recover_missing_horizon_components()'s docs). Runs after RFV/restriction processing (so
  # synthesized rows average already-fully-processed sibling columns) and before validation (so
  # the fuller dataset is what gets validated/cached).
  if (verbose) log_message("INFO", "Checking for components missing horizon data", category = "Download")
  all_components <- fetch_ssurgo_all_components_working(mukey_list, verbose = verbose)
  recovery_result <- recover_missing_horizon_components(ssurgo_data, all_components, verbose = verbose)
  ssurgo_data <- recovery_result$ssurgo_data
  components_missing_horizons <- recovery_result$components_missing_horizons
  components_recovered <- recovery_result$components_recovered

  # Step 8: Data validation using Module 8
  validation_results <- NULL
  if (validate_data) {
    if (verbose) log_message("INFO", "Validating downloaded data", category = "Download")
    validation_results <- validate_data_quality(
      ssurgo_data,
      required_columns = c("cokey", "hzdept_r", "hzdepb_r"),
      numeric_columns = paste0(gsub("_.*", "", properties), "_r")
    )

    if (validation_results$overall_quality$score < 0.7) {
      log_message("WARN", "Downloaded data has quality issues", category = "Download")
    }
  }

  # Step 9: Cache data if enabled
  cache_info <- list(cache_hit = FALSE, cache_enabled = FALSE)
  if (!is.null(cache_dir)) {
    if (verbose) log_message("INFO", "Caching data", category = "Download")
    cache_info <- cache_ssurgo_data(
      data = ssurgo_data,
      mu = mu,
      aoi_wkt = aoi_wkt,
      properties = properties,
      include_restrictions = include_restrictions,
      cache_dir = cache_dir,
      verbose = verbose
    )
  }

  # Step 10: Generate comprehensive metadata
  end_time <- Sys.time()
  metadata <- create_download_metadata(
    start_time = start_time,
    end_time = end_time,
    aoi_wkt = aoi_wkt,
    properties = properties,
    include_restrictions = include_restrictions,
    mukey_count = length(mukey_list),
    data_rows = nrow(ssurgo_data),
    unique_cokeys = length(unique(ssurgo_data$cokey %||% character(0))),
    spatial_result = spatial_result,
    validation_results = validation_results,
    components_missing_horizons = components_missing_horizons,
    components_recovered = components_recovered
  )

  if (verbose) {
    log_message("INFO", "Download Complete", category = "Download")
    log_message("INFO", paste("Duration:", round(metadata$duration_seconds, 2), "seconds"), category = "Download")
    log_message("INFO", paste("Rows downloaded:", metadata$data_rows), category = "Download")
  }

  # Return comprehensive results
  return(list(
    ssurgo_data = ssurgo_data,
    mu = mu,
    metadata = metadata,
    validation_results = validation_results,
    cache_info = cache_info,
    components_missing_horizons = components_missing_horizons,
    components_recovered = components_recovered
  ))
}

# ==============================================================================
# SSURGO-SPECIFIC WORKING CODE INTEGRATION FUNCTIONS
# ==============================================================================

#' Create SSURGO Property Lookup Table (Working Version)
#'
#' Creates lookup table exactly as used in the working soil_infill_functions.R
#'
#' @return Data frame with property mappings for low, representative, and high values
#'
#' @export
create_ssurgo_property_lookup_working <- function() {
  # Replace with:
  data.frame(
    Property = c("sandtotal", "claytotal", "silttotal", "dbovendry", "ph1to1h2o",
                 "cec7", "om", "wthirdbar", "wfifteenbar"),
    SSURGO_Label_Low = c("sandtotal_l", "claytotal_l", "silttotal_l", "dbovendry_l",
                         "ph1to1h2o_l", "cec7_l", "om_l", "wthirdbar_l", "wfifteenbar_l"),
    SSURGO_Label_Rep = c("sandtotal_r", "claytotal_r", "silttotal_r", "dbovendry_r",
                         "ph1to1h2o_r", "cec7_r", "om_r", "wthirdbar_r", "wfifteenbar_r"),
    SSURGO_Label_High = c("sandtotal_h", "claytotal_h", "silttotal_h", "dbovendry_h",
                          "ph1to1h2o_h", "cec7_h", "om_h", "wthirdbar_h", "wfifteenbar_h"),
    stringsAsFactors = FALSE
  )
}

#' Process AOI and Get Map Unit Keys (Working Version)
#'
#' Handles spatial processing exactly as in working code
#'
#' @param aoi_wkt Well-Known Text representation of area of interest
#' @param verbose Logical; provide progress messages
#'
#' @return List with processed AOI, map units, and mukey list
process_aoi_and_get_mukeys_working <- function(aoi_wkt, verbose = FALSE) {

  # Convert AOI WKT to spatial object and reproject (exact working logic)
  tryCatch({
    if (verbose) log_message("DEBUG", "Converting WKT to spatial object", category = "Spatial")

    aoi <- terra::vect(aoi_wkt, crs = "epsg:4326")
    aoi <- terra::project(aoi, "epsg:5070")
    aoi <- terra::as.polygons(terra::ext(aoi), crs = "epsg:5070")

    if (verbose) {
      bbox <- terra::ext(aoi)
      log_message("DEBUG", paste("AOI bounds (EPSG:5070): X:", round(bbox[1]), "to", round(bbox[2]),
                                 ", Y:", round(bbox[3]), "to", round(bbox[4])), category = "Spatial")
    }

  }, error = function(e) {
    stop(paste("Error processing AOI WKT:", e$message))
  })

  # Fetch SSURGO map unit keys (exact working logic)
  tryCatch({
    if (verbose) log_message("DEBUG", "Fetching SSURGO map unit keys", category = "Spatial")

    mu <- soilDB::mukey.wcs(aoi = aoi, db = "gssurgo")
    mukey_list <- unique(terra::values(mu))
    mukey_list <- mukey_list[!is.na(mukey_list)]

    if (verbose) log_message("DEBUG", paste("Retrieved", length(mukey_list), "unique map unit keys"), category = "Spatial")

  }, error = function(e) {
    stop(paste("Error fetching mukeys:", e$message))
  })

  return(list(
    aoi = aoi,
    mu = mu,
    mukey_list = mukey_list
  ))
}

#' Execute SSURGO Query (Working Version)
#'
#' Executes SSURGO database query exactly as in working code
#'
#' @param mukey_list Vector of map unit keys
#' @param properties Vector of properties to retrieve
#' @param ssurgo_lookup Property lookup table
#' @param include_restrictions Logical; include restriction data
#' @param verbose Logical; provide progress messages
#'
#' @return Data frame with SSURGO data
execute_ssurgo_query_working <- function(mukey_list, properties, ssurgo_lookup,
                                         include_restrictions = TRUE, verbose = FALSE) {

  # Select lookup table rows for requested properties
  props <- ssurgo_lookup[ssurgo_lookup$Property %in% properties, ]

  # Format mukeys for SQL query
  formatted_mukey <- paste0("(", paste0("'", as.integer(mukey_list), "'", collapse = ","), ")")

  # Construct property selection
  ssurgo_variables <- c(props$SSURGO_Label_Low, props$SSURGO_Label_Rep, props$SSURGO_Label_High)
  ssurgo_variables <- ssurgo_variables[!is.na(ssurgo_variables) & ssurgo_variables != ""]
  formatted_properties <- paste(ssurgo_variables, collapse = ", ")

  # Add restriction-related fields if requested
  if (include_restrictions) {
    if (verbose) log_message("DEBUG", "Including restriction and horizon suitability data", category = "Query")

    restriction_fields <- ",
      cht.texcl,
      cht.lieutex,
      ch.desgnmaster,
      cr.reskind,
      cr.resdept_l, cr.resdept_r, cr.resdept_h,
      cr.resthk_l, cr.resthk_r, cr.resthk_h"

    restriction_joins <- "
      LEFT JOIN corestrictions AS cr ON co.cokey = cr.cokey
      LEFT JOIN chtexture AS cht ON ch.chkey = cht.chtkey"
  } else {
    restriction_fields <- ""
    restriction_joins <- ""
  }

  # Construct the SQL query
  query <- sprintf(
    "SELECT
      co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
      ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
      chf.fragsize_r,
      %s%s
    FROM mapunit AS mu
      INNER JOIN component AS co ON mu.mukey = co.mukey
      INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
      LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey%s
    WHERE mu.mukey IN %s
    ORDER BY co.cokey, ch.hzdept_r ASC;",
    formatted_properties, restriction_fields, restriction_joins, formatted_mukey
  )

  if (verbose) log_message("DEBUG", "Submitting query to SSURGO database", category = "Query")

  # Execute query
  ssurgo_data <- tryCatch({
    soilDB::SDA_query(query)
  }, error = function(e) {
    log_message("ERROR", paste("Query failed:", e$message), category = "Query")
    log_message("DEBUG", "Failed query:", category = "Query")
    log_message("DEBUG", query, category = "Query")
    stop(paste("Error executing SSURGO query:", e$message))
  })

  if (is.null(ssurgo_data) || !is.data.frame(ssurgo_data) || nrow(ssurgo_data) == 0) {
    log_message("WARN", "No data returned from SSURGO query", category = "Query")
    return(data.frame())
  }

  if (verbose) {
    log_message("DEBUG", paste("Query successful:", nrow(ssurgo_data), "records returned"), category = "Query")
  }

  return(ssurgo_data)
}

#' Fetch All Components for a Set of Map Unit Keys (No Horizon Join)
#'
#' Companion query to `execute_ssurgo_query_working()`. That function's `INNER JOIN chorizon`
#' silently drops any real SSURGO component that has zero `chorizon` (horizon-level property) rows
#' in SDA - a real, verified data-completeness gap, especially common for minor components (e.g. a
#' 3%-share component whose horizon table was never populated). This function queries the
#' `component` table alone, with no `chorizon` join, so callers can detect exactly which components
#' `execute_ssurgo_query_working()`'s result is missing (see `recover_missing_horizon_components()`).
#'
#' @param mukey_list Vector of map unit keys (same input `execute_ssurgo_query_working()` takes).
#' @param verbose Logical; provide progress messages.
#' @return Data frame with one row per component: `mukey`, `cokey`, `compname`, `comppct_l`,
#'   `comppct_r`, `comppct_h`. An empty (0-row) data frame with these columns if the query returns
#'   nothing or fails outright - this is a best-effort enrichment step, not a hard requirement, so
#'   failures are logged and swallowed rather than raised (unlike `execute_ssurgo_query_working()`,
#'   whose own query failure is fatal to the whole download).
#' @keywords internal
fetch_ssurgo_all_components_working <- function(mukey_list, verbose = FALSE) {
  empty_result <- data.frame(mukey = character(0), cokey = character(0), compname = character(0),
                              comppct_l = numeric(0), comppct_r = numeric(0), comppct_h = numeric(0),
                              stringsAsFactors = FALSE)

  formatted_mukey <- paste0("(", paste0("'", as.integer(mukey_list), "'", collapse = ","), ")")

  query <- sprintf(
    "SELECT co.mukey, co.cokey, co.compname, co.comppct_l, co.comppct_r, co.comppct_h
     FROM component AS co
     WHERE co.mukey IN %s;",
    formatted_mukey
  )

  if (verbose) log_message("DEBUG", "Fetching component-only table (no chorizon join)", category = "Query")

  result <- tryCatch(
    soilDB::SDA_query(query),
    error = function(e) {
      log_message("WARN", paste("Component-only query failed:", e$message), category = "Query")
      NULL
    }
  )

  if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
    return(empty_result)
  }
  result
}

#' Synthesize a Representative Horizon Profile for a Component from AOI Siblings
#'
#' Given one component known to exist (via a real `comppct`) but missing all horizon data
#' (`target_component`), and the already fully processed horizon rows of OTHER cokeys within the
#' same `download_ssurgo_tabular()` call's `ssurgo_data` that share `target_component$compname`
#' (`sibling_horizons`), builds a synthetic horizon profile for `target_component$cokey` by
#' averaging sibling horizons within each `classify_genhz()` group.
#'
#' @section Alignment method:
#' Horizons are aligned by `classify_genhz(hzname)` group (`O`/`A`/`E`/`B`/`C`/`Cr`/`R`), not by
#' position or depth bin - robust to sibling profiles having different horizon counts or depths.
#' Sibling rows whose `hzname` doesn't classify (`classify_genhz()` returns `NA`) are dropped from
#' averaging entirely rather than guessed at - verified before choosing this: no depth-or-property-
#' based genhz fallback exists anywhere else in soilSIM to reuse (`classify_genhz()` itself and the
#' separate `aqp::generalizeHz()` usage in `R/depth-simulation.R` are both purely `hzname`-regex
#' matchers).
#'
#' @section Averaging rule per genhz group:
#' \itemize{
#'   \item Numeric columns (depths, thicknesses, every property column present): arithmetic mean
#'     across whichever sibling rows fall in that genhz group, `na.rm = TRUE`.
#'   \item `hzname`: the most frequent `hzname` string among the group's sibling rows (ties broken
#'     alphabetically) - reused verbatim, never invented, so it re-classifies to the same genhz
#'     group deterministically.
#'   \item Other character/logical columns (`texcl`, `desgnmaster`, `unsuitable_horizon`, etc.):
#'     most frequent non-`NA` value in the group; `NA` if every sibling value in the group is `NA`.
#'   \item `mukey`/`cokey`/`compname`/`comppct_l`/`comppct_r`/`comppct_h`: taken from
#'     `target_component`, never averaged from siblings - the missing component's own real identity
#'     and composition percentages are already known and correct.
#'   \item `chkey`: freshly generated (`"synth_<cokey>_<genhz>"`), unique and clearly synthetic.
#' }
#' A genhz group present in only one sibling still produces a row (mean of one = a copy); a group
#' present in zero siblings is simply absent from the output - never fabricated from nothing.
#'
#' @param target_component One-row data frame: `mukey`, `cokey`, `compname`, `comppct_l`,
#'   `comppct_r`, `comppct_h`.
#' @param sibling_horizons Data frame of horizon rows (same column shape as
#'   `download_ssurgo_tabular()`'s `ssurgo_data`) from one or more OTHER cokeys sharing
#'   `target_component$compname`.
#' @return Data frame of synthesized horizon rows for `target_component$cokey`, one row per genhz
#'   group produced, with `infill_method = "component_aoi_average"` and
#'   `component_synthesized = TRUE` on every row. A 0-row data frame if no sibling horizon was
#'   classifiable at all.
#' @keywords internal
synthesize_component_horizons_from_siblings <- function(target_component, sibling_horizons) {
  sibling_horizons$.genhz <- classify_genhz(sibling_horizons$hzname)
  sibling_horizons <- sibling_horizons[!is.na(sibling_horizons$.genhz), , drop = FALSE]
  if (nrow(sibling_horizons) == 0) {
    return(sibling_horizons[0, setdiff(names(sibling_horizons), ".genhz"), drop = FALSE])
  }

  identity_cols <- c("mukey", "cokey", "compname", "comppct_l", "comppct_r", "comppct_h",
                      "chkey", "hzname", ".genhz")
  avg_cols <- setdiff(names(sibling_horizons), identity_cols)
  numeric_cols <- avg_cols[vapply(sibling_horizons[avg_cols], is.numeric, logical(1))]
  other_cols <- setdiff(avg_cols, numeric_cols)

  # table()/names(table(...)) would coerce x to character (corrupting e.g. a logical
  # unsuitable_horizon column into "TRUE"/"FALSE" strings, breaking bind_rows() against the real
  # rows' logical column downstream) - index back into the original x instead, preserving type.
  mode_val <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    ux <- unique(x)
    counts <- vapply(ux, function(v) sum(x == v), integer(1))
    ux[[which.max(counts)]]
  }

  groups <- split(sibling_horizons, sibling_horizons$.genhz)
  rows <- lapply(names(groups), function(g) {
    grp <- groups[[g]]
    out <- target_component[1, c("mukey", "cokey", "compname", "comppct_l", "comppct_r", "comppct_h"), drop = FALSE]
    out$chkey <- paste0("synth_", target_component$cokey[1], "_", g)
    hz_tab <- sort(table(grp$hzname), decreasing = TRUE)
    out$hzname <- names(hz_tab)[1]
    for (col in numeric_cols) out[[col]] <- mean(as.numeric(grp[[col]]), na.rm = TRUE)
    for (col in other_cols) out[[col]] <- mode_val(grp[[col]])
    out
  })

  result <- dplyr::bind_rows(rows)
  result$infill_method <- "component_aoi_average"
  result$component_synthesized <- TRUE
  if ("hzdept_r" %in% names(result)) result <- result[order(result$hzdept_r), , drop = FALSE]
  rownames(result) <- NULL
  result
}

#' Recover SSURGO Components Missing All Horizon Data
#'
#' `execute_ssurgo_query_working()`'s `INNER JOIN chorizon` silently drops any real component that
#' has zero `chorizon` rows in SDA. This function finds such components (via `all_components`, a
#' component-only fetch with no chorizon join - see `fetch_ssurgo_all_components_working()`) and,
#' for each one, tries to synthesize a representative horizon profile by averaging OTHER cokeys in
#' the same AOI's `ssurgo_data` that share the missing component's `compname` (see
#' `synthesize_component_horizons_from_siblings()`). Components with no such AOI sibling are left
#' out of `ssurgo_data` (same as current behavior) but are returned in
#' `components_missing_horizons` instead of vanishing untraced.
#'
#' @param ssurgo_data The AOI's already fully processed horizon-grain data frame (i.e.
#'   `download_ssurgo_tabular()`'s Step 7 output - after RFV aggregation and restriction
#'   indicators, before `validate_data_quality()`).
#' @param all_components Data frame from `fetch_ssurgo_all_components_working()` (or shaped like
#'   it) - every real component for this AOI's `mukey_list`, independent of horizon data
#'   availability.
#' @param verbose Logical; provide progress messages.
#' @return `list(ssurgo_data = <ssurgo_data, with synthesized rows appended if any>,
#'   components_missing_horizons = <data frame: mukey/cokey/compname/comppct_l/r/h for components
#'   neither present in ssurgo_data nor resolvable via an AOI sibling>, components_recovered =
#'   <data frame: mukey/cokey/compname/n_sibling_cokeys_used for components that WERE
#'   synthesized>)`.
#' @keywords internal
recover_missing_horizon_components <- function(ssurgo_data, all_components, verbose = FALSE) {
  empty_missing <- data.frame(mukey = character(0), cokey = character(0), compname = character(0),
                               comppct_l = numeric(0), comppct_r = numeric(0), comppct_h = numeric(0),
                               stringsAsFactors = FALSE)
  empty_recovered <- data.frame(mukey = character(0), cokey = character(0), compname = character(0),
                                 n_sibling_cokeys_used = integer(0), stringsAsFactors = FALSE)

  if (is.null(all_components) || nrow(all_components) == 0) {
    return(list(ssurgo_data = ssurgo_data, components_missing_horizons = empty_missing,
                components_recovered = empty_recovered))
  }

  present_cokeys <- unique(as.character(ssurgo_data$cokey))
  missing <- all_components[!as.character(all_components$cokey) %in% present_cokeys, , drop = FALSE]
  if (nrow(missing) == 0) {
    return(list(ssurgo_data = ssurgo_data, components_missing_horizons = empty_missing,
                components_recovered = empty_recovered))
  }

  if (verbose) {
    log_message("WARN", paste(nrow(missing), "component(s) have comppct but zero chorizon rows",
                               "- attempting AOI-sibling recovery"), category = "ComponentRecovery")
  }

  unresolved_rows <- list()
  recovered_rows <- list()
  synthesized_list <- list()

  for (i in seq_len(nrow(missing))) {
    target <- missing[i, , drop = FALSE]
    siblings <- ssurgo_data[!is.na(ssurgo_data$compname) & ssurgo_data$compname == target$compname &
                               as.character(ssurgo_data$cokey) != as.character(target$cokey), , drop = FALSE]
    if (nrow(siblings) == 0) {
      unresolved_rows[[length(unresolved_rows) + 1]] <- target
      next
    }
    synth <- synthesize_component_horizons_from_siblings(target, siblings)
    if (nrow(synth) == 0) {
      unresolved_rows[[length(unresolved_rows) + 1]] <- target
      next
    }
    synthesized_list[[length(synthesized_list) + 1]] <- synth
    recovered_rows[[length(recovered_rows) + 1]] <- data.frame(
      mukey = target$mukey, cokey = target$cokey, compname = target$compname,
      n_sibling_cokeys_used = length(unique(siblings$cokey)), stringsAsFactors = FALSE
    )
  }

  if (length(synthesized_list) > 0) {
    # Real chkey values from SDA are integer; synthesize_component_horizons_from_siblings()
    # deliberately generates a "synth_<cokey>_<genhz>" STRING chkey (unique, clearly synthetic) -
    # dplyr::bind_rows() refuses to combine an integer column with a character one, so coerce the
    # real rows' chkey to character first. chkey is only ever used as an identifier/grouping key
    # elsewhere in the package (e.g. aggregate_rock_fragment_volume_working()'s group_by(chkey)),
    # never arithmetically, so this is safe.
    if ("chkey" %in% names(ssurgo_data)) ssurgo_data$chkey <- as.character(ssurgo_data$chkey)
    ssurgo_data <- dplyr::bind_rows(ssurgo_data, dplyr::bind_rows(synthesized_list))
    # bind_rows() fills structurally-absent columns with NA, not the real rows' correct defaults -
    # backfill explicitly so ensure_infilling_columns() (R/data-infilling.R) never sees NA here.
    if (!"infill_method" %in% names(ssurgo_data)) ssurgo_data$infill_method <- ""
    ssurgo_data$infill_method[is.na(ssurgo_data$infill_method)] <- ""
    if (!"component_synthesized" %in% names(ssurgo_data)) ssurgo_data$component_synthesized <- FALSE
    ssurgo_data$component_synthesized[is.na(ssurgo_data$component_synthesized)] <- FALSE
  }

  list(
    ssurgo_data = ssurgo_data,
    components_missing_horizons = if (length(unresolved_rows) > 0) dplyr::bind_rows(unresolved_rows) else empty_missing,
    components_recovered = if (length(recovered_rows) > 0) dplyr::bind_rows(recovered_rows) else empty_recovered
  )
}

#' Aggregate Rock Fragment Volume (Working Version)
#'
#' Aggregates rock fragment volume exactly as in working code
#'
#' @param ssurgo_data SSURGO data frame
#' @param verbose Logical; provide progress messages
#'
#' @return Data frame with aggregated RFV values
aggregate_rock_fragment_volume_working <- function(ssurgo_data, verbose = FALSE) {

  # Check if RFV columns exist
  rfv_columns <- c("rfv_l", "rfv_r", "rfv_h")
  available_rfv_columns <- intersect(rfv_columns, names(ssurgo_data))

  if (length(available_rfv_columns) == 0) {
    if (verbose) log_message("DEBUG", "No RFV columns found, skipping RFV aggregation", category = "Processing")
    return(ssurgo_data)
  }

  if (!all(c("chkey", "fragsize_r") %in% names(ssurgo_data))) {
    if (verbose) log_message("DEBUG", "Required columns (chkey, fragsize_r) not found, skipping RFV aggregation", category = "Processing")
    return(ssurgo_data)
  }

  if (verbose) log_message("DEBUG", paste("Aggregating RFV using columns:", paste(available_rfv_columns, collapse = ", ")), category = "Processing")

  # Convert RFV columns to numeric for available columns
  for (col in available_rfv_columns) {
    ssurgo_data[[col]] <- as.numeric(ssurgo_data[[col]])
  }

  # Create aggregation data with only available columns
  aggregation_columns <- c("chkey", "fragsize_r", available_rfv_columns)

  # Aggregate RFV by horizon for available columns
  rfv_sum <- ssurgo_data |>
    dplyr::select(dplyr::all_of(aggregation_columns)) |>
    dplyr::filter(!is.na(fragsize_r)) |>  # Remove rows with missing fragsize_r
    dplyr::group_by(chkey, fragsize_r) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(available_rfv_columns),
                    ~ as.numeric(dplyr::first(.x)),
                    .names = "{.col}"),
      .groups = 'drop'
    ) |>
    dplyr::group_by(chkey) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(available_rfv_columns),
                    ~ as.numeric(sum(.x, na.rm = TRUE)),
                    .names = "{.col}"),
      .groups = 'drop'
    )

  # Remove existing RFV columns and join with aggregated data
  ssurgo_data <- ssurgo_data |>
    dplyr::select(-dplyr::any_of(available_rfv_columns)) |>
    dplyr::left_join(rfv_sum, by = "chkey") |>
    dplyr::distinct()

  if (verbose) {
    log_message("DEBUG", paste("Aggregated RFV for", nrow(rfv_sum), "horizons"), category = "Processing")
    if ("rfv_r" %in% available_rfv_columns) {
      avg_rfv <- mean(rfv_sum$rfv_r, na.rm = TRUE)
      log_message("DEBUG", paste("Average RFV:", round(avg_rfv, 1), "%"), category = "Processing")
    }
  }

  return(ssurgo_data)
}


#' Add Restriction Indicators (Working Version)
#'
#' Adds restriction analysis exactly as in working code plus unsuitable_horizon flag
#'
#' @param ssurgo_data SSURGO data frame with restriction data
#' @param verbose Logical; provide progress messages
#'
#' @return Data frame with restriction indicators and unsuitable horizon flags
add_restriction_indicators_working <- function(ssurgo_data, verbose = FALSE) {

  if (verbose) log_message("DEBUG", "Analyzing restriction data and horizon suitability", category = "Processing")

  # Check what property columns actually exist in the data
  available_property_columns <- names(ssurgo_data)[grepl("_r$", names(ssurgo_data))]

  # Define the complete list of properties we might want to check
  potential_required_properties <- c("dbovendry_r", "rfv_r", "claytotal_r", "ph1to1h2o_r",
                                     "sandtotal_r", "silttotal_r", "cec7_r", "om_r",
                                     "wthirdbar_r", "wfifteenbar_r")

  # Only check properties that actually exist in the data
  required_properties_to_check <- intersect(potential_required_properties, available_property_columns)

  if (verbose) {
    log_message("DEBUG", paste("Available property columns:", length(available_property_columns)), category = "Processing")
    log_message("DEBUG", paste("Checking for missing values in:", paste(required_properties_to_check, collapse = ", ")), category = "Processing")
  }

  # Apply the restriction logic with flexible property checking
  ssurgo_data <- ssurgo_data |>
    dplyr::mutate(
      # Check for missing required properties (only those that exist)
      missing_required_properties = if(length(required_properties_to_check) > 0) {
        rowSums(dplyr::across(dplyr::all_of(required_properties_to_check), is.na)) > 0
      } else {
        FALSE  # If no properties to check, assume none are missing
      },

      # Calculate effective restriction top depth (resdept_l preferred)
      restriction_top = dplyr::coalesce(resdept_l, resdept_r),

      # Horizon is within or below the restriction depth
      hz_below_restriction = !is.na(hzdept_r) & !is.na(restriction_top) & hzdept_r >= restriction_top,

      # Flag for bedrock restrictions
      has_reskind_restriction = !is.na(reskind) &
        grepl("lithic|paralithic|bedrock|fragipan|duripan|petrocalcic|petrogypsic|cemented|indurated|densic", reskind, ignore.case = TRUE) &
        missing_required_properties & hz_below_restriction,

      # Flag for cemented horizons using texture class and in lieu texture
      is_cemented = (!is.na(texcl) & grepl("cemented|indurated", texcl, ignore.case = TRUE)) |
        (!is.na(lieutex) & grepl("cemented|indurated|bedrock|pavement", lieutex, ignore.case = TRUE)),

      # Combine horizon name patterns for restriction detection
      horizon_suggests_restriction = dplyr::case_when(
        is.na(hzname) & is.na(desgnmaster) ~ FALSE,
        # Match master designation R
        desgnmaster == "R" ~ TRUE,
        # Match Cr, Cd, Cx in hzname with or without numeric prefix
        !is.na(hzname) & grepl("(^|[0-9]+)(Cr|Cd|Cx)", hzname) ~ TRUE,
        # Match restrictive descriptors
        !is.na(hzname) & grepl("\\bfragipan\\b|\\bduripan\\b|\\bpetrocalcic\\b|\\bpetrogypsic\\b|\\bcemented\\b|\\bindurated\\b", hzname, ignore.case = TRUE) ~ TRUE,
        # Cemented horizon by suffix
        !is.na(hzname) & grepl("m$", hzname) & nchar(hzname) > 1 ~ TRUE,
        TRUE ~ FALSE
      ),

      # Organic horizon detection
      organic_horizon = dplyr::case_when(
        is.na(hzname) & is.na(desgnmaster) ~ FALSE,
        # Master designation is O
        desgnmaster == "O" ~ TRUE,
        # hzname starts with O or numeric prefix + O
        !is.na(hzname) & grepl("(^|[0-9]+)O", hzname) ~ TRUE,
        TRUE ~ FALSE
      ),

      # Texture-based restriction indicators
      texture_suggests_restriction = dplyr::case_when(
        (!is.na(texcl) & grepl("cemented|indurated", texcl, ignore.case = TRUE)) ~ TRUE,
        (!is.na(lieutex) & grepl("cemented|indurated|bedrock|pavement", lieutex, ignore.case = TRUE)) ~ TRUE,
        TRUE ~ FALSE
      ),

      # Overall restriction flag for easy filtering
      is_potentially_restrictive = has_reskind_restriction |
        is_cemented |
        horizon_suggests_restriction |
        texture_suggests_restriction
    )

  # Use Module 8 unsuitable horizon detection (assuming this function exists)
  ssurgo_data$unsuitable_horizon <- tryCatch({
    is_unsuitable(ssurgo_data, hzname_col = "hzname")
  }, error = function(e) {
    if (verbose) log_message("WARN", paste("is_unsuitable function not available, using basic detection:", e$message), category = "Processing")

    # Fallback unsuitable detection if Module 8 function not available
    (!is.na(ssurgo_data$hzname) & grepl("^O|^R|Cr|Cd|cemented|indurated", ssurgo_data$hzname, ignore.case = TRUE)) |
      (!is.na(ssurgo_data$organic_horizon) & ssurgo_data$organic_horizon) |
      (!is.na(ssurgo_data$is_potentially_restrictive) & ssurgo_data$is_potentially_restrictive)
  })

  # Report on unsuitable horizon detection
  if (verbose) {
    unsuitable_count <- sum(ssurgo_data$unsuitable_horizon, na.rm = TRUE)
    organic_count <- sum(ssurgo_data$organic_horizon, na.rm = TRUE)
    restrictive_count <- sum(ssurgo_data$is_potentially_restrictive, na.rm = TRUE)
    missing_props_count <- sum(ssurgo_data$missing_required_properties, na.rm = TRUE)

    log_message("DEBUG", paste("Unsuitable horizon analysis: Total unsuitable:", unsuitable_count), category = "Processing")
    log_message("DEBUG", paste("Organic horizons:", organic_count), category = "Processing")
    log_message("DEBUG", paste("Restrictive horizons:", restrictive_count), category = "Processing")
    log_message("DEBUG", paste("Horizons with missing properties:", missing_props_count), category = "Processing")
  }

  return(ssurgo_data)
}


# ==============================================================================
# SSURGO-SPECIFIC INPUT VALIDATION
# ==============================================================================

#' Validate Download Inputs (SSURGO-specific)
#'
#' SSURGO-specific input validation that uses Module 8 general validation functions
#' for parameter validation, logging, and error handling, with custom SSURGO-specific
#' business logic for WKT geometry and property validation.
#'
#' @param aoi_wkt Area of interest in WKT format
#' @param properties Vector of property names
#' @param include_restrictions Logical for restriction data inclusion
#' @param strict_geometry Whether to use strict geometry validation
#' @param max_area_deg2 Maximum allowed area in square degrees
#'
#' @return List with validation results including Module 8 metadata
#'
#' @export
validate_download_inputs_ssurgo <- function(aoi_wkt,
                                            properties,
                                            include_restrictions = TRUE,
                                            strict_geometry = TRUE,
                                            max_area_deg2 = 100) {

  log_message("INFO", "Starting SSURGO download input validation", category = "Validation")

  # Define parameter specifications for Module 8 validation
  param_specs <- list(
    aoi_wkt = list(
      type = "character",
      required = TRUE,
      description = "Well-Known Text (WKT) format geometry string"
    ),
    properties = list(
      type = "character",
      required = TRUE,
      description = "Vector of SSURGO property names"
    ),
    include_restrictions = list(
      type = "logical",
      required = FALSE,
      description = "Flag to include restriction data"
    ),
    strict_geometry = list(
      type = "logical",
      required = FALSE,
      description = "Enable strict geometry validation"
    ),
    max_area_deg2 = list(
      type = "numeric",
      required = FALSE,
      range = c(0.001, 1000),
      description = "Maximum allowed area in square degrees"
    )
  )

  # Use Module 8 parameter validation
  log_message("DEBUG", "Validating basic parameter structure", category = "Validation")

  param_validation <- validate_parameters(
    parameters = list(
      aoi_wkt = aoi_wkt,
      properties = properties,
      include_restrictions = include_restrictions,
      strict_geometry = strict_geometry,
      max_area_deg2 = max_area_deg2
    ),
    parameter_specs = param_specs,
    strict_mode = TRUE
  )

  # Initialize validation results with Module 8 output
  validation_results <- list(
    valid = param_validation$valid,
    errors = param_validation$errors,
    warnings = param_validation$warnings,
    corrected_parameters = param_validation$corrected_parameters,
    validation_metadata = list(
      timestamp = Sys.time(),
      module = "SSURGO Download Validation",
      version = "2.0.0"
    )
  )

  # Early return if basic parameter validation failed
  if (!param_validation$valid) {
    log_message("ERROR", paste("Basic parameter validation failed:",
                               paste(param_validation$errors, collapse = "; ")),
                category = "Validation")
    return(validation_results)
  }

  log_message("DEBUG", "Basic parameter validation passed", category = "Validation")

  # Use corrected parameters for further validation
  validated_params <- param_validation$corrected_parameters

  # ============================================================================
  # SSURGO-SPECIFIC VALIDATIONS
  # ============================================================================

  # 1. WKT Geometry Validation with SSURGO-specific parameters
  log_message("DEBUG", "Validating WKT geometry", category = "Validation")

  geometry_validation <- tryCatch({
    validate_wkt_geometry(
      wkt_string = validated_params$aoi_wkt,
      crs = "epsg:4326",                    # SSURGO uses WGS84 geographic coordinates
      validation_context = "geographic",    # Lat/lon coordinate validation
      area_limits = list(
        min = 0.0001,                       # ~11m x 11m at equator (minimum meaningful area)
        max = validated_params$max_area_deg2, # User-configurable maximum (default 100 deg^2)
        max_action = "warn"                 # Large areas trigger warnings, not errors
      ),
      complexity_limits = list(
        max_vertices = 1000,                # Performance threshold for complex polygons
        max_parts = 100                     # Multi-part geometry limit
      ),
      bounds_check = list(
        xmin = -180, xmax = 180,           # Valid longitude range
        ymin = -90, ymax = 90              # Valid latitude range
      ),
      strict_mode = validated_params$strict_geometry,
      return_geometry = FALSE             # Don't return geometry object (not needed for validation)
    )
  }, error = function(e) {
    handle_workflow_error(
      error = e,
      context = "WKT geometry validation",
      recovery_action = "warn"
    )
    list(valid = FALSE, errors = paste("WKT validation error:", e$message), warnings = character(0))
  })

  # Merge geometry validation results
  validation_results$errors <- c(validation_results$errors, geometry_validation$errors)
  validation_results$warnings <- c(validation_results$warnings, geometry_validation$warnings)

  if (!geometry_validation$valid) {
    validation_results$valid <- FALSE
  }

  # Add geometry metadata - adjust for generalized function structure
  if ("geometry_stats" %in% names(geometry_validation)) {
    validation_results$validation_metadata$geometry_stats <- list(
      bbox = geometry_validation$geometry_stats$bbox,
      area_deg2 = geometry_validation$geometry_stats$area_value,
      area_units = geometry_validation$geometry_stats$area_units,
      geometry_type = geometry_validation$geometry_stats$geometry_type,
      n_vertices = geometry_validation$geometry_stats$n_vertices,
      complexity_score = geometry_validation$geometry_stats$complexity_validation$complexity_stats$complexity_score %||% NA
    )
  }

  # 2. SSURGO Property Validation using simple generalized function
  log_message("DEBUG", "Validating SSURGO properties", category = "Validation")

  property_validation <- tryCatch({
    validate_properties_with_synonyms(
      properties = validated_params$properties,
      property_lookup = "ssurgo",              # Use built-in SSURGO property definitions
      strict_mode = validated_params$strict_geometry, # Reuse strict_mode flag
      max_invalid_pct = 50,                    # Allow up to 50% invalid properties
      performance_threshold = 10               # Warn if >10 properties (SSURGO can be slow)
    )
  }, error = function(e) {
    handle_workflow_error(
      error = e,
      context = "SSURGO property validation",
      recovery_action = "warn"
    )
    list(valid = TRUE, errors = character(0), warnings = paste("Property validation warning:", e$message))
  })

  # Merge property validation results
  validation_results$errors <- c(validation_results$errors, property_validation$errors)
  validation_results$warnings <- c(validation_results$warnings, property_validation$warnings)

  if (!property_validation$valid) {
    validation_results$valid <- FALSE
  }

  # Add property metadata if available
  if (!is.null(property_validation$property_stats)) {
    prop_stats <- property_validation$property_stats
    validation_results$validation_metadata$property_stats <- list(
      total_requested = prop_stats$total_requested,
      valid_properties = prop_stats$total_valid %||% prop_stats$valid_properties,
      invalid_properties = prop_stats$invalid_properties,
      validation_method = prop_stats$validation_method,
      synonym_matches = prop_stats$synonym_matches %||% 0,
      invalid_property_names = prop_stats$invalid_property_names
    )
  }

  # Add property metadata
  if ("property_stats" %in% names(property_validation)) {
    validation_results$validation_metadata$property_stats <- property_validation$property_stats
  }

  # 3. Final Validation Status
  validation_results$valid <- validation_results$valid &&
    geometry_validation$valid &&
    property_validation$valid &&
    length(validation_results$errors) == 0

  # Create summary message using Module 8 patterns
  if (validation_results$valid) {
    if (length(validation_results$warnings) > 0) {
      summary_message <- paste("Validation passed with", length(validation_results$warnings), "warnings")
      log_message("WARN", summary_message, category = "Validation")
    } else {
      summary_message <- "Validation passed successfully"
      log_message("INFO", summary_message, category = "Validation")
    }
  } else {
    summary_message <- paste("Validation failed with", length(validation_results$errors), "errors")
    log_message("ERROR", summary_message, category = "Validation")
  }

  validation_results$message <- summary_message

  # Log summary statistics
  log_message("DEBUG", paste("Validation summary - Errors:", length(validation_results$errors),
                             ", Warnings:", length(validation_results$warnings)),
              category = "Validation")

  return(validation_results)
}



#' Enhanced SSURGO Download Input Validation (Wrapper)
#'
#' Convenience wrapper that provides additional configuration options
#' and integrates with Module 8 configuration management
#'
#' @param aoi_wkt Area of interest in WKT format
#' @param properties Vector of property names
#' @param include_restrictions Logical for restriction data inclusion
#' @param config_file Optional configuration file path
#' @param log_validation Whether to enable detailed logging
#'
#' @return Enhanced validation results
#'
#' @export
validate_download_inputs_ssurgo_with_config <- function(aoi_wkt,
                                                     properties,
                                                     include_restrictions = TRUE,
                                                     config_file = NULL,
                                                     log_validation = TRUE) {

  log_message("INFO", "Starting enhanced SSURGO validation with configuration",
              category = "EnhancedValidation")

  # Load configuration using Module 8
  if (!is.null(config_file)) {
    config <- tryCatch({
      load_configuration(config_file, validate_config = TRUE)
    }, error = function(e) {
      log_message("WARN", paste("Config load failed, using defaults:", e$message),
                  category = "EnhancedValidation")
      get_default_configuration("validation")
    })
    validation_config <- config$validation
  } else {
    validation_config <- get_default_configuration("validation")
  }

  # Setup logging if requested
  if (log_validation && !is.null(validation_config$generate_plots) && validation_config$generate_plots) {
    setup_logging(log_level = "DEBUG", include_timestamp = TRUE)
  }

  # Extract validation parameters from config with safe defaults
  strict_geometry <- validation_config$strict_geometry %||% TRUE
  max_area <- validation_config$max_area_deg2 %||% 100

  log_message("INFO", "Starting enhanced SSURGO validation with configuration",
              category = "EnhancedValidation")

  # Initialize validation results
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    corrected_parameters = list(),
    validation_metadata = list(
      timestamp = Sys.time(),
      module = "SSURGO Download Validation",
      version = "2.0.0"
    )
  )

  # 1. Basic parameter validation using Module 8
  param_specs <- list(
    aoi_wkt = list(
      type = "character",
      required = TRUE,
      description = "Well-Known Text (WKT) format geometry string"
    ),
    properties = list(
      type = "character",
      required = TRUE,
      description = "Vector of SSURGO property names"
    ),
    include_restrictions = list(
      type = "logical",
      required = FALSE,
      description = "Flag to include restriction data"
    )
  )

  param_validation <- validate_parameters(
    parameters = list(
      aoi_wkt = aoi_wkt,
      properties = properties,
      include_restrictions = include_restrictions
    ),
    parameter_specs = param_specs,
    strict_mode = TRUE
  )

  if (!param_validation$valid) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(validation_results$errors, param_validation$errors)
    return(validation_results)
  }

  # 2. WKT Geometry Validation
  log_message("DEBUG", "Validating WKT geometry", category = "Validation")

  geometry_validation <- validate_wkt_geometry(
    wkt_string = aoi_wkt,
    crs = "epsg:4326",
    validation_context = "geographic",
    area_limits = list(
      min = 0.0001,
      max = max_area,
      max_action = "warn"
    ),
    complexity_limits = list(
      max_vertices = 1000,
      max_parts = 100
    ),
    bounds_check = list(
      xmin = -180, xmax = 180,
      ymin = -90, ymax = 90
    ),
    strict_mode = strict_geometry,
    return_geometry = FALSE
  )

  # Merge geometry validation results
  validation_results$errors <- c(validation_results$errors, geometry_validation$errors)
  validation_results$warnings <- c(validation_results$warnings, geometry_validation$warnings)

  if (!geometry_validation$valid) {
    validation_results$valid <- FALSE
  }

  # Store geometry metadata
  validation_results$validation_metadata$geometry_stats <- geometry_validation$geometry_stats

  # 3. Property Validation
  log_message("DEBUG", "Validating SSURGO properties", category = "Validation")

  property_validation <- validate_properties_with_synonyms(
    properties = properties,
    property_lookup = "ssurgo",
    strict_mode = strict_geometry,
    max_invalid_pct = 50,
    performance_threshold = 10
  )

  # Merge property validation results
  validation_results$errors <- c(validation_results$errors, property_validation$errors)
  validation_results$warnings <- c(validation_results$warnings, property_validation$warnings)

  if (!property_validation$valid) {
    validation_results$valid <- FALSE
  }

  # Store property metadata
  validation_results$validation_metadata$property_stats <- property_validation$property_stats

  # 4. Final validation status
  validation_results$valid <- validation_results$valid &&
    geometry_validation$valid &&
    property_validation$valid &&
    length(validation_results$errors) == 0

  # Create summary message
  if (validation_results$valid) {
    if (length(validation_results$warnings) > 0) {
      summary_message <- paste("Validation passed with", length(validation_results$warnings), "warnings")
      log_message("WARN", summary_message, category = "Validation")
    } else {
      summary_message <- "Validation passed successfully"
      log_message("INFO", summary_message, category = "Validation")
    }
  } else {
    summary_message <- paste("Validation failed with", length(validation_results$errors), "errors")
    log_message("ERROR", summary_message, category = "Validation")
  }

  validation_results$message <- summary_message

  # Generate enhanced validation report
  if (!is.null(validation_config$generate_plots) && validation_config$generate_plots) {
    validation_results$validation_report <- generate_validation_report_ssurgo(validation_results)
  }

  return(validation_results)
}

#' Generate SSURGO Validation Report
#'
#' Creates a summary report of validation results
#'
#' @param validation_results Validation results from main function
#'
#' @return Validation report
#'
#' @export
generate_validation_report_ssurgo <- function(validation_results) {

  tryCatch({
    report <- list(
      summary = list(
        status = if(validation_results$valid) "PASSED" else "FAILED",
        total_errors = length(validation_results$errors),
        total_warnings = length(validation_results$warnings),
        validation_timestamp = validation_results$validation_metadata$timestamp
      ),

      details = list(
        errors = validation_results$errors,
        warnings = validation_results$warnings
      )
    )

    # Add geometry report with safe access
    if (!is.null(validation_results$validation_metadata$geometry_stats)) {
      geom_stats <- validation_results$validation_metadata$geometry_stats

      # Safe numeric conversion with defaults
      area_deg2 <- tryCatch({
        as.numeric(geom_stats$area_deg2)
      }, error = function(e) {
        log_message("DEBUG", "Could not convert area to numeric, using default", category = "Validation")
        return(1.0)
      })

      bbox <- tryCatch({
        as.numeric(geom_stats$bbox)
      }, error = function(e) {
        log_message("DEBUG", "Could not convert bbox to numeric, using default", category = "Validation")
        return(c(-180, -90, 180, 90))
      })

      report$geometry_summary <- list(
        area_deg2 = round(area_deg2, 4),
        geometry_type = geom_stats$geometry_type %||% "POLYGON",
        bbox_str = paste("(", paste(round(bbox, 4), collapse = ", "), ")")
      )
    }

    # Add property report with safe access
    if (!is.null(validation_results$validation_metadata$property_stats)) {
      prop_stats <- validation_results$validation_metadata$property_stats

      report$property_summary <- list(
        total_requested = prop_stats$total_requested %||% 0,
        valid_count = prop_stats$valid_properties %||% 0,
        invalid_count = prop_stats$invalid_properties %||% 0,
        success_rate = round((prop_stats$valid_properties %||% 0) /
                               max(1, prop_stats$total_requested %||% 1) * 100, 1)
      )
    }

    log_message("INFO", paste("Validation report generated - Status:", report$summary$status),
                category = "ReportGeneration")

    return(report)

  }, error = function(e) {
    log_message("ERROR", paste("Report generation error:", e$message), category = "ReportGeneration")

    # Return minimal report on error
    return(list(
      summary = list(
        status = "ERROR",
        total_errors = 1,
        total_warnings = 0,
        validation_timestamp = Sys.time()
      ),
      details = list(
        errors = paste("Report generation failed:", e$message),
        warnings = character(0)
      )
    ))
  })
}

# ==============================================================================
# CACHING FUNCTIONS (SSURGO-SPECIFIC)
# ==============================================================================

#' Check SSURGO Data Cache (Enhanced)
#'
#' Enhanced cache checking with metadata validation.
#'
#' @param aoi_wkt Area of interest WKT
#' @param properties Vector of properties
#' @param include_restrictions Logical for restriction data
#' @param cache_dir Cache directory path
#' @param max_age_days Maximum age of cached data in days
#' @param verbose Logical; provide progress messages
#'
#' @return List with cache check results or NULL if no valid cache
#'
#' @export
check_ssurgo_cache <- function(aoi_wkt, properties, include_restrictions,
                                        cache_dir, max_age_days = 30, verbose = FALSE) {

  if (!dir.exists(cache_dir)) {
    if (verbose) log_message("DEBUG", paste("Cache directory does not exist:", cache_dir), category = "Cache")
    return(NULL)
  }

  # Generate cache key
  cache_key <- generate_ssurgo_cache_key(aoi_wkt, properties, include_restrictions)
  cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))

  if (!file.exists(cache_file)) {
    if (verbose) log_message("DEBUG", "No cache file found for this request", category = "Cache")
    return(NULL)
  }

  # Check cache age
  cache_age_days <- as.numeric(difftime(Sys.time(), file.mtime(cache_file), units = "days"))

  if (cache_age_days > max_age_days) {
    if (verbose) {
      log_message("DEBUG", paste("Cache file is too old (", round(cache_age_days, 1),
                                 " days >", max_age_days, "days)"), category = "Cache")
    }
    return(NULL)
  }

  # Try to load cached data with validation
  tryCatch({
    cached_result <- readRDS(cache_file)

    # Validate cached data structure
    if (!all(c("data", "mu", "metadata") %in% names(cached_result))) {
      if (verbose) log_message("DEBUG", "Invalid cache file structure", category = "Cache")
      return(NULL)
    }

    if (verbose) {
      log_message("DEBUG", paste("Found valid cache file (", round(cache_age_days, 1), "days old)"), category = "Cache")
      log_message("DEBUG", paste("Cached data:", nrow(cached_result$data), "rows"), category = "Cache")
    }

    return(list(
      cache_hit = TRUE,
      cache_file = cache_file,
      timestamp = file.mtime(cache_file),
      age_days = cache_age_days,
      data = cached_result$data,
      mu = cached_result$mu,
      metadata = cached_result$metadata
    ))

  }, error = function(e) {
    if (verbose) log_message("WARN", paste("Failed to load cache file:", e$message), category = "Cache")
    return(NULL)
  })
}

#' Cache SSURGO Data (Enhanced)
#'
#' Enhanced data caching with comprehensive metadata.
#'
#' @param data Downloaded SSURGO data
#' @param mu Spatial map unit data
#' @param aoi_wkt Original AOI WKT
#' @param properties Requested properties
#' @param include_restrictions Whether restrictions were included
#' @param cache_dir Cache directory path
#' @param compress Logical; compress cached data
#' @param verbose Logical; provide progress messages
#'
#' @return List with cache operation results
#'
#' @export
cache_ssurgo_data <- function(data, mu, aoi_wkt, properties, include_restrictions,
                                       cache_dir, compress = TRUE, verbose = FALSE) {

  # Create cache directory if needed using Module 8 utilities
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    if (verbose) log_message("DEBUG", paste("Created cache directory:", cache_dir), category = "Cache")
  }

  # Generate cache key and file path
  cache_key <- generate_ssurgo_cache_key(aoi_wkt, properties, include_restrictions)
  cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))

  # Prepare comprehensive cache data
  cache_data <- list(
    data = data,
    mu = mu,
    metadata = list(
      cache_timestamp = Sys.time(),
      aoi_wkt = aoi_wkt,
      properties = properties,
      include_restrictions = include_restrictions,
      data_rows = nrow(data),
      unique_cokeys = length(unique(data$cokey %||% character(0))),
      r_version = R.version.string
    )
  )

  # Save to cache with error handling
  tryCatch({
    saveRDS(cache_data, file = cache_file, compress = compress)

    if (verbose) {
      file_size_mb <- round(file.size(cache_file) / 1024^2, 2)
      log_message("DEBUG", paste("Cached data to:", basename(cache_file), "(", file_size_mb, "MB)"), category = "Cache")
    }

    return(list(
      cached = TRUE,
      cache_file = cache_file,
      cache_key = cache_key,
      file_size_bytes = file.size(cache_file)
    ))

  }, error = function(e) {
    log_message("WARN", paste("Failed to cache data:", e$message), category = "Cache")
    return(list(cached = FALSE, error = e$message))
  })
}

# ==============================================================================
# UTILITY FUNCTIONS (SSURGO-SPECIFIC)
# ==============================================================================

#' Generate SSURGO Cache Key
#'
#' Creates a unique cache key based on download parameters.
#'
#' @param aoi_wkt Area of interest WKT
#' @param properties Vector of properties
#' @param include_restrictions Logical for restriction data
#'
#' @return Character string cache key
generate_ssurgo_cache_key <- function(aoi_wkt, properties, include_restrictions) {

  # Create hash components
  aoi_hash <- digest::digest(aoi_wkt)
  props_hash <- digest::digest(sort(properties))
  restrictions_flag <- ifelse(include_restrictions, "with_restrictions", "no_restrictions")

  # Combine into cache key
  cache_key <- paste0("ssurgo_", substr(aoi_hash, 1, 8), "_",
                      substr(props_hash, 1, 8), "_", restrictions_flag)

  # Ensure valid filename
  cache_key <- gsub("[^A-Za-z0-9_-]", "_", cache_key)

  return(cache_key)
}

#' Create Download Metadata
#'
#' Creates comprehensive metadata for download operations.
#'
#' @param start_time Download start time
#' @param end_time Download end time
#' @param aoi_wkt Area of interest WKT
#' @param properties Requested properties
#' @param include_restrictions Whether restrictions were included
#' @param mukey_count Number of map unit keys
#' @param data_rows Number of data rows returned
#' @param unique_cokeys Number of unique component keys
#' @param spatial_result Spatial processing results
#' @param validation_results Data validation results
#' @param components_missing_horizons Data frame or `NULL` - components with a real `comppct` but
#'   no AOI sibling to recover a profile from (see `recover_missing_horizon_components()`).
#' @param components_recovered Data frame or `NULL` - components successfully synthesized from an
#'   AOI sibling's averaged profile (see `recover_missing_horizon_components()`).
#'
#' @return List with comprehensive metadata
create_download_metadata <- function(start_time, end_time, aoi_wkt, properties, include_restrictions,
                                     mukey_count, data_rows, unique_cokeys, spatial_result, validation_results,
                                     components_missing_horizons = NULL, components_recovered = NULL) {

  list(
    # Timing information
    download_start = start_time,
    download_end = end_time,
    duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),

    # Request parameters
    aoi_wkt = aoi_wkt,
    properties_requested = properties,
    include_restrictions = include_restrictions,

    # Data summary
    mukey_count = mukey_count,
    data_rows = data_rows,
    unique_cokeys = unique_cokeys,
    components_missing_horizons_count = if (!is.null(components_missing_horizons)) nrow(components_missing_horizons) else NA_integer_,
    components_recovered_count = if (!is.null(components_recovered)) nrow(components_recovered) else NA_integer_,

    # Spatial information
    aoi_area_m2 = if (!is.null(spatial_result$aoi)) {
      as.numeric(terra::expanse(spatial_result$aoi))
    } else {
      NA
    },

    # Data quality
    validation_passed = if (!is.null(validation_results)) {
      validation_results$overall_quality$score > 0.7
    } else {
      NA
    },

    # System information
    download_timestamp = Sys.time(),
    r_version = R.version.string,
    data_source = "NRCS SSURGO via Soil Data Access (SDA)"
  )
}

#' Download and Prepare SSURGO Data (Workflow Convenience Wrapper)
#'
#' Thin convenience wrapper around [download_ssurgo_tabular()] that always
#' requests restriction data and validation, then filters the result to
#' `max_depth` - a one-call entry point for the common "download, ready for
#' simulation" path.
#'
#' @param aoi_wkt Area of interest in WKT format.
#' @param properties Character vector of SSURGO properties to download.
#' @param max_depth Maximum depth (cm); horizons deeper than this are dropped.
#'   Pass `NULL` to skip depth filtering.
#' @param cache_dir Optional cache directory (see [download_ssurgo_tabular()]).
#' @param verbose Logical; provide progress messages.
#'
#' @return The same list [download_ssurgo_tabular()] returns, with `ssurgo_data`
#'   depth-filtered and a `preparation_metadata` element added.
#'
#' @export
download_and_prepare_ssurgo <- function(aoi_wkt,
                                        properties = c("clay", "sand", "silt", "db", "ph", "cec", "rfv", "w3b", "w15b"),
                                        max_depth = 250,
                                        cache_dir = NULL,
                                        verbose = FALSE) {

  if (verbose) log_message("INFO", "SSURGO DOWNLOAD AND PREPARATION", category = "Workflow")

  # Download comprehensive data
  download_result <- download_ssurgo_tabular(
    aoi_wkt = aoi_wkt,
    properties = properties,
    include_restrictions = TRUE,
    cache_dir = cache_dir,
    validate_data = TRUE,
    verbose = verbose
  )

  # Filter to maximum depth using Module 8 utilities
  if (!is.null(max_depth) && "hzdepb_r" %in% names(download_result$ssurgo_data)) {
    initial_rows <- nrow(download_result$ssurgo_data)

    download_result$ssurgo_data <- download_result$ssurgo_data |>
      dplyr::filter(hzdepb_r <= max_depth)

    final_rows <- nrow(download_result$ssurgo_data)

    if (verbose) {
      log_message("INFO", paste("Filtered to max depth", max_depth, "cm:", initial_rows, "->", final_rows, "rows"), category = "Workflow")
    }
  }

  # Add preparation metadata
  download_result$preparation_metadata <- list(
    max_depth_applied = max_depth,
    preparation_timestamp = Sys.time(),
    ready_for_simulation = TRUE,
    compatible_with_infill_functions = TRUE
  )

  return(download_result)
}

