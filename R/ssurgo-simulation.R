#' @title Raster SSURGO Percentile Prior: Mukey Lookup and Monte Carlo Draws
#'
#' @description The SSURGO half of the raster fusion prior/likelihood pipeline (see
#'   `R/raster-fusion.R`): given a `terra::SpatVector` AOI, rasterizes SSURGO map units,
#'   Monte Carlo-simulates per-component soil properties (reusing already-ported
#'   `sim_component_comp()`/`simulate_cokey_generalized()` from `R/property-simulation.R`),
#'   aggregates to a requested depth window, and rasterizes per-mukey percentiles for
#'   `fuse_property_adaptive()` to consume. Ported from
#'   `code_ref/reanalysis-platform/{ssurgo_simulation.R, ssurgo_prior.R}`.
#'
#'   Reuses `R/kssl-reference-correlations.R`'s already-built-in, already-shipped
#'   `.kssl_property_matrices()`/`.kssl_texture_matrices()` genhz-keyed correlation matrices
#'   directly as `simulate_cokey_generalized()`'s `correlation_matrices`/
#'   `txt_correlation_matrices` arguments, rather than re-reading the raw
#'   `data/global_cor_matrices.rds`/`data/global_cor_texture_matrices.rds` files at runtime the
#'   way the source's `load_global_correlation_matrices()` did (those files already back
#'   `R/sysdata.rda`, built once via `data-raw/build_kssl_reference_correlations.R` - no new
#'   asset-loading code needed). Similarly reuses `classify_genhz()` (already exported from
#'   `kssl-reference-correlations.R`) in place of the source's `aqp::generalizeHz()` call.
#'
#'   `download_ssurgo_tabular()`'s (`R/ssurgo-acquisition.R`) default `properties` and always-included
#'   base columns (`comppct_l/r/h`, all `_l/_r/_h` horizon triplets) already match
#'   `sim_component_comp()`/`simulate_cokey_generalized()`'s expected SSURGO-stem input vocabulary
#'   directly - the only translation table genuinely needed is `property_to_sim_column()`, mapping
#'   a caller-facing property id to `simulate_cokey_generalized()`'s short-code *output* column
#'   name.
#' @name ssurgo_simulation
NULL

#' Fetch a Raster of SSURGO Map Unit Keys for an AOI
#'
#' @param aoi_vect A `terra::SpatVector` (projected, e.g. EPSG:5070) - the AOI footprint itself,
#'   not widened to its bounding box (unlike `R/ssurgo-acquisition.R`'s
#'   `process_aoi_and_get_mukeys_working()`, which bbox-widens for its tabular by-mukey-list
#'   query - inappropriate here, where over-fetching/over-rasterizing the AOI's full bbox would
#'   waste both a `mukey.wcs()` grid and a `SDA_spatialQuery()` polygon fetch).
#' @return A single-layer, categorical (factor) `terra::SpatRaster` named `"mukey"`, or `NULL` if
#'   either the grid or the polygon fetch returns nothing (including when `soilDB::SDA_spatialQuery()`
#'   errors outright, e.g. because its `sf` dependency isn't loadable in this session - a real,
#'   encountered failure mode, not hypothetical).
#' @export
fetch_ssurgo_mukey_raster <- function(aoi_vect) {
  mu <- soilDB::mukey.wcs(aoi = aoi_vect, db = "gssurgo")
  if (is.null(mu)) return(NULL)

  ssurgo_p <- tryCatch(
    soilDB::SDA_spatialQuery(aoi_vect, what = "mupolygon", db = "SSURGO", geomIntersection = TRUE),
    error = function(e) {
      warning(sprintf("fetch_ssurgo_mukey_raster(): SDA_spatialQuery() failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(ssurgo_p) || nrow(ssurgo_p) == 0) return(NULL)
  ssurgo_p <- terra::project(ssurgo_p, terra::crs(mu))

  mukey_raster <- terra::rasterize(ssurgo_p, mu, field = "mukey")
  names(mukey_raster) <- "mukey"
  terra::as.factor(mukey_raster)
}

#' Infill Missing Soil Property Values
#'
#' Loops `infill_soil_property()` (`R/data-infilling.R`, which already special-cases `"rfv"`
#' internally) over the standard SSURGO property set.
#'
#' @param df A horizon data frame, as returned by `download_ssurgo_tabular()`.
#' @return `df` with missing values infilled where possible.
#' @export
infill_soil_data <- function(df) {
  properties <- c("sandtotal", "claytotal", "silttotal", "dbovendry", "wthirdbar",
                   "wfifteenbar", "ph1to1h2o", "cec7", "om")
  for (p in properties) {
    if (paste0(p, "_r") %in% names(df)) {
      df <- infill_soil_property(df, p)
    }
  }
  if ("rfv_r" %in% names(df)) {
    df <- infill_soil_property(df, "rfv")
  }
  df
}

#' Depth-Trend GP Adjustment for a Single Cokey's Data
#'
#' The per-cokey unit of work `maybe_adjust_soil_data_depth_trend()` maps over every cokey,
#' factored out so it can be dispatched to parallel workers unchanged.
#'
#' @param cokey_data Simulation data for a single cokey.
#' @param properties Character vector of property column names to adjust.
#' @param min_depths Minimum distinct depths required to attempt GP fitting.
#' @return `cokey_data`, depth-trend-adjusted where possible.
#' @keywords internal
adjust_one_cokey_depth_trend <- function(cokey_data, properties, min_depths) {
  n_depths <- length(unique(cokey_data$hzdept_r[!is.na(cokey_data$hzdept_r)]))
  if (n_depths < min_depths) {
    return(cokey_data)
  }
  apply_local_gp_adjustments(cokey_data, properties = properties, min_depths = min_depths)
}

#' Depth-Trend GP Adjustment, Guarded by `GPfit` Availability
#'
#' Applies `apply_local_gp_adjustments()` (`R/multivariate-adjustment.R` - fits its own local GP
#' per property from each cokey's own within-simulation depth trend, no pre-supplied GP models
#' needed) per cokey, when `GPfit` is installed and a cokey has enough distinct depths. Cokeys
#' with fewer than `min_depths` distinct depths pass through unadjusted, exactly as the original
#' per-cokey guard did.
#'
#' Each cokey's GP fitting is completely independent of every other cokey's, so this step is
#' embarrassingly parallel - profiling on a real AOI showed it as the dominant cost of the whole
#' SSURGO simulation pipeline (see `apply_local_gp_adjustments()`/`fit_local_gp_model_single()`),
#' so for AOIs with many cokeys, `parallel = TRUE` can give a further speedup roughly proportional
#' to available cores on top of the sequential-path optimizations already applied there.
#'
#' @param sim_long Long-format simulated data with `cokey`, `hzdept_r`, and property columns.
#' @param properties Character vector of property column names to adjust.
#' @param min_depths Minimum distinct depths required to attempt GP fitting (default 2, matching
#'   the source's `length(unique_depths) >= 2` guard).
#' @param parallel Logical; if `TRUE`, process cokeys across multiple worker processes via the
#'   `parallel` package (default `FALSE` - sequential, matching prior behavior exactly). Mirrors
#'   the Windows-cluster/`mclapply` pattern already used by `multivariate-adjustment.R`'s
#'   `process_cokeys_parallel()` - falls back to sequential processing if the parallel setup
#'   itself errors.
#' @param n_cores Number of worker processes to use when `parallel = TRUE` (default
#'   `max(1, parallel::detectCores() - 1)`).
#' @return `sim_long`, depth-trend-adjusted where possible.
#' @export
maybe_adjust_soil_data_depth_trend <- function(sim_long, properties, min_depths = 2,
                                                parallel = FALSE, n_cores = NULL) {
  if (!requireNamespace("GPfit", quietly = TRUE)) {
    warning("GPfit not installed - skipping depth-trend GP adjustment.")
    return(sim_long)
  }

  properties <- intersect(properties, names(sim_long))
  if (length(properties) == 0) {
    return(sim_long)
  }

  if (!parallel) {
    return(
      sim_long |>
        dplyr::group_by(cokey) |>
        dplyr::group_modify(function(cokey_data, ...) {
          adjust_one_cokey_depth_trend(cokey_data, properties, min_depths)
        }) |>
        dplyr::ungroup()
    )
  }

  cokey_groups <- split(sim_long, sim_long$cokey)
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  adjusted <- tryCatch({
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))

      # Load the package (and therefore all its Imports, including dplyr and GPfit) on each
      # worker, and propagate the current logging option - mirrors the same fix already applied
      # to monte-carlo.R's run_parallel_simulation() and multivariate-adjustment.R's
      # process_cokeys_parallel().
      parallel::clusterEvalQ(cl, library(soilSIM))
      current_log_cfg <- getOption("soil_workflow_log_config")
      parallel::clusterCall(cl, function(cfg) options(soil_workflow_log_config = cfg), current_log_cfg)

      parallel::parLapply(cl, cokey_groups, adjust_one_cokey_depth_trend,
                           properties = properties, min_depths = min_depths)
    } else {
      parallel::mclapply(cokey_groups, adjust_one_cokey_depth_trend,
                          properties = properties, min_depths = min_depths, mc.cores = n_cores)
    }
  }, error = function(e) {
    handle_workflow_error(e, "Parallel depth-trend adjustment", "warn")
    lapply(cokey_groups, adjust_one_cokey_depth_trend, properties = properties, min_depths = min_depths)
  })

  dplyr::bind_rows(adjusted)
}

#' Thickness-Weighted Mean of Simulated Properties over a Depth Window
#'
#' For each replicate (a unique combination of `replicate_cols`), averages every property column
#' across that replicate's horizons, weighted by how much of each horizon's thickness overlaps
#' `[top_depth, bottom_depth]`. Horizons with zero overlap contribute nothing.
#'
#' @param sim_long Long-format simulated data with `hzdept_r`/`hzdepb_r` and `replicate_cols`.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param property_cols Character vector of property columns to aggregate.
#' @param replicate_cols Columns identifying one replicate (default `c("mukey", "cokey",
#'   "simulation_number")`).
#' @return One row per replicate, with `top`/`bottom` set to the requested window and each
#'   property column replaced by its overlap-weighted mean (`NA` if the replicate has no overlap
#'   with the window).
#' @export
aggregate_depth_window_by_replicate <- function(sim_long, top_depth, bottom_depth, property_cols,
                                                 replicate_cols = c("mukey", "cokey", "simulation_number")) {
  overlap <- pmax(0, pmin(sim_long$hzdepb_r, bottom_depth) - pmax(sim_long$hzdept_r, top_depth))
  sim_long <- sim_long[overlap > 0, , drop = FALSE]
  overlap <- overlap[overlap > 0]

  if (nrow(sim_long) == 0) {
    return(sim_long[0, ])
  }

  sim_long$.weight <- overlap

  sim_long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(replicate_cols))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(property_cols), ~ stats::weighted.mean(.x, w = .weight, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(top = top_depth, bottom = bottom_depth)
}

#' Map a Property Id to its Simulated-Data Column Name
#'
#' @param property_id One of `"ph"`, `"bulk_density"`, `"soc"`, `"cec"`, `"clay"`, `"sand"`,
#'   `"silt"`, `"rock_fragments"`.
#' @return The corresponding column name in `simulate_cokey_generalized()`'s output.
#' @export
property_to_sim_column <- function(property_id) {
  mapping <- c(
    ph = "ph", bulk_density = "db", soc = "soc", cec = "cec",
    clay = "clay_total", sand = "sand_total", silt = "silt_total",
    rock_fragments = "rfv"
  )
  if (!property_id %in% names(mapping)) {
    stop(sprintf("property_to_sim_column(): no simulated column mapping for property '%s'.", property_id))
  }
  unname(mapping[property_id])
}

#' Monte Carlo-Simulate SSURGO Property Draws for an AOI
#'
#' Orchestrates the full SSURGO percentile-prior pipeline for one AOI/depth-window: fetch (cached)
#' tabular SSURGO data, infill missing values, derive `genhz`, simulate component composition
#' (`sim_component_comp()`) and join it onto horizons, simulate correlated properties per cokey
#' (`simulate_cokey_generalized()`, using the KSSL reference correlation matrices), optionally
#' remove organic horizons and apply depth-trend GP adjustment, then aggregate to the requested
#' depth window.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param n_mc Number of triangular draws `sim_component_comp()` uses per component (default 1000).
#' @param parallel,n_cores Passed through to `maybe_adjust_soil_data_depth_trend()`'s
#'   `parallel`/`n_cores` - the depth-trend GP adjustment step is this function's dominant cost
#'   for AOIs with many cokeys, and each cokey's GP fitting is independent of every other cokey's.
#'   Default `parallel = FALSE` matches prior behavior exactly.
#' @return A data frame, one row per `mukey`/`cokey`/`simulation_number` replicate, with simulated
#'   property columns aggregated over the depth window - or `NULL` if the tabular fetch fails.
#' @export
simulate_ssurgo_mapunit_draws <- function(aoi_vect, top_depth, bottom_depth, n_mc = 1000,
                                           parallel = FALSE, n_cores = NULL) {
  # download_ssurgo_tabular()/process_aoi_and_get_mukeys_working() (R/ssurgo-acquisition.R)
  # always assume their aoi_wkt argument is lon/lat EPSG:4326 (hardcoded there), regardless of
  # aoi_vect's actual CRS - extracting WKT directly from an already-projected aoi_vect (e.g.
  # EPSG:5070, as this pipeline's other functions expect) would silently mislabel projected
  # meter coordinates as degrees, producing an invalid request. Reproject to EPSG:4326 first.
  aoi_wkt <- terra::geom(terra::project(aoi_vect, "epsg:4326"), wkt = TRUE)[1]

  cache_key <- build_cache_key(aoi_vect, "ssurgo_tabular", top_depth, bottom_depth, "ssurgo_tabular")
  hz_data <- cache_get(cache_key)

  if (is.null(hz_data)) {
    # download_ssurgo_tabular() returns list(ssurgo_data=, mu=, metadata=, ...), not a bare
    # data frame - unwrap it here.
    download_result <- tryCatch(
      download_ssurgo_tabular(aoi_wkt, cache_dir = NULL, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(download_result) || is.null(download_result$ssurgo_data) || nrow(download_result$ssurgo_data) == 0) {
      return(NULL)
    }
    hz_data <- download_result$ssurgo_data
    cache_set(cache_key, "ssurgo_tabular", hz_data)
  }

  hz_data <- infill_soil_data(hz_data)

  if (any(grepl("O", hz_data$hzname))) {
    hz_data <- remove_organic_layer(hz_data)
  }

  hz_data$genhz <- classify_genhz(hz_data$hzname)

  component_data <- sim_component_comp(hz_data, n_simulations = n_mc)
  hz_data <- dplyr::left_join(
    hz_data, component_data[, c("cokey", "sim_comppct")],
    by = "cokey"
  )

  property_matrices <- .kssl_property_matrices()
  texture_matrices <- .kssl_texture_matrices()

  sim_results <- lapply(unique(hz_data$cokey), function(ck) {
    sim_cokey <- hz_data[hz_data$cokey == ck, , drop = FALSE]
    tryCatch(
      simulate_cokey_generalized(sim_cokey, property_matrices, texture_matrices),
      error = function(e) {
        message("simulate_ssurgo_mapunit_draws(): skipping cokey ", ck, ": ", e$message)
        NULL
      }
    )
  })
  # dplyr::bind_rows() (not base do.call(rbind, ...)) - one cokey can come back with
  # texture columns (sand_total/silt_total/clay_total) while another lacks them entirely
  # (e.g. every one of its rows had a texture triplet with a leftover NA - see
  # simulate_cokey_generalized()'s own per-row texture tryCatch()). Base rbind() requires
  # identical columns across all inputs and would otherwise drop the WHOLE mukey's data by
  # erroring here. bind_rows() unions columns, filling NA where a given cokey has none.
  sim_long <- dplyr::bind_rows(sim_results)

  if (is.null(sim_long) || nrow(sim_long) == 0) {
    return(NULL)
  }

  property_cols <- intersect(SSURGO_SIM_PROPERTY_COLUMNS, names(sim_long))
  sim_long <- maybe_adjust_soil_data_depth_trend(sim_long, property_cols,
                                                  parallel = parallel, n_cores = n_cores)

  aggregate_depth_window_by_replicate(
    sim_long, top_depth, bottom_depth, property_cols,
    replicate_cols = c("mukey", "cokey", "simulation_number")
  )
}

#' @rdname simulate_ssurgo_mapunit_draws
#' @export
SSURGO_SIM_PROPERTY_COLUMNS <- c("db", "wr_3b", "wr_15b", "rfv", "ph", "cec", "soc",
                                  "sand_total", "silt_total", "clay_total")

#' Merge Per-Mukey Percentile Values onto a Mukey Raster
#'
#' @param mukey_raster A categorical (factor) mukey `terra::SpatRaster` from
#'   `fetch_ssurgo_mukey_raster()`.
#' @param percentile_by_mukey A data frame with a `mukey` column and one column per percentile
#'   (e.g. `P05`, `P25`, `P50`, `P75`, `P95`).
#' @return A named list of single-layer `terra::SpatRaster`s, one per percentile column.
#' @export
rasterize_mukey_percentiles <- function(mukey_raster, percentile_by_mukey) {
  # terra::cats()'s first column is the raster's internal numeric factor ID (must stay first -
  # levels<- requires it); the SECOND column holds the actual category labels (here, the mukey
  # codes) as CHARACTER, not numeric - both confirmed empirically (as.factor()/cats()'s exact
  # column shape isn't documented precisely enough to assume). dplyr::left_join() (not base
  # merge(), which reorders the `by` column to the front and breaks the "ID must stay first" rule)
  # preserves that column order.
  cats_table <- terra::cats(mukey_raster)[[1]]
  mukey_col <- names(cats_table)[2]

  percentile_cols <- setdiff(names(percentile_by_mukey), "mukey")
  percentile_by_mukey[[mukey_col]] <- as.character(percentile_by_mukey$mukey)

  merged <- dplyr::left_join(cats_table, percentile_by_mukey[, c(mukey_col, percentile_cols)], by = mukey_col)
  levels(mukey_raster) <- merged

  value_stack <- terra::catalyze(mukey_raster)

  stats::setNames(
    lapply(percentile_cols, function(p) value_stack[[p]]),
    percentile_cols
  )
}

#' Compute Per-Mukey Percentile Rasters from Already-Simulated SSURGO Draws
#'
#' The quantile/rasterize half of [fetch_ssurgo_percentiles()], factored out so a single shared
#' [simulate_ssurgo_mapunit_draws()] result can be reused across multiple properties instead of
#' resimulating once per property - `simulate_cokey_generalized()` already simulates every
#' recognized property jointly in one pass per cokey, so a caller that needs several properties
#' from the same AOI/depth window (e.g. `run_stage1_fusion_group()`'s texture members) only needs
#' to run the (expensive) simulation once and call this per property afterward.
#'
#' @param mukey_raster A categorical (factor) mukey `terra::SpatRaster` from
#'   `fetch_ssurgo_mukey_raster()`.
#' @param draws A data frame from `simulate_ssurgo_mapunit_draws()`, for the same AOI/depth
#'   window `mukey_raster` covers.
#' @param property_id One of `property_to_sim_column()`'s recognized ids.
#' @param probs Percentile probabilities to compute (default `c(0.05, 0.25, 0.5, 0.75, 0.95)`).
#' @return `list(values = <named list of percentile-value SpatRasters>, probs = probs)`, or `NULL`
#'   if `property_id`'s simulated column isn't present in `draws`.
#' @keywords internal
percentiles_from_draws <- function(mukey_raster, draws, property_id,
                                    probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  sim_col <- property_to_sim_column(property_id)
  if (!sim_col %in% names(draws)) return(NULL)

  percentile_names <- sprintf("P%02d", round(probs * 100))
  by_mukey <- draws |>
    dplyr::group_by(mukey) |>
    dplyr::summarise(q = list(stats::quantile(.data[[sim_col]], probs = probs, na.rm = TRUE)), .groups = "drop")
  q_mat <- do.call(rbind, by_mukey$q)
  colnames(q_mat) <- percentile_names
  percentile_by_mukey <- cbind(mukey = by_mukey$mukey, as.data.frame(q_mat))

  values <- rasterize_mukey_percentiles(mukey_raster, percentile_by_mukey)

  list(values = values, probs = probs)
}

#' Fetch SSURGO Percentile-Value Rasters for an AOI
#'
#' The top-level SSURGO "prior" entry point for `R/raster-fusion.R`'s `fuse_property_adaptive()`:
#' rasterizes map units, Monte Carlo-simulates the requested property, and returns per-cell
#' percentile-value rasters in the `list(values=, probs=)` shape `fuse_property_adaptive()` expects.
#'
#' @param aoi_vect A `terra::SpatVector` AOI.
#' @param property_id One of `property_to_sim_column()`'s recognized ids.
#' @param top_depth,bottom_depth Numeric depth window bounds in cm.
#' @param probs Percentile probabilities to compute (default `c(0.05, 0.25, 0.5, 0.75, 0.95)`).
#' @param n_mc Number of triangular draws per component (default 1000).
#' @return `list(values = <named list of percentile-value SpatRasters>, probs = probs)`, or `NULL`
#'   if the mukey raster or the Monte Carlo draws are unavailable for this AOI.
#' @export
fetch_ssurgo_percentiles <- function(aoi_vect, property_id, top_depth, bottom_depth,
                                      probs = c(0.05, 0.25, 0.5, 0.75, 0.95), n_mc = 1000) {
  mukey_raster <- fetch_ssurgo_mukey_raster(aoi_vect)
  if (is.null(mukey_raster)) return(NULL)

  draws <- simulate_ssurgo_mapunit_draws(aoi_vect, top_depth, bottom_depth, n_mc)
  if (is.null(draws) || nrow(draws) == 0) return(NULL)

  percentiles_from_draws(mukey_raster, draws, property_id, probs)
}
