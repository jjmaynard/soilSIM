#' @title Disk Cache for Raster Fusion Fetch Results
#'
#' @description `build_cache_key()`/`cache_get()`/`cache_set()`/`CACHE_TTL_SECONDS` are called by
#'   `run_stage1_fusion()`/`run_stage1_fusion_group()` (`R/raster-fusion.R`) and
#'   `simulate_ssurgo_mapunit_draws()` (`R/ssurgo-simulation.R`) but were genuinely undefined
#'   anywhere in the source bundle they were ported from (`code_ref/reanalysis-platform/`) -
#'   confirmed by a repo-wide grep; `HANDOFF_NOTES.md` lists them as externals "assumed to already
#'   exist" in that bundle's *target* project, which was never soilSIM.
#'
#'   Rather than inventing a new scheme or adding a dependency (`memoise`/`R.cache`), this adapts
#'   the already-working, already-tested disk-RDS cache pattern from
#'   `R/ssurgo-acquisition.R`'s `generate_ssurgo_cache_key()`/`check_ssurgo_cache()`/
#'   `cache_ssurgo_data()` (`digest`-based keying, age/TTL invalidation via file mtime) -
#'   generalized here to the AOI/property-id/depth/kind shape the raster fusion code expects.
#'
#'   Cache files live under `tools::R_user_dir("soilSIM", "cache")` (the standard R >= 4.0
#'   per-package cache location) rather than a caller-supplied directory, since
#'   `run_stage1_fusion()`'s own calls (`cache_get(key, CACHE_TTL_SECONDS)`,
#'   `cache_set(key, "ssurgo", data)`) don't thread a directory argument through at all.
#' @name raster_cache
NULL

#' @rdname raster_cache
#' @export
CACHE_TTL_SECONDS <- 30 * 24 * 3600

#' @keywords internal
raster_cache_dir <- function() {
  tools::R_user_dir("soilSIM", "cache")
}

#' Build a cache key for one AOI/property/depth-window/kind combination
#'
#' @param aoi_vect A `terra::SpatVector` (single-feature AOI).
#' @param id A property id or composition-group name (e.g. `"ph"`, `"texture"`).
#' @param top_depth,bottom_depth Numeric depth bounds in cm.
#' @param kind A short string distinguishing what's cached for this key (e.g. `"ssurgo"`,
#'   `"solus"`, `"texture_group"`, `"posterior"`).
#' @return A character string, safe for use as a filename.
#' @export
build_cache_key <- function(aoi_vect, id, top_depth, bottom_depth, kind) {
  aoi_wkt <- terra::geom(aoi_vect, wkt = TRUE)[1]
  aoi_hash <- digest::digest(aoi_wkt)

  key <- paste0("raster_", substr(aoi_hash, 1, 12), "_", id, "_",
                top_depth, "_", bottom_depth, "_", kind)
  gsub("[^A-Za-z0-9_-]", "_", key)
}

#' Retrieve a cached value if present and not older than `ttl_seconds`
#'
#' @param key A cache key from `build_cache_key()`.
#' @param ttl_seconds Maximum cache age in seconds before a hit is treated as stale.
#' @return The cached value, or `NULL` on a cache miss/stale entry/read failure.
#' @export
cache_get <- function(key, ttl_seconds = CACHE_TTL_SECONDS) {
  cache_file <- file.path(raster_cache_dir(), paste0(key, ".rds"))

  if (!file.exists(cache_file)) {
    return(NULL)
  }

  age_seconds <- as.numeric(difftime(Sys.time(), file.mtime(cache_file), units = "secs"))
  if (age_seconds > ttl_seconds) {
    return(NULL)
  }

  tryCatch(readRDS(cache_file), error = function(e) NULL)
}

#' Store a value in the cache under `key`
#'
#' @param key A cache key from `build_cache_key()`.
#' @param kind Unused beyond documenting intent at call sites (matches the original bundle's
#'   3-argument `cache_set(key, kind, value)` calling convention) - the value is stored keyed
#'   only by `key`, since `build_cache_key()` already encodes `kind`.
#' @param value The R object to cache (anything `saveRDS()` can serialize, including
#'   `terra::SpatRaster` objects via their `wrap()`-compatible representation - see Known
#'   limitation below).
#' @return `TRUE` on success, `FALSE` on a write failure (never errors).
#'
#' @section Known limitation:
#' `terra::SpatRaster`/`SpatVector` objects hold an external pointer to in-memory/on-disk GDAL
#' state that does not survive a plain `saveRDS()`/`readRDS()` round-trip as-is. Callers caching
#' raster values must `terra::wrap()` them first and `terra::unwrap()` on read
#' (`cache_get()`/`cache_set()` do not do this automatically, since not every cached value is a
#' raster - e.g. Monte Carlo draw data frames are not).
#' @export
cache_set <- function(key, kind, value) {
  cache_dir <- raster_cache_dir()
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  cache_file <- file.path(cache_dir, paste0(key, ".rds"))

  tryCatch({
    saveRDS(value, cache_file)
    TRUE
  }, error = function(e) FALSE)
}

#' Wrap/unwrap the `SpatRaster`s inside a `list(values = <list of SpatRaster>, probs = )`
#' percentile structure for caching - see `cache_set()`'s `@section Known limitation:`.
#' @param x A `list(values=, probs=)` structure, as returned by `fetch_ssurgo_percentiles()`/
#'   `fetch_solus_percentiles()`.
#' @keywords internal
wrap_percentile_list <- function(x) {
  x$values <- lapply(x$values, terra::wrap)
  x
}
#' @rdname wrap_percentile_list
#' @keywords internal
unwrap_percentile_list <- function(x) {
  x$values <- lapply(x$values, terra::unwrap)
  x
}

#' Wrap/unwrap every `SpatRaster` found anywhere inside an arbitrarily-nested list
#'
#' A generic counterpart to `wrap_percentile_list()`/`unwrap_percentile_list()`, for result
#' structures whose shape isn't the fixed `list(values=, probs=)` percentile form - e.g.
#' `run_stage1_fusion()`'s return value, which nests `SpatRaster`s at `prior$values`,
#' `likelihood$values`, and inside `posterior` (whose own shape varies by `dist`). Recurses into
#' every list element, replacing each `SpatRaster` with `terra::wrap()`'s `PackedSpatRaster`
#' representation (or the reverse via `terra::unwrap()`) and leaving everything else untouched -
#' directly addresses `cache_set()`'s own `@section Known limitation:` for callers caching a
#' full fusion result rather than a raw percentile fetch.
#'
#' @param x A list, arbitrarily nested, that may contain `SpatRaster` objects at any depth.
#' @return `x`, with every `SpatRaster` replaced by its wrapped/unwrapped equivalent.
#' @export
wrap_nested_rasters <- function(x) {
  if (inherits(x, "SpatRaster")) {
    terra::wrap(x)
  } else if (is.list(x)) {
    lapply(x, wrap_nested_rasters)
  } else {
    x
  }
}
#' @rdname wrap_nested_rasters
#' @export
unwrap_nested_rasters <- function(x) {
  if (inherits(x, "PackedSpatRaster")) {
    terra::unwrap(x)
  } else if (is.list(x)) {
    lapply(x, unwrap_nested_rasters)
  } else {
    x
  }
}
