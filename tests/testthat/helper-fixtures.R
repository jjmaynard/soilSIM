# Shared synthetic-data builders for soilSIM tests.

#' Build one synthetic horizon row with l/r/h triplets for the requested properties
#'
#' @param properties Character vector of property names (no `_l`/`_r`/`_h` suffix).
#' @param ... Named arguments of the form `<property>_l`, `<property>_r`,
#'   `<property>_h` overriding the defaults below.
make_horizon_row <- function(properties = c("dbovendry", "sandtotal", "claytotal", "silttotal"),
                              hzdept_r = 0, hzdepb_r = 20, genhz = "A", cokey = "1", mukey = "1", ...) {
  overrides <- list(...)
  row <- data.frame(cokey = cokey, mukey = mukey, hzdept_r = hzdept_r, hzdepb_r = hzdepb_r,
                     genhz = genhz, hzname = "A", stringsAsFactors = FALSE)

  defaults <- list(
    dbovendry = c(l = 1.2, r = 1.4, h = 1.6),
    sandtotal = c(l = 30, r = 40, h = 50),
    claytotal = c(l = 15, r = 20, h = 25),
    silttotal = c(l = 30, r = 40, h = 50),
    ph1to1h2o = c(l = 5.5, r = 6.0, h = 6.5),
    cec7 = c(l = 8, r = 12, h = 18),
    om = c(l = 1, r = 2, h = 4),
    rfv = c(l = 5, r = 10, h = 20)
  )

  for (p in properties) {
    d <- defaults[[p]] %||% c(l = 1, r = 2, h = 3)
    row[[paste0(p, "_l")]] <- overrides[[paste0(p, "_l")]] %||% d[["l"]]
    row[[paste0(p, "_r")]] <- overrides[[paste0(p, "_r")]] %||% d[["r"]]
    row[[paste0(p, "_h")]] <- overrides[[paste0(p, "_h")]] %||% d[["h"]]
  }
  row
}

#' Build a synthetic multi-horizon soil_data data frame
#'
#' @param n_horizons Number of horizon rows to generate.
#' @param properties Character vector of property names.
make_soil_data <- function(n_horizons = 20, properties = c("dbovendry", "sandtotal", "claytotal", "silttotal")) {
  set.seed(123)
  rows <- lapply(seq_len(n_horizons), function(i) {
    genhz <- sample(c("A", "B", "C"), 1)
    jitter <- stats::rnorm(1, 0, 1)
    make_horizon_row(
      properties = properties,
      hzdept_r = (i - 1) * 10, hzdepb_r = i * 10,
      genhz = genhz, cokey = as.character(((i - 1) %/% 3) + 1), mukey = "1",
      dbovendry_r = 1.4 + jitter * 0.05
    )
  })
  do.call(rbind, rows)
}

#' Build a synthetic soil_data data frame with ONLY texture (sand/clay/silt) properties
make_texture_only_soil_data <- function(n_horizons = 20) {
  make_soil_data(n_horizons, properties = c("sandtotal", "claytotal", "silttotal"))
}

#' Build a synthetic soil_data data frame with real, varying hzname values and NO genhz column
#'
#' `make_horizon_row()`/`make_soil_data()` always set a `genhz` column (and a
#' fixed `hzname = "A"`), so neither can exercise
#' `estimate_property_correlations()`'s `classify_genhz(hzname)`
#' auto-derivation path (only active when
#' `correlation_fallback = "kssl_global"` is requested AND `genhz` is
#' absent). This builder cycles through real hzname strings that
#' `classify_genhz()` maps to A/B/C and never sets a `genhz` column.
#' `n_horizons` defaults to keep each hzname group under the default
#' `min_group_n = 5`, so empirical per-group estimation fails and the KSSL
#' fallback path is actually exercised.
#'
#' @param n_horizons Number of horizon rows to generate.
#' @param properties Character vector of property names (should include at
#'   least 2 with a KSSL mapping - see `.kssl_property_name_map` - to
#'   produce a non-identity fallback matrix; default covers 3).
make_soil_data_hzname_only <- function(n_horizons = 12,
                                        properties = c("dbovendry", "ph1to1h2o", "cec7")) {
  hznames <- c("Ap", "Bt1", "Bt2", "Bw", "C1", "C2")
  # Base _r value per property (matching make_horizon_row()'s own defaults),
  # so each can be jittered below - without this, only dbovendry varies and
  # every other property is a constant column, which Hmisc::rcorr() can't
  # compute a correlation for (undefined variance).
  base_r <- c(dbovendry = 1.4, sandtotal = 40, claytotal = 20, silttotal = 40,
              ph1to1h2o = 6.0, cec7 = 12, om = 2, rfv = 10)
  set.seed(456)
  rows <- lapply(seq_len(n_horizons), function(i) {
    overrides <- lapply(properties, function(p) {
      jitter <- stats::rnorm(1, 0, 1)
      (base_r[[p]] %||% 2) + jitter * (base_r[[p]] %||% 2) * 0.05
    })
    names(overrides) <- paste0(properties, "_r")

    row <- do.call(make_horizon_row, c(
      list(properties = properties, hzdept_r = (i - 1) * 10, hzdepb_r = i * 10,
           cokey = as.character(((i - 1) %/% 3) + 1), mukey = "1"),
      overrides
    ))
    row$hzname <- hznames[((i - 1) %% length(hznames)) + 1]
    row$genhz <- NULL
    row
  })
  do.call(rbind, rows)
}

#' Build a synthetic multi-horizon data frame with hzdept/hzdepb l/r/h depth
#' triplets, for `R/depth-simulation.R`'s plain-data-frame functions
#' (`simulate_soil_profile_top_down()`/`_bottom_up()`/`_thickness()`,
#' `infill_missing_depth_variability()`).
#'
#' @param cokey,mukey,compname Single values applied to every horizon (one profile).
make_depth_horizon_data <- function(cokey = "1", mukey = "1", compname = "testseries") {
  data.frame(
    cokey = cokey, mukey = mukey, compname = compname,
    hzname = c("A", "Bt", "Cr"),
    hzdept_l = c(0, 15, 45), hzdept_r = c(0, 20, 50), hzdept_h = c(0, 25, 55),
    hzdepb_l = c(15, 45, 95), hzdepb_r = c(20, 50, 100), hzdepb_h = c(25, 55, 105),
    stringsAsFactors = FALSE
  )
}

#' Build a synthetic `aqp::SoilProfileCollection` for
#' `R/depth-simulation.R`'s SoilProfileCollection-level functions
#' (`simulate_and_perturb_soil_profiles()` and friends).
#'
#' Adds a constant `sim_comppct` column since `simulate_and_perturb_soil_profiles()`
#' derives its simulation count from it (see that function's `@section Known
#' limitation:` - normally populated by a composition-simulation step this
#' package doesn't yet implement).
#'
#' @param sim_comppct Number of perturbed replicates to simulate per profile.
make_soil_profile_collection <- function(sim_comppct = 5) {
  df <- make_depth_horizon_data()
  df$sim_comppct <- sim_comppct
  df$id <- df$compname
  aqp::depths(df) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(df) <- "hzname"
  df
}

#' Build a synthetic data frame shaped exactly like
#' `soilDB::ROSETTA(..., include.sd = TRUE)`'s output (plus a caller-added
#' `layerID`), for offline `simulate_vg_aws()` testing (`R/aws-simulation.R`).
#' `alpha`/`npar` are in log10 space, matching `simulate_vg_aws()`'s
#' `10^(...)` back-transform.
make_rosetta_shaped_data <- function() {
  data.frame(
    compname = c("testseries", "testseries"),
    hzname = c("A", "Bt"),
    layerID = c("testseries_0", "testseries_20"),
    alpha = c(-1.5, -1.8), sd_alpha = c(0.1, 0.1),
    npar = c(0.15, 0.12), sd_npar = c(0.02, 0.02),
    theta_r = c(0.06, 0.08), sd_theta_r = c(0.01, 0.01),
    theta_s = c(0.42, 0.40), sd_theta_s = c(0.02, 0.02),
    stringsAsFactors = FALSE
  )
}

#' Build a synthetic multi-horizon data frame with ROSETTA's expected input
#' variable names (`sand_total`/`silt_total`/`clay_total`/
#' `bulk_density_third_bar`/`water_retention_third_bar`/
#' `water_retention_15_bar`) plus `compname`/`hzdept_r`/`hzdepb_r`/`cokey`,
#' for the live `calculate_aws_df()` test (`R/aws-simulation.R`).
make_aws_texture_data <- function() {
  data.frame(
    cokey = c("1", "1"), compname = c("testseries", "testseries"),
    hzdept_r = c(0, 20), hzdepb_r = c(20, 50),
    sand_total = c(40, 35), silt_total = c(40, 40), clay_total = c(20, 25),
    bulk_density_third_bar = c(1.35, 1.42),
    water_retention_third_bar = c(0.28, 0.30),
    water_retention_15_bar = c(0.12, 0.14),
    stringsAsFactors = FALSE
  )
}

#' Build a small synthetic percentile-value `terra::SpatRaster` stack, for
#' `R/distribution-fitting-raster.R`/`R/raster-fusion.R` tests (the first
#' `SpatRaster`-native code in soilSIM - no prior fixture pattern to reuse).
#'
#' Every layer is spatially constant (same value in every cell) by default,
#' so results can be cross-checked cell-by-cell against the equivalent
#' scalar functions in `R/distributions.R`/`R/bayesian-updating.R`.
#'
#' @param values_by_prob Named numeric vector, e.g. `c(P5 = 2, P50 = 10, P95 = 18)`.
#' @param nrow,ncol Raster dimensions (small - these are cheap correctness
#'   fixtures, not performance tests).
#' @return A named list of single-layer SpatRasters, one per `values_by_prob` entry.
make_percentile_rasters <- function(values_by_prob, nrow = 2, ncol = 2) {
  ncells <- nrow * ncol
  stats::setNames(
    lapply(values_by_prob, function(v) terra::rast(nrows = nrow, ncols = ncol, vals = rep(v, ncells))),
    names(values_by_prob)
  )
}

#' Build synthetic component-level data for `sim_component_comp()` testing
#' (`R/property-simulation.R`) - `mukey`/`cokey`/`compname`/`comppct_l/r/h`,
#' one row per component.
make_component_data <- function() {
  data.frame(
    mukey = c("1", "1"), cokey = c("1", "2"), compname = c("compA", "compB"),
    comppct_l = c(45, 25), comppct_r = c(55, 35), comppct_h = c(65, 45),
    stringsAsFactors = FALSE
  )
}

#' Build synthetic data for `recover_missing_horizon_components()`/
#' `synthesize_component_horizons_from_siblings()` testing (`R/ssurgo-acquisition.R`).
#'
#' `ssurgo_data`: two sibling cokeys ("10","11") sharing `compname = "compA"`, each with an A and
#' a Bt horizon (deliberately different exact depths/values so averaging is verifiable), plus one
#' unrelated cokey ("12", `compB`) that should never be touched by recovery logic aimed at compA.
#'
#' `all_components`: the full component-only table (no chorizon join) - includes all three cokeys
#' above (redundant with `ssurgo_data`, as a real `fetch_ssurgo_all_components_working()` result
#' would be) PLUS a 4th component, cokey "20" under a DIFFERENT mukey ("2"), also named "compA",
#' with a real `comppct` but zero horizon rows in `ssurgo_data` - the gap this feature recovers.
#' Deliberately placed under a different mukey than its siblings to verify sibling search is
#' AOI-wide (matching `compname` anywhere in the AOI's `ssurgo_data`), not restricted to one mukey.
make_component_recovery_fixture <- function() {
  sibling1 <- data.frame(
    mukey = "1", cokey = "10", compname = "compA",
    comppct_l = 40, comppct_r = 50, comppct_h = 60,
    hzname = c("A", "Bt"), hzdept_r = c(0, 20), hzdepb_r = c(20, 60),
    sandtotal_r = c(40, 30), claytotal_r = c(15, 25),
    unsuitable_horizon = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  sibling2 <- data.frame(
    mukey = "1", cokey = "11", compname = "compA",
    comppct_l = 40, comppct_r = 50, comppct_h = 60,
    hzname = c("A", "Bt"), hzdept_r = c(0, 18), hzdepb_r = c(18, 55),
    sandtotal_r = c(44, 34), claytotal_r = c(19, 29),
    unsuitable_horizon = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  compB <- data.frame(
    mukey = "1", cokey = "12", compname = "compB",
    comppct_l = 10, comppct_r = 15, comppct_h = 20,
    hzname = "A", hzdept_r = 0, hzdepb_r = 25,
    sandtotal_r = 55, claytotal_r = 10,
    unsuitable_horizon = FALSE,
    stringsAsFactors = FALSE
  )
  ssurgo_data <- dplyr::bind_rows(sibling1, sibling2, compB)

  all_components <- data.frame(
    mukey = c("1", "1", "1", "2"),
    cokey = c("10", "11", "12", "20"),
    compname = c("compA", "compA", "compB", "compA"),
    comppct_l = c(40, 40, 10, 3), comppct_r = c(50, 50, 15, 5), comppct_h = c(60, 60, 20, 7),
    stringsAsFactors = FALSE
  )

  list(ssurgo_data = ssurgo_data, all_components = all_components)
}

#' Build a genhz-keyed list of small, positive-definite correlation matrices
#' for `simulate_cokey_generalized()` testing (`R/property-simulation.R`),
#' covering `db`/`ph` plus the `ilr1`/`ilr2` texture pseudo-properties (present
#' so the same fixture works whether or not the texture path is exercised -
#' `simulate_cokey_generalized()` only subsets the columns it actually needs).
make_property_correlation_matrices <- function(genhz = c("A", "B")) {
  vars <- c("db", "ph", "ilr1", "ilr2")
  m <- diag(4)
  dimnames(m) <- list(vars, vars)
  m["db", "ph"] <- m["ph", "db"] <- 0.3
  m["ilr1", "ilr2"] <- m["ilr2", "ilr1"] <- 0.2
  stats::setNames(lapply(genhz, function(g) m), genhz)
}

#' Build a genhz-keyed list of 3x3 sand/silt/clay correlation matrices for
#' `simulate_cokey_generalized()`'s texture step (`R/property-simulation.R`).
make_texture_correlation_matrices <- function(genhz = c("A", "B")) {
  m <- matrix(c(1, -0.4, -0.4, -0.4, 1, -0.3, -0.4, -0.3, 1), nrow = 3)
  stats::setNames(lapply(genhz, function(g) m), genhz)
}

#' Build a synthetic single-cokey `sim_cokey` data frame for
#' `simulate_cokey_generalized()` testing, with `genhz`, `sim_comppct`, and
#' `_l/_r/_h` triplets for `dbovendry`/`ph1to1h2o` (and, when `texture =
#' TRUE`, `sandtotal`/`silttotal`/`claytotal`).
make_sim_cokey_data <- function(texture = FALSE, sim_comppct = 20) {
  row <- make_horizon_row(
    properties = if (texture) c("dbovendry", "ph1to1h2o", "sandtotal", "silttotal", "claytotal") else c("dbovendry", "ph1to1h2o"),
    genhz = "A", cokey = "1", mukey = "1"
  )
  row$compname <- "testseries"
  row$sim_comppct <- sim_comppct
  row
}
