# ---------------------------------------------------------------------------
# raster-fusion-bridge.R : A.3 zonal_distribution_from_posterior(),
#                          A.4 remarginalize_ensemble_to_posterior(),
#                          A.5 remarginalized_awc()
# All offline: synthetic mukey rasters + synthetic run_stage1_fusion()-shaped posteriors.
# ---------------------------------------------------------------------------

.default_probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)

# a 4x4 factor mukey raster: top half mukey "100", bottom half mukey "200"
.mk_raster <- function() {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
                   vals = rep(c(100L, 100L, 200L, 200L), each = 4))
  names(r) <- "mukey"
  terra::as.factor(r)
}

# posterior with spatially-CONSTANT percentile layers at qnorm(probs, mean, sd)
.const_posterior <- function(template, mean = 10, sd = 2, probs = .default_probs) {
  vals <- stats::qnorm(probs, mean, sd)
  list(percentiles = stats::setNames(
    lapply(vals, function(v) terra::setValues(template, rep(v, terra::ncell(template)))),
    paste0("P", round(probs * 100))
  ))
}

# posterior whose P5/P50/P95 differ by mukey region (for the zonal test)
.regional_posterior <- function(mk) {
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  layer <- function(v100, v200) {
    terra::setValues(tmpl, ifelse(terra::values(mk)[, 1] == 100, v100, v200))
  }
  list(percentiles = list(
    P5  = layer(1, 4),
    P50 = layer(2, 5),
    P95 = layer(3, 6)
  ))
}

# synthetic ensemble: 2 mukeys, correlated wr_3b/wr_15b + rfv, 2 windows
.synth_ensemble <- function(mk, n = 40) {
  set.seed(1)
  one_window <- function(bump) {
    mkcode <- rep(c("100", "200"), each = n)
    wr3 <- c(stats::runif(n, 25, 40), stats::runif(n, 20, 35)) + bump
    wr15 <- wr3 * 0.45 + stats::rnorm(2 * n, 0, 0.5)   # strong rank dependence on wr3
    data.frame(
      mukey = mkcode,
      cokey = "1",
      simulation_number = rep(seq_len(n), 2),
      wr_3b = wr3, wr_15b = wr15, rfv = c(stats::runif(n, 0, 15), stats::runif(n, 5, 25)),
      top = 0, bottom = 5
    )
  }
  dbw <- list("0-5" = one_window(0), "5-15" = one_window(1.5))
  extract_mukey_joint_ensemble(
    aoi_vect = NULL, depth_windows = list(c(0, 5), c(5, 15)),
    draws_by_window = dbw, mukey_raster = mk,
    properties = c("wr_3b", "wr_15b", "rfv")
  )
}

# synthetic ensemble for the Saxton-Rawls AWC path: sand/silt/clay/db/soc/rfv, 2 windows
.synth_ensemble_sr <- function(mk, n = 40) {
  set.seed(3)
  one_window <- function(bump) {
    sand <- c(stats::runif(n, 35, 55), stats::runif(n, 20, 40)) + bump
    clay <- c(stats::runif(n, 12, 25), stats::runif(n, 20, 35))
    data.frame(
      mukey = rep(c("100", "200"), each = n), cokey = "1",
      simulation_number = rep(seq_len(n), 2),
      sand_total = sand, clay_total = clay, silt_total = pmax(0, 100 - sand - clay),
      db = c(stats::runif(n, 1.3, 1.5), stats::runif(n, 1.2, 1.45)),
      soc = c(stats::runif(n, 0.5, 1.5), stats::runif(n, 1.0, 2.5)),
      rfv = c(stats::runif(n, 0, 10), stats::runif(n, 5, 20)),
      top = 0, bottom = 5
    )
  }
  dbw <- list("0-5" = one_window(0), "5-15" = one_window(2))
  extract_mukey_joint_ensemble(
    aoi_vect = NULL, depth_windows = list(c(0, 5), c(5, 15)),
    draws_by_window = dbw, mukey_raster = mk,
    properties = c("sand_total", "silt_total", "clay_total", "db", "soc", "rfv")
  )
}

test_that("zonal_distribution_from_posterior() reduces each mukey to its own low/rep/high", {
  mk <- .mk_raster()
  post <- .regional_posterior(mk)
  z <- zonal_distribution_from_posterior(post, mk, probs = c(0.05, 0.5, 0.95))

  expect_setequal(names(z), c("100", "200"))
  expect_equal(unname(z[["100"]]), c(1, 2, 3))
  expect_equal(unname(z[["200"]]), c(4, 5, 6))
  expect_named(z[["100"]], c("low", "rep", "high"))
})

test_that("zonal_distribution_from_posterior() rejects a mean path / non-matching probs", {
  mk <- .mk_raster()
  post <- .regional_posterior(mk)
  expect_error(zonal_distribution_from_posterior(post, mk, probs = c(0.5, 0.05, 0.95)), "ascending")
  expect_error(zonal_distribution_from_posterior(post, mk, probs = c(0.05, 0.5, 0.90)), "exact posterior percentile")
})

test_that("remarginalize_ensemble_to_posterior() recovers the posterior marginal per pixel", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble(mk)

  pbpw <- list(wr_3b = list("0-5" = .const_posterior(tmpl, mean = 30, sd = 3)))
  out <- remarginalize_ensemble_to_posterior(ens, pbpw, n_out = 30, summarize = TRUE,
                                             probs = c(0.1, 0.5, 0.9))

  got <- out$percentiles$wr_3b$`0-5`
  # posterior is spatially constant -> every pixel's re-marginalized quantiles match the
  # piecewise-linear inverse-CDF of the posterior knots at the requested probs
  knot_p <- .default_probs
  knot_v <- stats::qnorm(knot_p, 30, 3)
  expect_equal(as.numeric(terra::global(got$P10, "mean", na.rm = TRUE)),
               stats::approx(knot_p, knot_v, xout = 0.1, rule = 2)$y, tolerance = 0.5)
  expect_equal(as.numeric(terra::global(got$P50, "mean", na.rm = TRUE)),
               stats::approx(knot_p, knot_v, xout = 0.5, rule = 2)$y, tolerance = 0.3)
  # spatially constant (top vs bottom mukey identical, since posterior is constant)
  expect_lt(as.numeric(terra::global(got$P50, "sd", na.rm = TRUE)), 1e-6)
})

test_that("remarginalize_ensemble_to_posterior() preserves the cross-property rank copula", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble(mk)

  pbpw <- list(
    wr_3b  = list("0-5" = .const_posterior(tmpl, mean = 30, sd = 4)),
    wr_15b = list("0-5" = .const_posterior(tmpl, mean = 13, sd = 2))
  )
  out <- remarginalize_ensemble_to_posterior(ens, pbpw, n_out = 30, summarize = FALSE)
  n_kept <- out$n_kept

  # source Spearman for mukey 100 (first n_kept rows)
  src <- ens$by_mukey[["100"]]$windows[["0-5"]][seq_len(n_kept), ]
  src_rho <- suppressWarnings(stats::cor(src[, "wr_3b"], src[, "wr_15b"], method = "spearman"))

  # transformed Spearman at a mukey-100 pixel (cell 1)
  a <- as.numeric(terra::values(out$ensemble$wr_3b$`0-5`)[1, ])
  b <- as.numeric(terra::values(out$ensemble$wr_15b$`0-5`)[1, ])
  trans_rho <- suppressWarnings(stats::cor(a, b, method = "spearman"))

  expect_equal(trans_rho, src_rho, tolerance = 1e-6)
})

test_that("remarginalize_ensemble_to_posterior() gives NA (not fabricated values) for a mukey with no data for a property", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble(mk)
  # wipe mukey 200's wr_3b for the 0-5 window (e.g. every cokey's sim failed there)
  ens$by_mukey[["200"]]$windows[["0-5"]][, "wr_3b"] <- NA_real_

  pbpw <- list(
    wr_3b  = list("0-5" = .const_posterior(tmpl, mean = 30, sd = 4)),
    wr_15b = list("0-5" = .const_posterior(tmpl, mean = 13, sd = 2))
  )
  out <- remarginalize_ensemble_to_posterior(ens, pbpw, n_out = 30, summarize = FALSE)

  v <- terra::values(out$ensemble$wr_3b$`0-5`)
  mk2 <- mk; levels(mk2) <- NULL
  mkv <- terra::values(mk2)[, 1]                       # bottom-half cells are mukey 200
  # base rank()'s default na.last = TRUE would have ranked the NAs 1..n and produced finite values;
  # na.last = "keep" keeps them NA -> those cells are NA.
  expect_true(all(is.na(v[mkv == 200, ])))
  expect_false(any(is.na(v[mkv == 100, ])))            # mukey 100 untouched
})

.mk_post <- function(tmpl, m, s) {
  list("0-5" = .const_posterior(tmpl, m, s), "5-15" = .const_posterior(tmpl, m, s))
}

test_that("saxton_rawls_raster() matches calculate_saxton_rawls_single() cell-for-cell", {
  tmpl <- terra::rast(nrows = 2, ncols = 2, vals = 1)
  cases <- list(
    c(45, 18, 37, 1.4, 5, 1.0), c(25, 32, 43, 1.3, 15, 2.5), c(70, 8, 22, 1.6, 0, 0.5),
    c(50, 25, 40, 1.4, 5, 1.0),   # sand+silt+clay = 115 -> renormalization branch
    c(30, 20, 35, 1.35, 8, 2.0)   # = 85 -> renormalization branch
  )
  for (cc in cases) {
    single <- calculate_saxton_rawls_single(cc[1], cc[2], cc[3], cc[4], cc[5], cc[6])
    r <- saxton_rawls_raster(terra::setValues(tmpl, cc[1]), terra::setValues(tmpl, cc[2]),
                             terra::setValues(tmpl, cc[3]), terra::setValues(tmpl, cc[4]),
                             terra::setValues(tmpl, cc[5]), terra::setValues(tmpl, cc[6]))
    expect_equal(terra::values(r$fc)[1], single$field_capacity, tolerance = 0.05,
                 info = paste(cc, collapse = ","))
    expect_equal(terra::values(r$wp)[1], single$wilting_point, tolerance = 0.05,
                 info = paste(cc, collapse = ","))
    expect_lt(terra::values(r$wp)[1], terra::values(r$fc)[1])
  }
})

test_that("remarginalized_awc(method = 'saxton_rawls') returns clamped per-pixel AWC (default path)", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)

  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )
  res <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.05, 0.5, 0.95))  # method defaults

  expect_named(res$awc_cm, c("P5", "P50", "P95"))
  expect_setequal(res$windows_used, c("0-5", "5-15"))
  for (lyr in res$awc_cm) expect_true(all(terra::values(lyr) >= 0, na.rm = TRUE))
  # 0-15 cm of a loam: AWC ~ 0.12-0.18 vol-frac diff over 15 cm ~ 1.5-3 cm - loose sanity band
  med <- as.numeric(terra::global(res$awc_cm$P50, "mean", na.rm = TRUE))
  expect_gt(med, 0.5); expect_lt(med, 6)
})

test_that("remarginalized_awc(method = 'direct') uses supplied wr_3b/wr_15b posteriors", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble(mk)
  pbpw <- list(wr_3b = .mk_post(tmpl, 30, 3), wr_15b = .mk_post(tmpl, 13, 2), rfv = .mk_post(tmpl, 8, 3))

  res <- remarginalized_awc(ens, pbpw, method = "direct", n_out = 25, probs = c(0.05, 0.5, 0.95))
  med <- as.numeric(terra::global(res$awc_cm$P50, "mean", na.rm = TRUE))
  expect_gt(med, 1); expect_lt(med, 4)   # (30-13)/100 * 15 * (1-.08) ~ 2.3 cm
})

test_that("remarginalized_awc() errors when required properties are missing", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens_wr <- .synth_ensemble(mk)                    # has wr_3b/wr_15b, not texture
  expect_error(remarginalized_awc(ens_wr, list(wr_3b = .mk_post(tmpl, 30, 3))),
               "sand_total/silt_total/clay_total/db")               # saxton_rawls default
  expect_error(remarginalized_awc(ens_wr, list(wr_3b = .mk_post(tmpl, 30, 3)), method = "direct"),
               "wr_3b/wr_15b")
})

test_that("remarginalized_awc(tile_rows=) is numerically identical to the whole-grid run", {
  r <- terra::rast(nrows = 12, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 12,
                   vals = rep(c(100L, 200L), each = 36))
  names(r) <- "mukey"; mk <- terra::as.factor(r)
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)

  pbpw <- list(sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
               clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08))

  whole <- remarginalized_awc(ens, pbpw, n_out = 20)
  tiled <- remarginalized_awc(ens, pbpw, n_out = 20, tile_rows = 4)

  expect_equal(tiled$n_tiles, 3)
  expect_equal(whole$n_tiles, 1)
  for (nm in names(whole$awc_cm)) {
    expect_equal(terra::values(tiled$awc_cm[[nm]]), terra::values(whole$awc_cm[[nm]]), tolerance = 1e-9)
  }
})

# ---------------------------------------------------------------------------
# S2 - remarginalized_awc(restriction_depth=): bedrock/restriction-depth AWC truncation
# (MULTI_PROPERTY_FUSION_PLAN.md task S2). Windows are "0-5" (top=0,bottom=5) and "5-15"
# (top=5,bottom=15) throughout, from .synth_ensemble_sr()'s fixed depth_windows.
# ---------------------------------------------------------------------------

test_that("remarginalized_awc(restriction_depth = NULL) is bit-identical to omitting the argument", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)
  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )

  set.seed(42); a <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.05, 0.5, 0.95))
  set.seed(42); b <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.05, 0.5, 0.95),
                                       restriction_depth = NULL)
  for (nm in names(a$awc_cm)) {
    expect_equal(terra::values(b$awc_cm[[nm]]), terra::values(a$awc_cm[[nm]]))
  }
})

test_that("remarginalized_awc(restriction_depth=) reduces AWC for a restriction straddling a window", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)
  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )
  # 10 cm sits inside the "5-15" window (top=5,bottom=15): that window's effective thickness is
  # 10-5=5 instead of its nominal 15-5=10 - straddling, partial credit, not full or zero.
  rd <- terra::setValues(tmpl, rep(10, terra::ncell(tmpl)))

  set.seed(42); unrestricted <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5))
  set.seed(42); restricted   <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5),
                                                   restriction_depth = rd)

  u <- as.numeric(terra::global(unrestricted$awc_cm$P50, "mean", na.rm = TRUE))
  r <- as.numeric(terra::global(restricted$awc_cm$P50, "mean", na.rm = TRUE))
  expect_lt(r, u)
  expect_gt(r, 0)  # partial credit, not zeroed out
})

test_that("remarginalized_awc(restriction_depth=) zeroes (not NAs) a window entirely below the restriction", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)
  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )
  # 3 cm is above the "5-15" window's top (5) entirely -> that window contributes exactly 0, not
  # NA (which would incorrectly blank the whole pixel's AWC, including the valid "0-5" window).
  rd <- terra::setValues(tmpl, rep(3, terra::ncell(tmpl)))

  set.seed(1); res <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5), restriction_depth = rd)
  expect_false(anyNA(terra::values(res$awc_cm$P50)))
  expect_true(all(terra::values(res$awc_cm$P50) >= 0, na.rm = TRUE))

  # rd=3 truncates BOTH windows (even "0-5" itself: pmax(0, pmin(5,3) - 0) = 3 cm, not its
  # nominal 5) - so the restricted AWC must be strictly less than the fully unrestricted run.
  set.seed(1); unrestricted <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5))
  expect_lt(as.numeric(terra::global(res$awc_cm$P50, "mean", na.rm = TRUE)),
            as.numeric(terra::global(unrestricted$awc_cm$P50, "mean", na.rm = TRUE)))
})

test_that("remarginalized_awc(restriction_depth = Inf) matches the unrestricted run", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)
  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )
  rd <- terra::setValues(tmpl, rep(Inf, terra::ncell(tmpl)))

  set.seed(7); unrestricted <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5))
  set.seed(7); restricted   <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5),
                                                  restriction_depth = rd)
  expect_equal(terra::values(restricted$awc_cm$P50), terra::values(unrestricted$awc_cm$P50),
              tolerance = 1e-9)
})

test_that("remarginalized_awc(restriction_depth=) treats NA cells as unrestricted (Inf), not propagated NA", {
  mk <- .mk_raster()
  tmpl <- terra::rast(mk); names(tmpl) <- "v"
  ens <- .synth_ensemble_sr(mk)
  pbpw <- list(
    sand_total = .mk_post(tmpl, 42, 4), silt_total = .mk_post(tmpl, 38, 4),
    clay_total = .mk_post(tmpl, 20, 3), db = .mk_post(tmpl, 1.4, 0.08),
    soc = .mk_post(tmpl, 1.2, 0.3), rfv = .mk_post(tmpl, 8, 3)
  )
  rd <- terra::setValues(tmpl, rep(NA_real_, terra::ncell(tmpl)))

  set.seed(9); unrestricted <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5))
  set.seed(9); na_restricted <- remarginalized_awc(ens, pbpw, n_out = 25, probs = c(0.5),
                                                   restriction_depth = rd)
  expect_false(anyNA(terra::values(na_restricted$awc_cm$P50)))
  expect_equal(terra::values(na_restricted$awc_cm$P50), terra::values(unrestricted$awc_cm$P50),
              tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# A.7 - live end-to-end (real AOI -> ensemble -> stage-1 fusion -> re-marginalize -> AWC)
# ---------------------------------------------------------------------------

test_that("end-to-end: extract_mukey_joint_ensemble -> run_stage1_fusion -> remarginalized_awc (live)", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  Sys.unsetenv("PROJ_LIB")   # documented terra/PROJ env quirk - see HANDOFF_NOTES.md / .onLoad()

  aoi <- terra::vect(
    "POLYGON((-121.652 36.610, -121.650 36.610, -121.650 36.612, -121.652 36.612, -121.652 36.610))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
                    "AOI has no valid CRS (documented PROJ_LIB issue) - not a soilSIM defect.")
  aoi <- terra::project(aoi, "epsg:5070")

  windows <- list(c(0, 5), c(5, 15))
  ens <- tryCatch(extract_mukey_joint_ensemble(aoi, windows), error = function(e) NULL)
  testthat::skip_if(is.null(ens), "live SDA/SSURGO simulation unavailable in this run.")

  # SOLUS100 has no water-retention variable, so the Saxton-Rawls path fuses texture + Db instead.
  cfg <- function(id, sv) list(id = id, solus_variable = sv, dist = "auto")
  props <- list(sand_total = "sandtotal", silt_total = "silttotal",
                clay_total = "claytotal", db = "dbovendry")
  post <- stats::setNames(vector("list", length(props)), names(props))
  for (nm in names(props)) {
    per_w <- list()
    for (w in windows) {
      r <- tryCatch(run_stage1_fusion(aoi, cfg(nm, props[[nm]]), w[1], w[2]), error = function(e) NULL)
      if (!is.null(r)) per_w[[paste0(w[1], "-", w[2])]] <- list(percentiles = r$posterior$percentiles)
    }
    post[[nm]] <- per_w
  }
  testthat::skip_if(any(vapply(post, length, integer(1)) == 0),
                    "live SOLUS fusion unavailable in this run.")

  res <- remarginalized_awc(ens, post, n_out = 100)
  expect_named(res$awc_cm, c("P5", "P25", "P50", "P75", "P95"))
  expect_s4_class(res$awc_cm$P50, "SpatRaster")
  v <- terra::values(res$awc_cm$P50)
  expect_true(all(v >= 0, na.rm = TRUE))
  expect_true(any(is.finite(v)))
  # P5 <= P50 <= P95 everywhere
  expect_true(all(terra::values(res$awc_cm$P5)  <= terra::values(res$awc_cm$P50) + 1e-6, na.rm = TRUE))
  expect_true(all(terra::values(res$awc_cm$P50) <= terra::values(res$awc_cm$P95) + 1e-6, na.rm = TRUE))
})
