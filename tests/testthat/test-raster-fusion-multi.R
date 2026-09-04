# Offline tests for run_stage1_fusion_multi() (MULTI_PROPERTY_FUSION_PLAN.md Change A).
# Every fetch/simulate is mocked - the point is orchestration (one simulation, per-leaf fusion,
# cache seeding, partial-failure handling), not the fusion math (covered in test-raster-fusion.R).

probs5 <- c(0.05, 0.25, 0.5, 0.75, 0.95)

.mp_pct <- function(v, n = 1) {
  stats::setNames(
    lapply(probs5, function(p) terra::rast(nrows = 1, ncols = n, vals = rep(stats::qnorm(p, v, 2), n))),
    paste0("P", round(probs5 * 100))
  )
}
.mp_prior <- function(v = 20) list(values = .mp_pct(v), probs = probs5)
.mp_solus <- function(v = 22) list(values = .mp_pct(v), probs = probs5)

.mp_mukey_raster <- function() {
  r <- terra::rast(nrows = 1, ncols = 1, vals = 900)
  names(r) <- "mukey"
  terra::as.factor(r)
}

.mp_aoi <- function() terra::vect(terra::ext(0, 1, 0, 1), crs = "EPSG:5070")

# A recorder for cache_set() calls: returns the mock fn + an environment holding the log.
.mp_cache_recorder <- function() {
  log <- new.env(parent = emptyenv())
  log$calls <- list()
  list(
    env = log,
    fn = function(key, kind, value) {
      log$calls[[length(log$calls) + 1]] <- list(key = key, kind = kind)
      invisible(TRUE)
    }
  )
}

test_that("run_stage1_fusion_multi() runs the SSURGO simulation exactly once for N x M leaves", {
  sim_n <- 0L
  windows <- list(c(0, 5), c(5, 15), c(15, 30))
  cfgs <- list(
    ph  = list(id = "ph",  solus_variable = "ph1to1h2o", dist = "normal"),
    db  = list(id = "db",  solus_variable = "dbovendry", dist = "normal"),
    soc = list(id = "soc", solus_variable = "soc",       dist = "normal")
  )

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(..., depth_windows = NULL) {
      sim_n <<- sim_n + 1L
      stats::setNames(lapply(depth_windows, function(w) data.frame(mukey = 1)),
                      vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), ""))
    },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(...) .mp_solus(),
    cache_get_valid_percentiles = function(...) NULL,
    cache_set = function(...) invisible(TRUE),
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    .package = "soilSIM"
  )

  res <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows)

  expect_equal(sim_n, 1L)
  expect_named(res, c("ph", "db", "soc"))
  expect_named(res$ph, c("0-5", "5-15", "15-30"))
  expect_false(is.null(res$db[["5-15"]]$posterior))
  expect_true(all(c("prior", "likelihood", "posterior", "dist", "route") %in% names(res$ph[["0-5"]])))
})

test_that("run_stage1_fusion_multi() seeds the per-(id, window) 'ssurgo' and 'solus' caches", {
  rec <- .mp_cache_recorder()
  windows <- list(c(0, 5), c(5, 15))
  cfgs <- list(ph = list(id = "ph", solus_variable = "ph1to1h2o", dist = "normal"))

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(..., depth_windows = NULL) {
      stats::setNames(lapply(depth_windows, function(w) data.frame(mukey = 1)),
                      vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), ""))
    },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(...) .mp_solus(),
    cache_get_valid_percentiles = function(...) NULL,
    cache_set = rec$fn,
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    .package = "soilSIM"
  )

  run_stage1_fusion_multi(.mp_aoi(), cfgs, windows)

  kinds_by_key <- vapply(rec$env$calls, function(c) paste(c$key, c$kind), character(1))
  expect_true("ph_0_5_ssurgo ssurgo" %in% kinds_by_key)
  expect_true("ph_5_15_ssurgo ssurgo" %in% kinds_by_key)
  expect_true("ph_0_5_solus solus" %in% kinds_by_key)
  expect_true("ph_5_15_solus solus" %in% kinds_by_key)
})

test_that("run_stage1_fusion_multi() validates property_configs and composition-group inputs", {
  aoi <- .mp_aoi()
  w <- list(c(0, 5))

  expect_error(run_stage1_fusion_multi(aoi, list(list(id = "ph")), w), "named list")
  expect_error(
    run_stage1_fusion_multi(aoi, list(ph = list(id = "clay", solus_variable = "claytotal")), w),
    "\\$id must match"
  )
  expect_error(
    run_stage1_fusion_multi(aoi, list(clay = list(id = "clay", solus_variable = "claytotal",
                                                  composition_group = "texture")), w),
    "composition_groups.*is NULL"
  )
  expect_error(
    run_stage1_fusion_multi(
      aoi,
      list(clay = list(id = "clay", solus_variable = "claytotal", composition_group = "texture")),
      w,
      composition_groups = list(texture = list(members = c("clay", "sand", "silt")))
    ),
    "member\\(s\\) missing.*sand"
  )
  expect_error(run_stage1_fusion_multi(aoi, list(ph = list(id = "ph")), list(c(5, 5))), "bottom > top")
})

test_that("run_stage1_fusion_multi() returns NULL leaves on a per-leaf SOLUS failure and keeps the rest", {
  windows <- list(c(0, 5))
  cfgs <- list(
    ph = list(id = "ph", solus_variable = "ph1to1h2o", dist = "normal"),
    db = list(id = "db", solus_variable = "dbovendry", dist = "normal")
  )

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(..., depth_windows = NULL) {
      stats::setNames(lapply(depth_windows, function(w) data.frame(mukey = 1)),
                      vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), ""))
    },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(aoi_vect, solus_variable, ...) {
      if (identical(solus_variable, "dbovendry")) NULL else .mp_solus()
    },
    cache_get_valid_percentiles = function(...) NULL,
    cache_set = function(...) invisible(TRUE),
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    .package = "soilSIM"
  )

  res <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows)
  expect_false(is.null(res$ph[["0-5"]]))
  expect_null(res$db[["0-5"]])
})

test_that("run_stage1_fusion_multi() is quiet by default and honors simplify = TRUE", {
  windows <- list(c(0, 5))
  cfgs <- list(ph = list(id = "ph", solus_variable = "ph1to1h2o", dist = "normal"))

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(..., depth_windows = NULL) {
      stats::setNames(lapply(depth_windows, function(w) data.frame(mukey = 1)),
                      vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), ""))
    },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(...) .mp_solus(),
    cache_get_valid_percentiles = function(...) NULL,
    cache_set = function(...) invisible(TRUE),
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    .package = "soilSIM"
  )

  expect_silent(res <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows))

  res_s <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows, simplify = TRUE)
  expect_named(res_s$ph[["0-5"]], "percentiles")
  expect_true(all(grepl("^P", names(res_s$ph[["0-5"]]$percentiles))))
})

test_that("run_stage1_fusion_multi() skips the simulation when nothing needs it (all caches warm, no raw_draws)", {
  sim_n <- 0L
  windows <- list(c(0, 5))
  # prior_fusion_method = "percentile" -> resolve_want_raw_draws() FALSE for dist = "normal"
  cfgs <- list(ph = list(id = "ph", solus_variable = "ph1to1h2o", dist = "normal",
                         prior_fusion_method = "percentile"))

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(...) { sim_n <<- sim_n + 1L; list() },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(...) .mp_solus(),
    cache_get_valid_percentiles = function(...) .mp_prior(),  # everything warm
    cache_set = function(...) invisible(TRUE),
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    .package = "soilSIM"
  )

  res <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows)
  expect_equal(sim_n, 0L)
  expect_false(is.null(res$ph[["0-5"]]$posterior))
})

test_that("run_stage1_fusion_multi() distributes a compositional group's members and fuses the group once", {
  fuse_n <- 0L
  windows <- list(c(0, 5), c(5, 15))
  comp <- list(texture = list(members = c("clay", "sand", "silt")))
  cfgs <- list(
    clay = list(id = "clay", solus_variable = "claytotal", composition_group = "texture"),
    sand = list(id = "sand", solus_variable = "sandtotal", composition_group = "texture"),
    silt = list(id = "silt", solus_variable = "silttotal", composition_group = "texture"),
    ph   = list(id = "ph",   solus_variable = "ph1to1h2o", dist = "normal")
  )

  fake_member <- function(id) list(
    posterior = list(percentiles = .mp_pct(15)),
    dist = "texture_ilr", route = "closed_form_ilr_group",
    route_detail = NULL, n_fallback_cells = 0
  )

  testthat::local_mocked_bindings(
    fetch_ssurgo_mukey_raster = function(...) .mp_mukey_raster(),
    simulate_ssurgo_mapunit_draws = function(..., depth_windows = NULL) {
      stats::setNames(lapply(depth_windows, function(w) data.frame(mukey = 1)),
                      vapply(depth_windows, function(w) paste0(w[[1]], "-", w[[2]]), ""))
    },
    percentiles_from_draws = function(...) .mp_prior(),
    fetch_solus_percentiles = function(...) .mp_solus(),
    cache_get_valid_percentiles = function(...) NULL,
    cache_set = function(...) invisible(TRUE),
    build_cache_key = function(aoi_vect, id, top, bottom, kind) paste(id, top, bottom, kind, sep = "_"),
    mukey_draws_lookup = function(...) NULL,
    stage1_fuse_texture_group_from_fetched = function(fetched, ...) {
      fuse_n <<- fuse_n + 1L
      stats::setNames(lapply(fetched, function(f) fake_member(f$id)), vapply(fetched, `[[`, "", "id"))
    },
    .package = "soilSIM"
  )

  res <- run_stage1_fusion_multi(.mp_aoi(), cfgs, windows, composition_groups = comp)

  expect_equal(fuse_n, 2L)  # once per window, not once per member
  expect_equal(res$clay[["0-5"]]$dist, "texture_ilr")
  expect_equal(res$sand[["5-15"]]$dist, "texture_ilr")
  expect_true(all(c("prior", "likelihood", "posterior") %in% names(res$ph[["0-5"]])))  # standalone unaffected
})
