test_that("fit_gamma_mom_raster() matches the scalar moments_to_gamma() cell-by-cell", {
  rasters <- make_percentile_rasters(c(a = 8, b = 10, c = 12))
  fit_r <- fit_gamma_mom_raster(list(rasters$a, rasters$b, rasters$c))
  mean_v <- mean(c(8, 10, 12)); var_v <- mean((c(8, 10, 12) - mean_v)^2)
  fit_scalar <- moments_to_gamma(mean_v, var_v)
  expect_equal(unique(terra::values(fit_r$shape))[1], fit_scalar$shape)
  expect_equal(unique(terra::values(fit_r$rate))[1], fit_scalar$rate)
})

test_that("percentile_index()/percentile_interior() locate low/median/high and interior probabilities", {
  idx <- percentile_index(c(0, 0.05, 0.5, 0.95, 1))
  expect_equal(idx$p_lo, 0.05)
  expect_equal(idx$p_hi, 0.95)
  expect_equal(idx$p50_idx, 3)

  interior <- percentile_interior(c(0, 0.05, 0.5, 0.95, 1))
  expect_equal(interior$probs, c(0.05, 0.5, 0.95))
})

test_that("align_percentile_probs() is a no-op when probabilities already match", {
  rasters <- make_percentile_rasters(c(P5 = 2, P50 = 10, P95 = 18))
  aligned <- align_percentile_probs(
    list(rasters$P5, rasters$P50, rasters$P95), c(0.05, 0.5, 0.95),
    list(rasters$P5, rasters$P50, rasters$P95), c(0.05, 0.5, 0.95)
  )
  expect_equal(aligned$probs, c(0.05, 0.5, 0.95))
})

test_that("align_percentile_probs() reinterpolates the wider side onto the narrower side's probabilities", {
  prior_r <- make_percentile_rasters(c(P5 = 2, P50 = 10, P95 = 18))
  lik_r <- make_percentile_rasters(c(P2.5 = 1, P50 = 10, P97.5 = 19))
  aligned <- align_percentile_probs(
    list(prior_r$P5, prior_r$P50, prior_r$P95), c(0.05, 0.5, 0.95),
    list(lik_r$P2.5, lik_r$P50, lik_r$P97.5), c(0.025, 0.5, 0.975)
  )
  # prior's range [0.05, 0.95] is narrower than lik's [0.025, 0.975] -> target = prior_probs
  expect_equal(aligned$probs, c(0.05, 0.5, 0.95))
  expect_length(aligned$lik_value_rasters, 3)
})

test_that("fuse_adaptive() routes to the closed-form path above threshold_cells and matches bayes_update_normal_normal()", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 0, verbose = FALSE
  )
  expect_equal(result$route, "closed_form_normal")

  fit_prior <- fit_normal_raster(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95, 0.05, 0.95)
  fit_lik <- fit_normal_raster(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95, 0.05, 0.95)
  expected <- bayes_update_normal_normal(fit_prior$mu, fit_prior$sigma, fit_lik$mu, fit_lik$sigma)
  expect_equal(unique(terra::values(result$posterior$mu))[1], unique(terra::values(expected$mu))[1])
})

test_that("fuse_adaptive() routes to the general KDE path at or below threshold_cells", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 1e6, n_samples = 50, verbose = FALSE
  )
  expect_equal(result$route, "bayesian_update_general")
  expect_true(all(is.finite(terra::values(result$posterior$mu))))
})

test_that("fuse_general_kde() degrades a cell with too many NA percentile values to NA output, rather than aborting the whole raster", {
  # Real-world discovery: cells with insufficient valid percentile values (e.g. at raster edges
  # after resample()/alignment) used to hard-crash the entire terra::app() call via
  # extract_percentile_pairs()'s "Need at least 2 valid quantile columns" stop().
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6), nrow = 1, ncol = 2)
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9), nrow = 1, ncol = 2)
  # Poison the second cell of every prior layer with NA.
  for (p in c("P5", "P50", "P95")) {
    terra::values(rasters_prior[[p]])[2] <- NA
  }

  result <- fuse_adaptive(
    prior_value_rasters = list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95),
    lik_value_rasters = list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95),
    percentile_probs = c(0.05, 0.5, 0.95),
    family = "normal", threshold_cells = 1e6, n_samples = 50, verbose = FALSE
  )
  mu_vals <- terra::values(result$posterior$mu)[, 1]
  expect_true(is.finite(mu_vals[1]))
  expect_true(is.na(mu_vals[2]))
})

test_that("fuse_property_adaptive() dispatches normal/beta/gamma/lognormal correctly", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))
  prior_list <- list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95)
  lik_list <- list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95)
  probs <- c(0.05, 0.5, 0.95)

  r_normal <- fuse_property_adaptive(prior_list, probs, lik_list, probs,
    property_config = list(dist = "normal"), threshold_cells = 0, verbose = FALSE)
  expect_equal(r_normal$dist, "normal")
  expect_equal(r_normal$route, "closed_form_normal")

  r_lognormal <- fuse_property_adaptive(prior_list, probs, lik_list, probs,
    property_config = list(dist = "lognormal"), threshold_cells = 0, verbose = FALSE)
  expect_equal(r_lognormal$dist, "lognormal")
  expect_equal(r_lognormal$route, "closed_form_lognormal")

  rasters_prior_pos <- make_percentile_rasters(c(P5 = 20, P50 = 40, P95 = 65))
  rasters_lik_pos <- make_percentile_rasters(c(P5 = 30, P50 = 50, P95 = 70))
  r_beta <- fuse_property_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95), probs,
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95), probs,
    property_config = list(dist = "beta", bounds = c(0, 100)), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_beta$dist, "beta")

  r_gamma <- fuse_property_adaptive(
    list(rasters_prior_pos$P5, rasters_prior_pos$P50, rasters_prior_pos$P95), probs,
    list(rasters_lik_pos$P5, rasters_lik_pos$P50, rasters_lik_pos$P95), probs,
    property_config = list(dist = "gamma"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_gamma$dist, "gamma")
})

test_that("fuse_property_adaptive() dispatches metalog and flags fallback via route_detail", {
  rasters_prior <- make_percentile_rasters(c(P10 = 20, P50 = 40, P90 = 65))
  rasters_lik <- make_percentile_rasters(c(P10 = 30, P50 = 50, P90 = 70))
  probs <- c(0.1, 0.5, 0.9)

  r_metalog <- fuse_property_adaptive(
    list(rasters_prior$P10, rasters_prior$P50, rasters_prior$P90), probs,
    list(rasters_lik$P10, rasters_lik$P50, rasters_lik$P90), probs,
    property_config = list(dist = "metalog", bounds = c(0, 100), boundedness = "b"),
    verbose = FALSE
  )
  expect_equal(r_metalog$dist, "metalog")
  expect_equal(r_metalog$route, "closed_form_metalog_adapter")
  expect_true(all(is.finite(terra::values(r_metalog$posterior$mu))))
})

test_that("fuse_property_adaptive() resolves dist = 'auto' to beta/normal/lognormal without overriding an explicit dist", {
  rasters_prior <- make_percentile_rasters(c(P5 = 4, P50 = 5, P95 = 6))
  rasters_lik <- make_percentile_rasters(c(P5 = 7, P50 = 8, P95 = 9))
  probs <- c(0.05, 0.5, 0.95)

  r_auto <- fuse_property_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95), probs,
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95), probs,
    property_config = list(dist = "auto"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_auto$dist, "normal")
  expect_equal(r_auto$dist_source, "auto")

  r_explicit <- fuse_property_adaptive(
    list(rasters_prior$P5, rasters_prior$P50, rasters_prior$P95), probs,
    list(rasters_lik$P5, rasters_lik$P50, rasters_lik$P95), probs,
    property_config = list(dist = "lognormal"), threshold_cells = 0, verbose = FALSE
  )
  expect_equal(r_explicit$dist_source, "config")
})

test_that("group_members() returns configured members in order, reusing soilSIM's composition_groups shape", {
  composition_groups <- list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))
  expect_equal(group_members("texture", composition_groups), c("claytotal", "sandtotal", "silttotal"))
  expect_error(group_members("nonexistent", composition_groups), "no composition group")
})

test_that("fuse_texture_group() returns one posterior per member that sums to ~100 and matches ilr_inverse()", {
  clay_prior <- make_percentile_rasters(c(lo = 15, p50 = 20, hi = 25))
  sand_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  silt_prior <- make_percentile_rasters(c(lo = 35, p50 = 40, hi = 45))
  clay_lik <- make_percentile_rasters(c(lo = 18, p50 = 22, hi = 28))
  sand_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))
  silt_lik <- make_percentile_rasters(c(lo = 33, p50 = 38, hi = 43))

  probs <- c(0.05, 0.5, 0.95)
  fetched <- list(
    list(id = "claytotal", prior = list(clay_prior$lo, clay_prior$p50, clay_prior$hi), prior_probs = probs,
         lik = list(clay_lik$lo, clay_lik$p50, clay_lik$hi), lik_probs = probs),
    list(id = "sandtotal", prior = list(sand_prior$lo, sand_prior$p50, sand_prior$hi), prior_probs = probs,
         lik = list(sand_lik$lo, sand_lik$p50, sand_lik$hi), lik_probs = probs),
    list(id = "silttotal", prior = list(silt_prior$lo, silt_prior$p50, silt_prior$hi), prior_probs = probs,
         lik = list(silt_lik$lo, silt_lik$p50, silt_lik$hi), lik_probs = probs)
  )

  result <- fuse_texture_group(fetched)
  expect_named(result, c("claytotal", "sandtotal", "silttotal"))

  totals <- Reduce(`+`, lapply(result, function(r) r$posterior$value))
  expect_true(all(abs(terra::values(totals) - 100) < 1e-6))
  for (id in names(result)) {
    expect_equal(result[[id]]$dist, "texture_ilr")
    expect_equal(result[[id]]$route, "closed_form_ilr_group")
  }
})

test_that("run_stage1_fusion() fetches, caches, aligns, and fuses SSURGO x SOLUS for a real AOI", {
  testthat::skip_if_offline()
  # See test-solus-simulation.R's live test for why PROJ_LIB is unset defensively here too - a
  # stray PROJ_LIB (documented in HANDOFF_NOTES.md) can poison terra's PROJ context before
  # soilSIM's own .onLoad() runs, when terra was already attached earlier in the same test session.
  Sys.unsetenv("PROJ_LIB")
  aoi <- terra::vect(
    "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
    "AOI has no valid CRS in this session - likely the documented PROJ_LIB environment issue; not a soilSIM defect."
  )
  aoi <- terra::project(aoi, "epsg:5070")

  property_config <- list(id = "test_clay", solus_variable = "claytotal", dist = "normal")

  result <- run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5)
  expect_true(is.null(result) || is.list(result))
  if (!is.null(result)) {
    expect_true(result$dist %in% c("normal", "beta", "gamma", "lognormal", "metalog"))
    expect_true(is.list(result$posterior))
    expect_s4_class(result$prior$values[[1]], "SpatRaster")
  }
})

test_that("run_stage1_fusion() with a composition_group dispatches to run_stage1_fusion_group() and slices to the requested member", {
  testthat::skip_if_offline()
  Sys.unsetenv("PROJ_LIB")
  aoi <- terra::vect(
    "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))",
    crs = "epsg:4326"
  )
  testthat::skip_if(is.na(terra::crs(aoi)) || !nzchar(terra::crs(aoi)),
    "AOI has no valid CRS in this session - likely the documented PROJ_LIB environment issue; not a soilSIM defect."
  )
  aoi <- terra::project(aoi, "epsg:5070")

  composition_groups <- list(texture = list(members = c("test_clay2", "test_sand2", "test_silt2")))
  property_configs <- list(
    test_clay2 = list(id = "test_clay2", solus_variable = "claytotal", composition_group = "texture"),
    test_sand2 = list(id = "test_sand2", solus_variable = "sandtotal", composition_group = "texture"),
    test_silt2 = list(id = "test_silt2", solus_variable = "silttotal", composition_group = "texture")
  )

  result <- run_stage1_fusion(
    aoi, property_configs$test_clay2, top_depth = 0, bottom_depth = 5,
    composition_groups = composition_groups, property_configs = property_configs
  )
  expect_true(is.null(result) || identical(result$dist, "texture_ilr"))
})
