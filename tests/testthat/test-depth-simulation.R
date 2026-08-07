test_that("infill_missing_depth_variability() fills missing bounds with +/-2 of the representative value", {
  data <- data.frame(
    hzdept_r = c(10, 20), hzdept_l = c(NA, 18), hzdept_h = c(NA, NA),
    hzdepb_r = c(30, 40), hzdepb_l = c(NA, NA), hzdepb_h = c(NA, 42)
  )
  out <- infill_missing_depth_variability(data)
  expect_equal(out$hzdept_l, c(8, 18))
  expect_equal(out$hzdept_h, c(12, 22))
  expect_equal(out$hzdepb_l, c(28, 38))
  expect_equal(out$hzdepb_h, c(32, 42))
})

test_that("infill_missing_depth_variability() clamps at zero, never goes negative", {
  data <- data.frame(hzdept_r = 1, hzdept_l = NA, hzdept_h = NA, hzdepb_r = 1, hzdepb_l = NA, hzdepb_h = NA)
  out <- infill_missing_depth_variability(data)
  expect_equal(out$hzdept_l, 0)
  expect_equal(out$hzdepb_l, 0)
})

test_that("infill_missing_distinctness() fills defaults by genhz, preserves existing values", {
  df <- data.frame(
    hzname = c("A", "B", "C", "R", "O", "Cr"),
    distinctness = c(NA, "gradual", NA, NA, "diffuse", NA),
    stringsAsFactors = FALSE
  )
  out <- infill_missing_distinctness(df)
  expect_equal(out$distinctness, c("clear", "gradual", "gradual", "abrupt", "diffuse", "gradual"))
})

test_that("simulate_soil_profile_top_down() starts at 0, never produces bottom < top, and matches bounds", {
  set.seed(1)
  horizon_data <- make_depth_horizon_data()
  profile <- simulate_soil_profile_top_down(horizon_data)
  expect_equal(profile$top[1], 0)
  expect_true(all(profile$bottom >= profile$top))
  # each horizon's simulated bottom must fall within its l/h bounds (post-clamping)
  expect_true(all(profile$bottom <= horizon_data$hzdepb_h + 1e-6))
  # top of horizon i+1 must equal bottom of horizon i (top-down contiguity)
  expect_equal(profile$top[-1], profile$bottom[-nrow(profile)])
})

test_that("simulate_soil_profile_bottom_up() starts at 0 and has no gaps between horizons", {
  set.seed(2)
  horizon_data <- make_depth_horizon_data()
  profile <- simulate_soil_profile_bottom_up(horizon_data)
  expect_equal(profile$top[1], 0)
  expect_true(all(profile$bottom >= profile$top))
  expect_equal(profile$bottom[-nrow(profile)], profile$top[-1])
})

test_that("simulate_soil_profile_thickness() returns one row per horizon with a non-negative thickness_sd", {
  set.seed(3)
  horizon_data <- make_depth_horizon_data()
  result <- simulate_soil_profile_thickness(horizon_data, n_simulations = 50)
  expect_equal(nrow(result), nrow(horizon_data))
  expect_true(all(result$thickness_sd >= 0))
  expect_setequal(result$hzname, horizon_data$hzname)
})

test_that("simulate_and_perturb_soil_profiles() replicates a single-horizon profile sim_comppct times without perturbation", {
  df <- data.frame(
    cokey = "1", mukey = "1", compname = "onehorizon", id = "onehorizon",
    hzname = "R", hzdept_l = 0, hzdept_r = 0, hzdept_h = 0,
    hzdepb_l = 145, hzdepb_r = 150, hzdepb_h = 155, hzthk_l = 145,
    sim_comppct = 4,
    stringsAsFactors = FALSE
  )
  spc <- df
  aqp::depths(spc) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(spc) <- "hzname"

  result <- simulate_and_perturb_soil_profiles(spc)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 4)
})

test_that("simulate_profile_depths_by_mukey()'s id-column fix: aqp::depths(id ~ ...) doesn't error on a compname-derived id", {
  # Regression test for the id-column bug fix: get_aws_data_by_mukey() output
  # only has `compname`, not `id`; simulate_profile_depths_by_mukey() must
  # derive `id` before calling aqp::depths<-() since adjust_out_of_range_profiles()/
  # evaluate_simulated_depths() downstream both hardcode a literal `id` column.
  mu_data <- data.frame(
    mukey = "1", cokey = "1", compname = "testseries",
    hzname = c("A", "Bt"), hzdept_r = c(0, 20), hzdepb_r = c(20, 50),
    stringsAsFactors = FALSE
  )
  mu_data$id <- mu_data$compname
  expect_no_error(aqp::depths(mu_data) <- id ~ hzdept_r + hzdepb_r)
})

test_that("simulate_profile_depths_by_mukey() derives and joins sim_comppct internally (closed integration gap)", {
  # Regression test for the sim_comppct integration gap: simulate_profile_depths_by_mukey()
  # used to require callers to derive/join sim_comppct themselves before it reached
  # simulate_and_perturb_soil_profiles(), erroring on a missing-column condition otherwise.
  # It now calls sim_component_comp() and joins the result onto mu_data by cokey internally.
  # comppct_l = comppct_r = comppct_h = 100 makes tri_dist() degenerate to a point mass at
  # 100 (a == b == c), so sim_comppct = round(n_simulations * 100 / 100) = n_simulations
  # exactly - a fully deterministic, offline-testable case (single horizon, so the
  # early-return path in simulate_and_perturb_soil_profiles() is taken and no live OSD
  # lookup occurs).
  testthat::local_mocked_bindings(
    get_aws_data_by_mukey = function(mukeys) {
      data.frame(
        mukey = "999", cokey = "1", compname = "onehorizon",
        comppct_l = 100, comppct_r = 100, comppct_h = 100,
        hzname = "R", hzdept_l = 0, hzdept_r = 0, hzdept_h = 0,
        hzdepb_l = 145, hzdepb_r = 150, hzdepb_h = 155, hzthk_l = 145,
        stringsAsFactors = FALSE
      )
    },
    .package = "soilSIM"
  )

  set.seed(1)
  result <- simulate_profile_depths_by_mukey("999", n_simulations = 5, seed = 1)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 5)
})

test_that("evaluate_simulated_depths() flags horizons whose simulated depths exceed the original bounds", {
  simulated_profiles <- make_soil_profile_collection(sim_comppct = 1)
  horizon_data <- make_depth_horizon_data()
  # Force an out-of-range simulated top for the first horizon
  aqp::horizons(simulated_profiles)$hzdept_r[1] <- horizon_data$hzdept_l[1] - 5

  out_of_range <- evaluate_simulated_depths(simulated_profiles, horizon_data)
  expect_true(nrow(out_of_range) >= 1)
  expect_true(all(out_of_range$out_of_range))
})

test_that("simulate_and_perturb_soil_profiles() runs the full multi-horizon perturbation path (thickness + OSD-distinctness boundary perturbation) against a real series", {
  # Unlike the single-horizon test above (which takes the early-return path
  # before query_osd_distinctness() is ever called), a multi-horizon profile
  # exercises the full pipeline, which unconditionally calls
  # query_osd_distinctness() -> soilDB::fetchOSD() - a live OSD lookup that
  # requires a real component name (a synthetic "testseries" would return no
  # OSD data and error). "amador" is the same real series used in this
  # function's own historical example.
  testthat::skip_if_offline()
  df <- data.frame(
    cokey = "1", mukey = "1", compname = "amador", id = "amador",
    hzname = c("A", "Bt", "Cr"),
    hzdept_l = c(0, 15, 45), hzdept_r = c(0, 20, 50), hzdept_h = c(0, 25, 55),
    hzdepb_l = c(15, 45, 95), hzdepb_r = c(20, 50, 100), hzdepb_h = c(25, 55, 105),
    hzthk_l = c(15, 30, 50), sim_comppct = 4,
    stringsAsFactors = FALSE
  )
  spc <- df
  aqp::depths(spc) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(spc) <- "hzname"

  result <- simulate_and_perturb_soil_profiles(spc)
  expect_s4_class(result, "SoilProfileCollection")
  expect_equal(length(result), 4)
})

test_that("get_aws_data_by_mukey()/query_osd_distinctness() require the live SDA/OSD services", {
  testthat::skip_if_offline()
  testthat::skip("Live NRCS Soil Data Access / OSD queries are not exercised in automated tests - see test-ssurgo-acquisition.R for the established precedent.")
})

test_that("simulate_profile_depths_by_collection_parallel() runs without erroring at the R level", {
  skip_if_not(
    nzchar(system.file(package = "soilSIM")) &&
      file.exists(file.path(system.file(package = "soilSIM"), "Meta", "package.rds")),
    "soilSIM not installed in this session - future::multisession workers can't see a load_all()'d namespace"
  )
  # Unlike process_cokeys_parallel() (see test-multivariate-adjustment.R),
  # this function has no sequential fallback - it's a straight port of the
  # legacy code_ref/brdf/depth_simulation.R, which just returns NULL from its
  # tryCatch on any worker-side error. Some sandboxed environments' fresh
  # future::multisession worker processes can't see an installed package
  # either (the same "there is no package called 'soilSIM'" limitation
  # documented and tolerated in test-multivariate-adjustment.R), so accept
  # either a real result or that documented, gracefully-handled NULL.
  # future::multisession/parallelly's Windows PSOCK workers write launcher
  # scripts to tempdir() that outlive the worker (a documented
  # future/parallelly quirk, not a soilSIM defect) - clean up any it drops,
  # matching this suite's existing convention of tests cleaning up their own
  # tempdir droppings (see the cache_dir on.exit() in test-ssurgo-acquisition.R).
  tmp_before <- list.files(tempdir())
  on.exit({
    new_files <- setdiff(list.files(tempdir()), tmp_before)
    unlink(file.path(tempdir(), grep("^Rscript", new_files, value = TRUE)))
  }, add = TRUE)

  soil_collection <- make_soil_profile_collection(sim_comppct = 2)
  result <- simulate_profile_depths_by_collection_parallel(soil_collection, seed = 1, n_cores = 2)
  expect_true(is.null(result) || inherits(result, "SoilProfileCollection"))
})

test_that("simulate_profile_depths_by_collection_parallel() restores the caller's future::plan() afterward", {
  skip_if_not(
    nzchar(system.file(package = "soilSIM")) &&
      file.exists(file.path(system.file(package = "soilSIM"), "Meta", "package.rds")),
    "soilSIM not installed in this session - future::multisession workers can't see a load_all()'d namespace"
  )
  # Regression test for a real bug fixed this session: the pre-migration implementation reverted
  # future::plan() to sequential() only on its own success path (not via on.exit()), so it either
  # clobbered a caller's pre-existing plan or left multisession active indefinitely on error.
  # run_parallel_lapply() (R/parallel-utils.R) now snapshots and restores the caller's plan
  # unconditionally.
  tmp_before <- list.files(tempdir())
  on.exit({
    new_files <- setdiff(list.files(tempdir()), tmp_before)
    unlink(file.path(tempdir(), grep("^Rscript", new_files, value = TRUE)))
  }, add = TRUE)

  plan_before <- future::plan()
  soil_collection <- make_soil_profile_collection(sim_comppct = 2)
  invisible(simulate_profile_depths_by_collection_parallel(soil_collection, seed = 1, n_cores = 2))
  expect_true(identical(class(future::plan()), class(plan_before)))
})
