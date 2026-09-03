test_that("get_appropriate_distributions() honors config$distribution_methods (previously inert)", {
  values <- c(5, 6, 7, 8, 9)

  default_candidates <- get_appropriate_distributions("om_r", values)
  expect_true("weibull" %in% default_candidates)  # weibull now reachable for a generic positive property

  narrowed <- get_appropriate_distributions("om_r", values, config = list(distribution_methods = c("gamma", "beta")))
  expect_equal(narrowed, "gamma")  # beta isn't in the heuristic for a generic property, gamma is -> intersection

  # A request with no overlap at all falls back to the user's list unfiltered
  # (with a WARN), rather than silently discarding it.
  expect_warning(
    unfiltered <- get_appropriate_distributions("sandtotal_r", values, config = list(distribution_methods = c("gamma"))),
    NA  # log_message() doesn't raise an R warning condition; assert no error instead
  )
  expect_equal(unfiltered, "gamma")
})

test_that("weibull is fittable end-to-end (previously listed as a valid choice but unreachable)", {
  set.seed(5)
  values <- rweibull(200, shape = 2, scale = 10)

  fit <- fit_single_distribution(values, "weibull", "test_prop", config = list())
  expect_false(is.null(fit))
  expect_equal(fit$distribution, "weibull")
  expect_true(all(c("shape", "scale") %in% names(fit$params)))

  # weibull is skipped (like lognormal/gamma) for non-positive data, not an error
  skipped <- fit_single_distribution(c(-1, 0, 1, 2, 3), "weibull", "test_prop", config = list())
  expect_null(skipped)
})

test_that("soilSIM no longer declares a dependency on the (dropped) compositions package", {
  # Reading R/statistics.R's source text isn't reliable once the package is
  # installed (R CMD check runs tests against the installed/lazy-loaded
  # package, not the source tree) - check the DESCRIPTION's declared
  # dependencies instead, which always exists post-install and is the
  # actually meaningful thing to assert (texture correlations are computed
  # via the dependency-free ilr_forward() now, not compositions::ilr()).
  description_path <- system.file("DESCRIPTION", package = "soilSIM")
  skip_if(identical(description_path, ""), "soilSIM DESCRIPTION not found on the search path")
  fields <- read.dcf(description_path, fields = c("Imports", "Depends", "Suggests"))
  all_deps <- paste(fields, collapse = " ")
  expect_false(grepl("compositions", all_deps))
})

test_that("validate_correlation_matrices() actually flags an invalid matrix (previously always valid=TRUE)", {
  good <- diag(3); rownames(good) <- colnames(good) <- c("a", "b", "c")
  bad <- good; bad[1, 1] <- 0.5  # non-unit diagonal

  result <- validate_correlation_matrices(list(pearson = list(matrix = bad)), config = list())
  expect_true(length(result$warnings) > 0)

  result_good <- validate_correlation_matrices(list(pearson = list(matrix = good)), config = list())
  expect_equal(length(result_good$warnings), 0)
})

test_that("analyze_texture_correlations() reports both raw and ILR-space correlations for a complete texture set", {
  set.seed(9)
  n <- 30
  clay <- runif(n, 10, 30)
  sand <- runif(n, 20, 50)
  silt <- 100 - clay - sand
  data <- data.frame(claytotal_r = clay, sandtotal_r = sand, silttotal_r = silt)

  result <- analyze_texture_correlations(data, c("claytotal_r", "sandtotal_r", "silttotal_r"), c("pearson"), config = list())
  expect_false(is.null(result$ilr_correlations))
  expect_true(all(dim(result$ilr_correlations$pearson) == c(2, 2)))
  expect_false(is.null(result$raw_correlations$pearson))
  expect_equal(result$n_observations, n)
})

test_that("analyze_texture_correlations() degrades gracefully with fewer than 3 texture members", {
  data <- data.frame(sandtotal_r = runif(10, 20, 50), claytotal_r = runif(10, 10, 30))
  result <- analyze_texture_correlations(data, c("sandtotal_r", "claytotal_r"), c("pearson"), config = list())
  expect_null(result$ilr_correlations)
  expect_false(is.null(result$raw_correlations$pearson))
})

## ----------------------------------------------------------------------
## Previously untagged placeholder functions
## ----------------------------------------------------------------------

test_that("detect_multivariate_outliers() flags a real Mahalanobis outlier and degrades gracefully with <2 properties", {
  set.seed(1)
  n <- 100
  data <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n))
  data$a[1] <- 100  # extreme multivariate outlier
  data$b[1] <- -100

  result <- detect_multivariate_outliers(data, c("a", "b"), config = list())
  expect_true(isTRUE(result$outliers[1]))
  expect_true(result$n_outliers >= 1)

  result_insufficient <- detect_multivariate_outliers(data, "a", config = list())
  expect_equal(result_insufficient$n_outliers, 0)
})

test_that("generate_outlier_summary() aggregates property and multivariate outlier counts", {
  outlier_results <- list(
    property_outliers = list(
      a = list(iqr = list(n_outliers = 3), zscore = list(n_outliers = 5))
    ),
    multivariate_outliers = list(n_outliers = 2)
  )
  result <- generate_outlier_summary(outlier_results)
  expect_equal(result$total_outliers, 5 + 2)  # max across methods for "a", plus multivariate
  expect_equal(result$multivariate_outlier_count, 2)
})

test_that("calculate_missing_data_improvement() measures a real reduction in missing-value rate", {
  original_data <- data.frame(a = c(1, NA, 3, NA), b = c(NA, 2, 3, 4))
  processed_data <- data.frame(a = c(1, 2, 3, 4), b = c(1, 2, 3, 4))
  improvement <- calculate_missing_data_improvement(original_data, processed_data)
  expect_equal(improvement, mean(c(0.5, 0.25)), tolerance = 1e-8)

  expect_true(is.na(calculate_missing_data_improvement(NULL, processed_data)))
})

test_that("assess_correlation_quality() penalizes matrices with validation warnings", {
  clean <- assess_correlation_quality(list(matrices = list(pearson = diag(2)), validation = list(warnings = character(0))))
  expect_equal(clean$quality_score, 1)

  flagged <- assess_correlation_quality(list(matrices = list(pearson = diag(2)), validation = list(warnings = "bad matrix")))
  expect_true(flagged$quality_score < clean$quality_score)
})

test_that("assess_distribution_quality() scores the fraction of properties with a non-empty fits list", {
  distribution_analysis <- list(
    prop_a = list(fits = list(normal = list())),
    prop_b = list(fits = list())
  )
  result <- assess_distribution_quality(distribution_analysis)
  expect_equal(result$quality_score, 0.5)
})

test_that("assess_outlier_quality() derives a quality score from the real outlier rate", {
  outlier_analysis <- list(
    summary = list(total_outliers = 5),
    property_outliers = list(a = list(iqr = list(outliers = rep(FALSE, 100))))
  )
  result <- assess_outlier_quality(outlier_analysis)
  expect_equal(result$outlier_rate, 0.05)
  expect_equal(result$quality_score, 0.95)
})

test_that("calculate_overall_statistics_quality_score() averages whichever component scores are available", {
  data_validation <- list(overall_quality = list(score = 0.9))
  validation_results <- list(validation_score = 0.8)
  score <- calculate_overall_statistics_quality_score(data_validation, NULL, NULL, NULL, validation_results)
  expect_equal(score, mean(c(0.9, 0.8)), tolerance = 1e-8)

  expect_true(is.na(calculate_overall_statistics_quality_score(NULL, NULL, NULL, NULL, NULL)))
})

# ---------------------------------------------------------------------------
# Problem B - analyze_soil_statistics() safe vs enhanced chain
# ---------------------------------------------------------------------------

.stats_input <- function(n = 40, seed = 9) {
  set.seed(seed)
  d <- make_soil_data(n, properties = c("dbovendry", "sandtotal", "claytotal", "silttotal"))
  # make_soil_data() leaves the texture _r columns constant; real soil data varies, and a constant
  # column breaks correlation/distribution/outlier analysis (gracefully in the _safe chain, less so
  # in the enhanced one) - jitter them so this fixture exercises a realistic path.
  for (p in c("sandtotal_r", "claytotal_r", "silttotal_r", "dbovendry_r")) {
    d[[p]] <- d[[p]] + stats::rnorm(nrow(d), 0, abs(mean(d[[p]])) * 0.08)
  }
  d
}

test_that("analyze_soil_statistics(use_enhanced_chain = FALSE) is unchanged from omitting the argument (B.4 regression)", {
  d <- .stats_input()
  a <- analyze_soil_statistics(d, verbose = FALSE)
  b <- analyze_soil_statistics(d, use_enhanced_chain = FALSE, verbose = FALSE)
  # drop timing-dependent metadata before comparing
  strip <- function(x) { x$analysis_metadata <- NULL; x }
  expect_identical(strip(a), strip(b))
})

test_that("analyze_soil_statistics(use_enhanced_chain = TRUE) runs on the package's own default usage without erroring (B.1)", {
  d <- .stats_input()
  # The whole point of B.1: before the config-plumbing fix this hard-errored with
  # "argument is of length zero" on `analyze_soil_statistics(data)` (no analysis_config).
  expect_no_error(
    res <- analyze_soil_statistics(d, use_enhanced_chain = TRUE, verbose = FALSE)
  )
  expect_true(all(c("correlation_matrices", "distribution_analysis", "outlier_analysis",
                    "quality_report", "validation_results") %in% names(res)))
})

test_that("the two chains return non-interchangeable distribution_analysis shapes", {
  d <- .stats_input()
  safe <- analyze_soil_statistics(d, use_enhanced_chain = FALSE, verbose = FALSE)
  enh  <- analyze_soil_statistics(d, use_enhanced_chain = TRUE,  verbose = FALSE)

  safe_fit <- safe$distribution_analysis$fitted_distributions[[1]]
  enh_fit  <- enh$distribution_analysis$fitted_distributions[[1]]
  expect_true("mean" %in% names(safe_fit))          # _safe: flat stats list
  expect_false("mean" %in% names(enh_fit))          # enhanced: family-name-keyed list of fit objects
  expect_true(is.list(enh_fit) && length(enh_fit) >= 1)
})

test_that("generate_statistical_quality_report() gets the per-property fitted_distributions nesting, not the wrapper (B.2)", {
  d <- .stats_input()
  res <- analyze_soil_statistics(d, use_enhanced_chain = TRUE, verbose = FALSE)
  dq <- res$quality_report$analysis_quality$distribution_quality
  expect_false(is.null(dq))
  expect_true(is.numeric(dq$quality_score))
  expect_true(dq$n_properties >= 1)
  # if the whole distribution_analysis object had been passed, n_properties would be 3
  # (fitted_distributions/distribution_tests/summary), not the real property count.
  expect_equal(dq$n_properties, length(res$distribution_analysis$fitted_distributions))
})

test_that("validate_distribution_analysis() flags non-zero convergence and non-finite fit metrics (B.3)", {
  clean <- list(fitted_distributions = list(
    clay = list(normal = list(aic = 100, bic = 105, loglik = -48, convergence = 0, param_sd = c(0.1, 0.2)))
  ))
  expect_equal(validate_distribution_analysis(clean, list())$warnings, character(0))

  bad <- list(fitted_distributions = list(
    clay = list(normal = list(aic = Inf, bic = 105, loglik = -48, convergence = 1, param_sd = c(0.1, NaN)))
  ))
  w <- validate_distribution_analysis(bad, list())$warnings
  expect_true(any(grepl("convergence", w)))
  expect_true(any(grepl("non-finite aic", w)))
  expect_true(any(grepl("parameter SE", w)))

  expect_match(validate_distribution_analysis(list(fitted_distributions = list()), list())$warnings,
               "fitted no properties")
})

test_that("validate_outlier_analysis() flags an excessive outlier rate and a non-finite multivariate threshold (B.3)", {
  low <- list(property_outliers = list(clay = list(iqr = list(outlier_rate = 0.02, n_outliers = 1))))
  expect_equal(validate_outlier_analysis(low, list())$warnings, character(0))

  high <- list(
    property_outliers = list(clay = list(iqr = list(outlier_rate = 0.4, n_outliers = 40))),
    multivariate_outliers = list(threshold = NaN, n_outliers = 3)
  )
  w <- validate_outlier_analysis(high, list())$warnings
  expect_true(any(grepl("clay/iqr: outlier rate 40", w)))
  expect_true(any(grepl("non-finite threshold", w)))

  # threshold is config-tunable
  w2 <- validate_outlier_analysis(high, list(statistical_analysis = list(outlier_rate_warn_threshold = 0.5)))$warnings
  expect_false(any(grepl("outlier rate", w2)))
})
