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
