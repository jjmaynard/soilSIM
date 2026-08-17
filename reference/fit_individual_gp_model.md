# Fit Individual GP Model

Enhanced version with better error handling and diagnostics

## Usage

``` r
fit_individual_gp_model(
  data,
  property,
  optimize_hyperparameters = TRUE,
  gp_control = c(20, 10, 2),
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Input data for the group

- property:

  Property name to model

- optimize_hyperparameters:

  Whether to optimize hyperparameters

- gp_control:

  Passed through to
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s
  `control` argument (population size / iteration counts for its
  internal hyperparameter search). Defaults to `c(20, 10, 2)`, far below
  `GP_fit()`'s own default `c(200*d, 80*d, 2*d)` - see
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)'s
  docs (`multivariate-adjustment.R`) for the empirical justification
  (identical fitted results, ~5x faster per call). Only used when
  `optimize_hyperparameters = TRUE` calls through to
  [`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)/
  [`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
  (many `GP_fit()` calls per fit - 5-fold CV x 3 correlation
  candidates + a refit), or by the direct `GP_fit()` call when
  `optimize_hyperparameters = FALSE`.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List containing GP model and diagnostics
