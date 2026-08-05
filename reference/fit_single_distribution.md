# Fit a single candidate distribution to a property's values

Modeled on the fitting style used in
BCQRF/code/R/distribution_fitting.R: a guarded
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) dependency,
a `tryCatch`-wrapped fit, and a structured return value rather than a
raw fit object.

## Usage

``` r
fit_single_distribution(values, dist_name, property_name, config)
```

## Arguments

- values:

  Numeric vector of observed property values.

- dist_name:

  One of "normal", "lognormal", "gamma", "beta" (the candidate names
  returned by get_appropriate_distributions()).

- property_name:

  Property name, used only for logging.

- config:

  Unused currently; kept for interface compatibility with callers.

## Value

A list with distribution/property_name/params/param_sd/loglik/aic/
bic/convergence/n/scale_info/fit_object, or NULL if fitting isn't
possible (too little data, wrong domain, or fitdist() failure).
