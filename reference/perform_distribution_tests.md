# Goodness-of-fit tests across a property's candidate distributions

Fits every candidate distribution (via get_appropriate_distributions() /
fit_single_distribution()) and runs fitdistrplus::gofstat() across the
successful fits, which returns KS/AD/CvM statistics (and, for the
non-censored continuous case here, approximate goodness-of-fit
decisions) for all candidates in a single call.

## Usage

``` r
perform_distribution_tests(values, property_name, config)
```

## Arguments

- values:

  Numeric vector of observed property values.

- property_name:

  Property name.

- config:

  Passed through to fit_single_distribution().

## Value

A list with `fits` (the individual fit_single_distribution() results,
best-AIC-first) and `gof` (a list of fitdistrplus::gofstat() results,
one per distinct data-scale group - e.g. beta fits are grouped
separately from unscaled fits - or NULL if no candidate fit succeeded).
