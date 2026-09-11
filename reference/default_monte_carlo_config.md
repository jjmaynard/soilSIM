# Get Monte Carlo Default Configuration

Returns default configuration optimized for Monte Carlo simulation

## Usage

``` r
default_monte_carlo_config(verbose = getOption("ssurgo.verbose", FALSE))
```

## Arguments

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Default Monte Carlo configuration

## See also

Other monte-carlo:
[`simulate_monte_carlo()`](https://jjmaynard.github.io/soilSIM/reference/simulate_monte_carlo.md)
