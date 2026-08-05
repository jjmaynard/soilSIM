# Generate Monte Carlo Realizations of Soil Properties

Master function for generating Monte Carlo realizations for validation,
error handling, logging, and data quality assessment.

## Usage

``` r
generate_monte_carlo_realizations(
  soil_data,
  properties,
  correlation_matrix = NULL,
  n_realizations = 1000,
  simulation_config = list(),
  validate_inputs = TRUE,
  parallel = FALSE,
  seed = NULL,
  observed_data = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- soil_data:

  Data frame with soil property data (must have \_l, \_r, \_h columns)

- properties:

  Character vector of properties to simulate

- correlation_matrix:

  Optional correlation matrix for properties

- n_realizations:

  Integer, number of Monte Carlo realizations (default = 1000)

- simulation_config:

  List of simulation configuration parameters

- validate_inputs:

  Logical, whether to validate inputs (default = TRUE)

- parallel:

  Logical, whether to use parallel processing (default = FALSE)

- seed:

  Integer, random seed for reproducibility

- observed_data:

  Optional named list keyed by property name (plus, for the texture
  group, `claytotal`/`sandtotal`/`silttotal` supplied together)
  providing a Bayesian-updating likelihood to fuse against each
  property's SSURGO-derived prior before simulating - see
  [`fuse_observed_data_into_priors()`](https://jjmaynard.github.io/soilSIM/reference/fuse_observed_data_into_priors.md)
  for the accepted shapes and behavior. `NULL` (default) skips fusion
  entirely, preserving today's prior-only behavior exactly.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List containing comprehensive simulation results

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic simulation
mc_results <- generate_monte_carlo_realizations(
  soil_data = soil_data,
  properties = c("sandtotal", "claytotal", "silttotal"),
  n_realizations = 1000
)

# Advanced simulation with configuration
mc_results <- generate_monte_carlo_realizations(
  soil_data = soil_data,
  properties = c("sandtotal", "claytotal", "dbovendry"),
  correlation_matrix = cor_matrix,
  n_realizations = 5000,
  simulation_config = list(
    max_depth = 200,
    constraint_rules = "strict",
    distribution_type = "triangular",
    quality_threshold = 0.8
  ),
  parallel = TRUE,
  seed = 12345
)
} # }
```
