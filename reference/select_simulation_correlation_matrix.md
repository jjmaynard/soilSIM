# Select a Correlation Matrix and Property Set for Single-Cokey Simulation

Chooses a single flat correlation matrix (and its associated property
set) to hand to
[`simulate_monte_carlo()`](https://jjmaynard.github.io/soilSIM/reference/simulate_monte_carlo.md)
for one cokey's data, from either a plain matrix or a genhz-keyed list
of matrices.

## Usage

``` r
select_simulation_correlation_matrix(
  cokey_data,
  correlation_matrices,
  txt_correlation_matrices
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey.

- correlation_matrices:

  A matrix, a genhz-keyed named list of matrices, or `NULL`.

- txt_correlation_matrices:

  Same shape as `correlation_matrices`, for texture properties; used as
  a fallback when `correlation_matrices` is unusable.

## Value

List with `properties` (character vector) and `correlation_matrix` (a
matrix or `NULL`, letting
[`simulate_monte_carlo()`](https://jjmaynard.github.io/soilSIM/reference/simulate_monte_carlo.md)
estimate one empirically).
