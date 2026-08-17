# Fuse one property's per-horizon prior against a shared observed-data likelihood

Fuse one property's per-horizon prior against a shared observed-data
likelihood

## Usage

``` r
fuse_one_property_prior(
  prior,
  likelihood,
  is_vector_likelihood,
  n_samples,
  supported_closed_form,
  prop,
  horizon_index
)
```

## Arguments

- prior:

  A `simulation_params[[i]][[prop]]` entry:
  `list(family=, fit=, source=)`.

- likelihood:

  The `observed_data[[prop]]` entry (raw vector or param list).

- is_vector_likelihood:

  Precomputed `TRUE`/`FALSE` for `likelihood`'s shape.

- n_samples:

  Sample size for the general route.

- supported_closed_form:

  Character vector of families with a conjugate route.

- prop, horizon_index:

  Used only for WARN message context.

## Value

A new `list(family=, fit=, source=)` posterior, or `NULL` to skip
(leaving the prior unchanged).
