# Summarize Distribution Setup

Real summary of `distributions` (a per-horizon list of per-property
distribution configs from `create_distribution_config()`/
`create_fallback_distribution()`, both of which carry a `valid` and a
`family` field): counts and per-family breakdown.

## Usage

``` r
summarize_distributions(distributions, properties, distribution_type)
```

## Arguments

- distributions:

  Per-horizon list of per-property distribution configs.

- properties:

  Unused (kept for interface compatibility).

- distribution_type:

  Fallback family label for entries with no `family`.

## Value

List with `total_distributions`, `valid_distributions`,
`distribution_type`, `family_counts`.
