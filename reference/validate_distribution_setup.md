# Validate Distribution Setup

Real validity-rate computation over `distributions` (see
[`summarize_distributions()`](https://jjmaynard.github.io/soilSIM/reference/summarize_distributions.md)
for the expected shape).

## Usage

``` r
validate_distribution_setup(distributions, properties, config)
```

## Arguments

- distributions:

  Per-horizon list of per-property distribution configs.

- properties:

  Unused (kept for interface compatibility).

- config:

  Unused (kept for interface compatibility).

## Value

List with `valid_distributions`, `total_distributions`, `validity_rate`.
