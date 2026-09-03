# Validate an Enhanced Distribution Analysis Result

Structural + fit-quality checks on
[`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md)'s
output, in the
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
/
[`validate_distribution_fidelity()`](https://jjmaynard.github.io/soilSIM/reference/validate_distribution_fidelity.md)
per-item-`tryCatch` style. Flags: a property with no successful family
fit; a nonzero `convergence` code; non-finite
`aic`/`bic`/`loglik`/`param_sd`. Never errors - returns collected
`warnings`. Also tolerates the `_safe` chain's flatter
`fitted_distributions[[prop]]` shape (`mean`/`sd`/...): with no
`aic`/`convergence` keys those checks are simply skipped.

## Usage

``` r
validate_distribution_analysis(distribution_analysis, config)
```

## Arguments

- distribution_analysis:

  Result of
  [`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md)
  (or `_safe`).

- config:

  Analysis configuration (unused today; kept for signature symmetry).

## Value

`list(warnings = <character>)`.
