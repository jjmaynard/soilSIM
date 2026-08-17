# Assess Distribution Analysis Quality

Real score based on the fraction of per-property distribution-fit
entries that produced a non-empty `fits` list.

## Usage

``` r
assess_distribution_quality(distribution_analysis)
```

## Arguments

- distribution_analysis:

  Named list keyed by property, each entry optionally with a `fits`
  field (see
  [`fit_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/fit_property_distributions.md)-style
  output).

## Value

List with `quality_score`, `n_properties`.
