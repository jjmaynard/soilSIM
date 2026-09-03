# Probabilities behind a `run_stage1_fusion()` posterior's percentile layers

Probabilities behind a
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
posterior's percentile layers

## Usage

``` r
posterior_percentile_probs(posterior)
```

## Arguments

- posterior:

  A
  [`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)-style
  `list(percentiles = <named P.. rasters>, ...)`.

## Value

Numeric probabilities in (0, 1), ascending, parsed from the `"P<pct>"`
layer names.
