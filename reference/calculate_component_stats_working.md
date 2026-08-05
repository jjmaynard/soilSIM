# Calculate Basic Component Statistics

Adds `is_major_component` (`comppct_r >= 15`) and, if `taxclname` is
present, `has_taxonomic_classification` flag columns.

## Usage

``` r
calculate_component_stats_working(df, verbose = FALSE)
```

## Arguments

- df:

  Component data frame.

- verbose:

  Logical; log a completion message.

## Value

`df` with the derived columns added.
