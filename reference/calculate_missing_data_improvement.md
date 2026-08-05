# Calculate Missing-Data Handling Improvement

Real reduction in missing-value rate (averaged across shared numeric
columns) between `original_data` and `processed_data`.

## Usage

``` r
calculate_missing_data_improvement(original_data, processed_data)
```

## Arguments

- original_data:

  Original input data.

- processed_data:

  Data after infilling/processing.

## Value

Numeric improvement in `[0, 1]` (0 = no improvement), or `NA` if not
computable.
