# Compute the Mode of a Numeric Vector

Compute the Mode of a Numeric Vector

## Usage

``` r
calculate_mode(x)
```

## Arguments

- x:

  A vector of values.

## Value

The most frequently-occurring value in `x`.

## Performance

`table(x)`/`factor(x)` build a full factor/hash table just to count
occurrences - [`Rprof()`](https://rdrr.io/r/utils/Rprof.html) profiling
on
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
(its only caller, invoked twice per horizon row) attributed 22% of that
function's total wall-clock to this call alone
(PERFORMANCE_IMPROVEMENT_PLAN.md Tier 4). `tabulate(match(x, ux))` on
pre-sorted unique values computes the same counts without the
factor-coercion overhead. Verified to return identical output (including
[`table()`](https://rdrr.io/r/base/table.html)'s implicit tie-break -
the smallest value among ties, since
[`table()`](https://rdrr.io/r/base/table.html) sorts unique values
ascending and [`which.max()`](https://rdrr.io/r/base/which.min.html)
returns the first maximum) across ties, negative values, singletons, and
random floating-point draws.
