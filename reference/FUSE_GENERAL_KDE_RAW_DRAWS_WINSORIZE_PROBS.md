# Default `winsorize_probs` `fuse_general_kde()` passes to `update_prior()` on its `raw_draws` branch only.

Default `winsorize_probs`
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
passes to
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
on its `raw_draws` branch only.

## Usage

``` r
FUSE_GENERAL_KDE_RAW_DRAWS_WINSORIZE_PROBS
```

## Format

An object of class `numeric` of length 2.

## Why this exists

[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
`stats::density(bw = "nrd0")` bandwidth choice is sensitive to outliers
a raw-draws resample can genuinely contain - confirmed via a dedicated
adversarial test: a right-skewed synthetic prior's raw resample
(realized range ~4-195) produced a measurably wider KDE bandwidth, and a
*less* accurate fused posterior mean on average across 8 seeds, than the
percentile-reconstruction route's resample (confined to ~7-64, the same
data's P5-P95 range). Clipping each side to its own P01-P99 range before
density estimation directly counters that outlier-driven inflation while
retaining the bulk of the real distribution's shape - only the most
extreme 1% on each tail is affected, which is exactly the region a
KDE-based fused mean is most vulnerable to.

[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
own standalone default is deliberately left `NULL` (unclipped) - this
constant only affects
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
`raw_draws` branch, which is the one path confirmed to need it; the
default percentile-reconstruction branch is already structurally bounded
near its own P5-P95 input range and is left untouched.
