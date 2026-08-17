# Estimate a correlation matrix robustly from data, with grouped and global fallbacks

Ports `code_ref/brdf/property_simulation.R`'s
[`simulate_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/simulate_soil_properties.md)
correlation pattern:
[`Hmisc::rcorr()`](https://rdrr.io/pkg/Hmisc/man/rcorr.html) per group
when the group has enough complete observations, else falling back;
forces symmetry; repairs non-positive-definite matrices via
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html). When
`group_var` is supplied, each qualifying group's empirical matrix is
combined into a single size-weighted average (mod05's simulation
architecture uses one correlation matrix for the whole run, not a
per-stratum matrix, so this folds stratification information in without
requiring a larger architecture change).

## Usage

``` r
estimate_correlation_matrix_robust(
  data,
  group_var = NULL,
  min_group_n = 5,
  global_fallback = NULL,
  group_fallback_matrices = NULL,
  global_fallback_method = "global_fallback"
)
```

## Arguments

- data:

  A data frame of numeric property columns (one row per horizon),
  optionally with a grouping column named by `group_var`.

- group_var:

  Optional column name to stratify by (e.g. `"genhz"`).

- min_group_n:

  Minimum complete-observation count required to trust a group's (or the
  pooled data's) empirical correlation over the fallback.

- global_fallback:

  Optional fallback correlation matrix; defaults to the identity matrix
  over the numeric columns found in `data`.

- group_fallback_matrices:

  Optional named list of correlation matrices, keyed by `group_var`
  value, used in place of dropping a group that fails to meet
  `min_group_n` (today's default behavior when this is `NULL`). A
  substituted group contributes to the same sample-size-weighted
  combination as successfully-estimated groups, with nominal weight
  `n_obs = min_group_n` (these matrices carry no real sample-size
  metadata, so `min_group_n` - the same bar local data must already
  clear to stand alone - is used as a legible, already-plumbed-through
  stand-in; any single successful empirical group already outweighs one
  substituted group, and multiple dominate further). A group with no
  empirical result AND no matching entry here is still dropped, exactly
  as before.

- global_fallback_method:

  Label for the final-tier `method` when `global_fallback` is used
  because no group (empirical or substituted) and no pooled fit
  succeeded. Defaults to `"global_fallback"` (today's exact string);
  callers supplying a non-identity `global_fallback` (e.g. a
  KSSL-derived prior) should override this so the returned `method`
  honestly distinguishes "genuinely uninformative" from "a real borrowed
  prior."

## Value

`list(matrix=, method=, n_obs=)`. `method` is one of
`"empirical_grouped"`, `"kssl_fallback_grouped"`,
`"empirical_grouped_kssl_blended"`, `"empirical_pooled"`, or
`global_fallback_method`'s value (`"global_fallback"` by default). The
two `kssl_*`-named values are generic outcomes of any
`group_fallback_matrices` use, not KSSL-specific logic in this function

- named for the fallback source this function was built to support
  (`kssl-reference-correlations.R`).
