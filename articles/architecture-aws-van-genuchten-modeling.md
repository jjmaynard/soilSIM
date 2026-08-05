# AWS / Van Genuchten Modeling

## Overview

This functional area covers available water storage (AWS) modeling:
evaluating the closed-form van Genuchten (1980) water-retention curve,
Monte Carlo-simulating available water holding capacity (AWHC) from
ROSETTA-derived pedotransfer parameters, and depth-slicing/summarizing
AWHC per soil component. It is implemented entirely in
`R/aws-simulation.R` and is self-contained - unlike
`R/depth-simulation.R`’s functions, none of the three functions here
depend on any other not-yet-ported helper
(e.g. [`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)).
It can be read, tested, and used in isolation from the rest of the
`soilSIM` migration.

## Core Functions

### 1. `van_genuchten()` - Water Retention Curve Evaluator

**Purpose**: Evaluates the closed-form van Genuchten (1980) equation for
volumetric water content at a given matric potential.

**Parameters**:

``` r

van_genuchten(h, alpha, n, theta_r, theta_s)
```

- `h` - Matric potential (cmH2O, typically negative).
- `alpha`, `n` - Van Genuchten shape parameters.
- `theta_r`, `theta_s` - Residual and saturated volumetric water
  content.

**Returns**: A numeric vector of volumetric water content, `theta(h)`.

**Algorithm/behavior**: The function computes `m <- 1 - (1 / n)` and
then evaluates the standard closed-form van Genuchten retention
equation:

    theta(h) = theta_r + (theta_s - theta_r) / (1 + |alpha * h|^n)^m

wrapped in [`as.numeric()`](https://rdrr.io/r/base/numeric.html) to
strip any attributes carried by vectorized inputs. It is a pure,
vectorized, side-effect-free function with no dependence on external
packages - just base R arithmetic - which is why it is reused directly
inside
[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)’s
row-wise Monte Carlo loop.

### 2. `simulate_vg_aws()` - Monte Carlo AWHC Simulation

**Purpose**: Monte Carlo-simulates field capacity (FC) and permanent
wilting point (PWP) volumetric water content - and their difference,
available water holding capacity (AWHC) - for each row (soil
layer/component) of van Genuchten parameters.

**Parameters**:

``` r

simulate_vg_aws(data, n_simulations = 100)
```

- `data` - A data frame with one row per soil layer/component, and
  columns `alpha`, `sd_alpha`, `npar`, `sd_npar`, `theta_r`,
  `sd_theta_r`, `theta_s`, `sd_theta_s`, `layerID` - exactly the shape
  produced by `soilDB::ROSETTA(..., include.sd = TRUE)` plus a
  caller-added `layerID` (see
  [`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)).
- `n_simulations` - Number of Monte Carlo draws per row (default `100`).

**Returns**: A named list (one element per row, keyed
`paste(layerID, i, sep = "_")`), each element a data frame of
`n_simulations` draws with columns `alpha`, `n`, `theta_r`, `theta_s`,
`sim_num`, `theta_fc`, `theta_pwp`, `AWHC`. Rows with any missing van
Genuchten parameter are skipped entirely (absent from the result list,
not present as `NA` rows).

**Algorithm/behavior**: For each row `i` of `data`, the function
extracts the mean/SD pair for each of `alpha`, `npar`, `theta_r`,
`theta_s`. If any of the four means is `NA`, the row is skipped via
[`next`](https://rdrr.io/r/base/Control.html). Otherwise it calls
`set.seed(123)` and draws `n_simulations` samples from
`stats::rnorm(mean, sd)` for each of the four parameters. The sampled
`alpha` and `n` values are back-transformed with `10^(...)` (treating
the sampled values as log10-space draws), while `theta_r`/`theta_s`
samples are used as-is. Fixed matric potentials for field capacity and
permanent wilting point are converted from kPa to cmH2O (1 kPa =
10.19716 cmH2O): `h_fc = -33 * 10.19716` and `h_pwp = -1500 * 10.19716`.
A per-draw data frame is built with the sampled/back-transformed
parameters plus `sim_num`, then
[`dplyr::rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html) +
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
calls
[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)
once at `h_fc` and once at `h_pwp` to get `theta_fc` and `theta_pwp` for
each draw; `AWHC <- theta_fc - theta_pwp` is computed as a plain
vectorized column after ungrouping. The resulting per-row data frame is
stored in the output list under key
`paste(data$layerID[i], i, sep = "_")`.

**Ported-as-is behavioral quirks** (documented in the source as
intentional, not bugs): - `set.seed(123)` is called **inside** the
per-row loop, so every row re-seeds and draws from an identically-seeded
random stream rather than an evolving stream across rows. This means the
`n_simulations` draws for row 1 and row 2 (before back-transformation)
are numerically identical in their underlying
[`rnorm()`](https://rdrr.io/r/stats/Normal.html) output, differing only
through each row’s own mean/SD. - Sampled `alpha`/`n` values are
back-transformed via `10^(...)`, i.e. the input `data$alpha`/`data$npar`
(and their SDs) are treated as already being in log10 space. This is a
standard technique for keeping van Genuchten shape parameters positive
after adding Gaussian noise, but whether
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)’s
actual reported values are meant to be interpreted this way is not
independently verified by this port - it is preserved exactly as it
existed in the legacy source.

### 3. `.aws_slab_mean()` - Internal `slab.fun` Helper (not exported)

**Purpose**: A minimal replacement for the no-longer-available
`aqp::mean_na()` (removed/renamed in current `aqp` versions), preserving
its old single-value-per-slab contract that
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)’s
current default `slab.fun` (`slab_function(method = "numeric")`) does
not provide (the current default returns quantile columns instead of a
single `value` column).

**Parameters**:

``` r

.aws_slab_mean(values, ...)
```

- `values` - Numeric vector of observations within one depth slab.
- `...` - Ignored; present because
  [`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)
  passes additional arguments positionally.

**Returns**: A single numeric value, `mean(values, na.rm = TRUE)`.

**Algorithm/behavior**: One-line pass-through to
`mean(values, na.rm = TRUE)`. It exists purely so
[`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)’s
call to
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)
yields a `value` column that
[`tidyr::pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
can consume, matching the shape the legacy code relied on from the
now-removed `aqp::mean_na()`.

### 4. `calculate_aws_df()` - Master AWS-by-Depth-Interval Function

**Purpose**: Runs
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)
on soil texture/bulk-density/water-retention inputs to derive van
Genuchten pedotransfer parameters, Monte Carlo-simulates AWHC per
horizon via
[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md),
then depth-slices and summarizes mean AWHC per component over standard
depth intervals via
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html).

**Parameters**:

``` r

calculate_aws_df(sim_data_df)
```

- `sim_data_df` - A data frame with one row per horizon, with columns
  `sand_total`, `silt_total`, `clay_total`, `bulk_density_third_bar`,
  `water_retention_third_bar`, `water_retention_15_bar` (ROSETTA’s
  expected input variable names - note these differ from the rest of
  `soilSIM`’s `sandtotal_r`/`claytotal_r`/etc. SSURGO-derived naming
  convention, since they’re ROSETTA’s own API contract), plus
  `compname`, `hzdept_r`, `hzdepb_r`, `cokey`.

**Returns**: A long-format data frame with one row per `cokey` per depth
slab actually spanned by that component’s horizons (slab boundaries per
`slab.structure = c(0, 5, 15, 30, 60, 100)`), with columns `cokey`,
`top`, `bottom`, `AWHC` (mean available water holding capacity over that
slab). The
[`tidyr::pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
step widens on `variable` (always `"AWHC"` here, since only one property
is slabbed), so it does not collapse depth slabs into columns - it
exists only to match
[`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)’s
long output shape to a plain `value`-column contract.

**Algorithm/behavior**: The function first guards on `httr` being
installed
([`requireNamespace("httr", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html)),
since
[`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)
uses it internally to make the live API call; if missing, it stops with
an explicit message. It then calls
`soilDB::ROSETTA(sim_data_df, vars = variables, v = "3", include.sd = TRUE)`
where `variables` is the fixed vector of the six ROSETTA input names,
requesting ROSETTA model version `"3"` with standard deviations
included. A `layerID` column is constructed as
`paste(compname, hzdept_r, sep = "_")` on `sim_data_df` and copied onto
`rosetta_data` so
[`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)
can key its result list. `simulate_vg_aws(rosetta_data)` is then called
to Monte Carlo-simulate AWHC per horizon. The resulting list of
per-horizon simulation data frames is passed through
[`Map()`](https://rdrr.io/r/base/funprog.html) alongside `sim_data_df`’s
`hzdept_r`, `hzdepb_r`, `compname`, and `cokey` vectors, attaching
`top`, `bottom`, `compname`, and `cokey` columns to each element; the
list is then flattened with
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
into `sim_aws_df`. `aqp::depths(sim_aws_df) <- cokey ~ top + bottom`
promotes the data frame to an `aqp` SoilProfileCollection keyed by
`cokey` with depths taken from the newly-attached `top`/`bottom` columns
(not the original `hzdept_r`/`hzdepb_r` names, which no longer exist on
this object).
`aqp::slab(sim_aws_df, fm = cokey ~ AWHC, slab.structure = c(0, 5, 15, 30, 60, 100), slab.fun = .aws_slab_mean)`
then computes the mean AWHC within each of the fixed depth slabs (0-5,
5-15, 15-30, 30-60, 60-100 cm) per `cokey`, using the local
`.aws_slab_mean()` in place of the removed `aqp::mean_na()`. Finally,
the slab output has its `contributing_fraction` column dropped and is
pivoted wider via
`tidyr::pivot_wider(names_from = variable, values_from = value)`, and
any remaining `NA` values are reassigned to `NA` (a no-op cleanup line
preserved from the original source) before the data frame is returned.

## Internal Connections

    calculate_aws_df(sim_data_df)
    ├── soilDB::ROSETTA(sim_data_df, vars, v = "3", include.sd = TRUE)   [live network call to handbook60.org]
    │   └── (internally uses httr::POST())
    ├── layerID construction (paste(compname, hzdept_r, sep = "_"))
    ├── simulate_vg_aws(rosetta_data, n_simulations = 100)
    │   ├── stats::rnorm() draws per row (set.seed(123) each iteration)
    │   └── van_genuchten(h_fc, alpha, n, theta_r, theta_s)
    │       van_genuchten(h_pwp, alpha, n, theta_r, theta_s)
    ├── Map() attaches top/bottom/compname/cokey to each simulated data frame
    ├── dplyr::bind_rows() flattens the list into sim_aws_df
    ├── aqp::depths(sim_aws_df) <- cokey ~ top + bottom
    ├── aqp::slab(sim_aws_df, fm = cokey ~ AWHC, slab.structure = ..., slab.fun = .aws_slab_mean)
    │   └── .aws_slab_mean(values, ...) -> mean(values, na.rm = TRUE)
    └── tidyr::pivot_wider(names_from = variable, values_from = value) -> returned data frame

## Dependencies

### External packages

    soilSIM (aws-simulation.R)
    ├── soilDB::ROSETTA()   - REQUIRED, live network call (POST to handbook60.org); derives van Genuchten
    │                          pedotransfer parameters (alpha, npar, theta_r, theta_s + SDs) from
    │                          texture/bulk-density/water-retention inputs
    ├── aqp::depths<-()     - promotes a data frame to a SoilProfileCollection keyed by cokey ~ top + bottom
    ├── aqp::slab()         - depth-slice summarization across fixed depth intervals
    ├── stats::rnorm()      - Monte Carlo parameter sampling
    ├── dplyr (rowwise, mutate, ungroup, bind_rows) - simulation and combination pipeline
    ├── tidyr::pivot_wider()- reshapes aqp::slab()'s long output to a value-per-row shape
    └── httr                - Suggests-guarded: not called directly by this file, but required
                               transitively because soilDB::ROSETTA() uses it internally for the
                               HTTP POST; calculate_aws_df() checks requireNamespace("httr", ...)
                               up front and stops with a clear message if it's absent

### soilSIM dependencies / consumers

This module is a leaf/standalone module within `soilSIM`. It does not
call into any other not-yet-ported `soilSIM` helper (unlike
`R/depth-simulation.R`, which depends on
[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)).
Upstream, component-level texture, bulk-density, and water-retention
data produced by SSURGO acquisition/property-simulation steps elsewhere
in `soilSIM` can be reshaped into the
`sand_total`/`silt_total`/`clay_total`/`bulk_density_third_bar`/`water_retention_third_bar`/`water_retention_15_bar`/`compname`/`hzdept_r`/`hzdepb_r`/`cokey`
shape
[`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)
expects, but no such wiring exists inside this file itself - the caller
is responsible for producing `sim_data_df` in the expected shape.

## Data Flow In/Out

**In**: A per-horizon data frame (`sim_data_df`) of component texture
and bulk-density/water-retention data keyed by `cokey`, `compname`,
`hzdept_r`, `hzdepb_r`, with the six ROSETTA-named property columns
(`sand_total`, `silt_total`, `clay_total`, `bulk_density_third_bar`,
`water_retention_third_bar`, `water_retention_15_bar`).

**Out**: A long-format data frame with one row per `cokey` per depth
slab (bounded by `slab.structure = c(0, 5, 15, 30, 60, 100)`), columns
`cokey`, `top`, `bottom`, `AWHC` - mean available water holding capacity
for that component over that depth slab. This shape is suitable for
direct join/merge back onto other per-`cokey` or per-depth-slab soil
property tables downstream.

## Known Limitations

- **Live network dependency**:
  [`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)
  requires live network access.
  [`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)
  POSTs to `https://www.handbook60.org/api/v1/rosetta/<version>` via
  [`httr::POST()`](https://httr.r-lib.org/reference/POST.html) rather
  than computing pedotransfer parameters locally. Calling this function
  without internet connectivity, or if `handbook60.org` is unreachable
  or down, will fail. There is no offline/local fallback path in this
  file.
- **Per-row `set.seed(123)` re-seeding**
  ([`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)):
  the random seed is reset inside the per-row loop rather than being set
  once before the loop, so each row’s Monte Carlo draws come from an
  identically-seeded random stream. This is preserved intentionally from
  the legacy source rather than “fixed,” since it may be a deliberate
  reproducibility choice and changing it would alter the simulated AWHC
  values every caller receives.
- **`10^(...)` back-transformation of `alpha`/`npar`**
  ([`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)):
  sampled `alpha` and `n` values are exponentiated as if the input
  means/SDs were supplied in log10 space. This is a standard technique
  for keeping van Genuchten shape parameters positive, but whether
  [`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)’s
  actual output values are intended to be interpreted this way has not
  been independently verified by this port - it is preserved as-is.
- **`.aws_slab_mean()` as a replacement for `aqp::mean_na()`**: current
  versions of `aqp` no longer export `mean_na()` (removed/renamed), and
  [`aqp::slab()`](https://ncss-tech.github.io/aqp/reference/slab.html)’s
  current default `slab.fun` returns quantile columns rather than a
  single `value` column. `.aws_slab_mean()` is a minimal local stand-in
  restoring the old single-value-per-slab contract; if `aqp`’s API
  changes further, this shim may need revisiting.
- **Silent row-skipping in
  [`simulate_vg_aws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_vg_aws.md)**:
  rows with any missing van Genuchten parameter mean are silently
  omitted from the returned list (no warning), so a horizon that ROSETTA
  failed to return usable parameters for simply disappears from
  downstream results rather than surfacing as an explicit `NA` row.

## Usage Example

``` r

library(soilSIM)

# sim_data_df: one row per horizon, in ROSETTA's expected input shape,
# plus compname/hzdept_r/hzdepb_r/cokey identifiers.
sim_data_df <- data.frame(
  cokey                      = c("12345678", "12345678"),
  compname                   = c("ExampleSoil", "ExampleSoil"),
  hzdept_r                   = c(0, 25),
  hzdepb_r                   = c(25, 60),
  sand_total                 = c(35, 30),
  silt_total                 = c(45, 40),
  clay_total                 = c(20, 30),
  bulk_density_third_bar     = c(1.35, 1.40),
  water_retention_third_bar  = c(0.28, 0.30),
  water_retention_15_bar     = c(0.12, 0.15)
)

# NOTE: requires live network access - calls soilDB::ROSETTA(), which POSTs
# to https://www.handbook60.org via httr::POST(). Will fail if offline or
# if handbook60.org is unreachable.
aws_by_depth <- calculate_aws_df(sim_data_df)

# aws_by_depth: long-format data frame with columns cokey, top, bottom, AWHC
# - one row per depth slab (0-5, 5-15, 15-30, 30-60, 60-100 cm) actually
#   spanned by the component's horizons.
head(aws_by_depth)
```
