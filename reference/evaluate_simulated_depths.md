# Evaluate Simulated Profile Depths

Evaluates whether the simulated soil profile depths fall outside the
specified range. Extracts the simulated top and bottom depths from a
`SoilProfileCollection`, merges these with the original horizon data,
and flags horizons where the simulated top is less than the lower bound
or the simulated bottom exceeds the upper bound.

## Usage

``` r
evaluate_simulated_depths(simulated_profiles, horizon_data)
```

## Arguments

- simulated_profiles:

  A SoilProfileCollection object containing simulated soil profiles.

- horizon_data:

  A data frame with original horizon data, including columns: -
  `hzname`: Horizon name. - `hzdept_l`: Lower bound for the top depth. -
  `hzdepb_h`: Upper bound for the bottom depth.

## Value

A data frame containing rows (horizons) where the simulated depths are
out of range.

## Examples

``` r
if (FALSE) { # \dontrun{
  out_of_range <- evaluate_simulated_depths(simulated_profiles, horizon_data)
  head(out_of_range)
} # }
```
