# Detect Distribution Anomalies

Runs IQR-based outlier detection (reusing the already-real
[`detect_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_outliers.md))
across `simulation_data`'s numeric property columns.

## Usage

``` r
detect_distribution_anomalies(simulation_data, distribution_checks)
```

## Arguments

- simulation_data:

  Simulation data.

- distribution_checks:

  Unused (no per-check configuration is currently defined upstream -
  reserved for future per-property checks).

## Value

List with `anomalies_detected` (count), `anomaly_types`,
`anomaly_severity`.
