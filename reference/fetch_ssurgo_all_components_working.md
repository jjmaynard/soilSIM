# Fetch All Components for a Set of Map Unit Keys (No Horizon Join)

Companion query to
[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md).
That function's `INNER JOIN chorizon` silently drops any real SSURGO
component that has zero `chorizon` (horizon-level property) rows in
SDA - a real, verified data-completeness gap, especially common for
minor components (e.g. a 3%-share component whose horizon table was
never populated). This function queries the `component` table alone,
with no `chorizon` join, so callers can detect exactly which components
[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)'s
result is missing (see
[`recover_missing_horizon_components()`](https://jjmaynard.github.io/soilSIM/reference/recover_missing_horizon_components.md)).

## Usage

``` r
fetch_ssurgo_all_components_working(mukey_list, verbose = FALSE)
```

## Arguments

- mukey_list:

  Vector of map unit keys (same input
  [`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md)
  takes).

- verbose:

  Logical; provide progress messages.

## Value

Data frame with one row per component: `mukey`, `cokey`, `compname`,
`comppct_l`, `comppct_r`, `comppct_h`. An empty (0-row) data frame with
these columns if the query returns nothing or fails outright - this is a
best-effort enrichment step, not a hard requirement, so failures are
logged and swallowed rather than raised (unlike
[`execute_ssurgo_query_working()`](https://jjmaynard.github.io/soilSIM/reference/execute_ssurgo_query_working.md),
whose own query failure is fatal to the whole download).
