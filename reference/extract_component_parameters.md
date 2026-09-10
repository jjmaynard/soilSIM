# Extract Component Distribution Parameters

Extracts `min`/`mode`/`max` for one component row from its
`comppct_l`/`comppct_r`/`comppct_h` triplet, filling missing
`comppct_l`/`comppct_h` from `comppct_r -/+ 2`.

## Usage

``` r
extract_component_parameters(component_row, config)
```

## Arguments

- component_row:

  One-row data frame/list for a single component.

- config:

  Unused (kept for interface compatibility).

## Value

`list(min=, mode=, max=)`.
