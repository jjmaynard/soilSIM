# Impute Rock Fragment Volume (RFV) Values for One Row

Row-level RFV imputation: texture-informed defaults when `rfv_r` is
missing, floor handling for near-zero values, and a simple
proportional-spread `_l/_h` estimate (+/-30%, clamped to `[0.05, 85]`)
for valid representative values.

## Usage

``` r
impute_rfv_values(row)
```

## Arguments

- row:

  A one-row data frame or list with `claytotal_r`/`sandtotal_r`/
  `silttotal_r`/`rfv_r` fields.

## Value

A one-row data frame with `rfv_l`/`rfv_r`/`rfv_h` set.
