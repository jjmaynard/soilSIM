# Shared minimal ggplot theme for diagnostic plots

No ggplot styling convention existed anywhere in modules/ prior to this
(grepped: zero ggplot()/geom\_\*() calls anywhere despite ggplot2 being
loaded in this file and mod06_gp_modeling.R) - this establishes one
rather than improvising ad hoc per plot below.

## Usage

``` r
theme_soil_diagnostics()
```
