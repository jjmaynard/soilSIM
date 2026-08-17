# Render an Arbitrary Value as Indented Markdown Bullet Lines

Recursively renders lists as nested bullets and atomic vectors as a
single formatted bullet; non-atomic, non-list values (e.g. a ggplot
object) are rendered by class name only, deliberately avoiding
[`print()`](https://rdrr.io/r/base/print.html)/[`format()`](https://rdrr.io/r/base/format.html)
on arbitrary objects (which could have side effects).

## Usage

``` r
render_value_as_lines(x, indent = 0)
```

## Arguments

- x:

  Value to render.

- indent:

  Current indentation depth.

## Value

Character vector of markdown lines.
