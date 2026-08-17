# Generate HTML Validation Report

Writes `report_content` to `output_file` as HTML. Uses
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
for richer output when `rmarkdown` is installed and `pandoc` is
available (both optional, `Suggests`-only dependencies); otherwise falls
back to a dependency-free `<pre>`-wrapped rendering that always works.

## Usage

``` r
generate_html_report(report_content, output_file)
```

## Arguments

- report_content:

  Report content from `create_report_content()`.

- output_file:

  Output file path (`.html`).

## Value

`TRUE` on success, `FALSE` on failure (logged as a warning).
