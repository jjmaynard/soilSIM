# Generate PDF Validation Report

Writes `report_content` to `output_file` as a PDF. Uses
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
for richer output when `rmarkdown`/`tinytex` are installed, `pandoc` is
available, and a working TeX distribution is detected (all optional,
`Suggests`-only dependencies); otherwise falls back to a dependency-free
plain-text rendering via
[`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html) that always
works.

## Usage

``` r
generate_pdf_report(report_content, output_file)
```

## Arguments

- report_content:

  Report content from `create_report_content()`.

- output_file:

  Output file path (`.pdf`).

## Value

`TRUE` on success, `FALSE` on failure (logged as a warning).
