# Generate Markdown Validation Report

Writes `report_content` to `output_file` as a plain markdown file
(dependency-free).

## Usage

``` r
generate_markdown_report(report_content, output_file)
```

## Arguments

- report_content:

  Report content from `create_report_content()`.

- output_file:

  Output file path (`.md`).

## Value

`TRUE` on success, `FALSE` on failure (logged as a warning).
