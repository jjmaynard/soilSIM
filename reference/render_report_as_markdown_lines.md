# Render Report Content as Markdown Lines

Dependency-free renderer that turns the nested list produced by
`create_report_content()` into plain markdown lines. Shared by
[`generate_markdown_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_markdown_report.md),
[`generate_html_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_html_report.md),
and
[`generate_pdf_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_pdf_report.md)'s
dependency-free fallback paths.

## Usage

``` r
render_report_as_markdown_lines(report_content)
```

## Arguments

- report_content:

  Nested list report content.

## Value

Character vector of markdown lines.
