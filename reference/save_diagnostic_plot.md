# Save a ggplot object to output_dir if one was provided

Save a ggplot object to output_dir if one was provided

## Usage

``` r
save_diagnostic_plot(plot_obj, filename, output_dir)
```

## Arguments

- plot_obj:

  A ggplot object to save.

- filename:

  File name (not path) to save the plot under.

- output_dir:

  Directory to save into, or `NULL` to skip saving.

## Value

The ggplot object itself (output_dir NULL or save failed) or the saved
file path (output_dir provided and save succeeded).
