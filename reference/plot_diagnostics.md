# Plot diagnostics for fitted model

Generate diagnostic plots for a fitted back-calculation model.

## Usage

``` r
plot_diagnostics(hiv_list, ntraces = 4)
```

## Arguments

- hiv_list:

  A list returned from the `run_backcalc` function containing model
  estimates and fit objects.

- ntraces:

  Number of parameter traces to plot for each diagnostic.

## Value

A ggplot2 object containing diagnostic plots.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
