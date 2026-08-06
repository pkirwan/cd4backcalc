# Box plots of mean squared error (MSE).

Box plots of mean squared error (MSE).

## Usage

``` r
mse_plot(
  hiv_lists,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  colours = NULL,
  sim_comparison = c("sampled", "generative"),
  years = NULL
)
```

## Arguments

- hiv_lists:

  A named list of model result lists, each returned from the `postproc`
  function. Names are used as model labels. If unnamed, models are
  labelled "Model 1", "Model 2", etc.

- quantity:

  The quantity to plot MSE for (default is `"incidence"`).
  Migration-model options are `"undiag_migration"`, `"ratio_abroad_uk"`,
  and `"mp"`.

- offset:

  The offset for text labels.

- lab:

  A logical indicating whether to add labels.

- colours:

  An optional named character vector of colours for each model. Names
  must match the names of `hiv_lists`. If `NULL` (default), ggplot's
  default discrete colour palette is used.

- sim_comparison:

  A string specifying which simulated values to compare against. Options
  are `"sampled"` (the stochastic realisation for the particular seed,
  the default) or `"generative"` (deterministic expected values).

- years:

  Optional numeric/integer vector specifying which years to include in
  bias/MSE calculations. If `NULL` (default), all available years are
  used.

## Value

A ggplot object of the MSE box plot.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
