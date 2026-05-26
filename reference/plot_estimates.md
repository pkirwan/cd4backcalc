# Plot estimated model quantities

Generate a ggplot2 object plotting estimated quantities over time,
optionally compared to the simulated values.

## Usage

``` r
plot_estimates(
  hiv_list,
  quantity,
  start_year = 0,
  chain_results = FALSE,
  hiv_list_2 = NULL,
  model_1_name = "Model 1",
  model_2_name = "Model 2",
  credible_intervals = TRUE,
  sim_comparison = c("sampled", "generative")
)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  A string specifying the quantity to plot. Options are "incidence",
  "undiag_prev", "undiag_migration", "all_migration", "diag_prob",
  "diagnoses", "migration_prob", "transfers", "transfers_detected",
  "ratio_abroad_uk", and "detect_prob".

- start_year:

  An integer specifying the starting year for the x-axis (default is 0).

- chain_results:

  A logical indicating whether to plot chain-by-chain estimates (default
  is FALSE).

- hiv_list_2:

  An optional second list of model estimates for comparison (default is
  NULL).

- model_1_name:

  A string specifying the name of the first model for the legend
  (default is "Model 1").

- model_2_name:

  A string specifying the name of the second model for the legend
  (default is "Model 2").

- credible_intervals:

  A logical indicating whether to plot 95% credible intervals (default
  is TRUE).

- sim_comparison:

  A string specifying which simulated values to compare against when
  `hiv_list` contains simulated data. Options are `"sampled"` (the
  stochastic realisation for the particular seed, the default) or
  `"generative"` (deterministic expected values).

## Value

A ggplot2 object.
