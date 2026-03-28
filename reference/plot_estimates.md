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
  model_2_name = "Model 2"
)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  A string specifying the quantity to plot. Options are "incidence",
  "undiag_prev", "undiag_migration", "all_migration", "diag_prob",
  "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and
  "detect_prob".

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

## Value

A ggplot2 object.
