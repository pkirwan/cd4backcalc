# Box plots of bias.

Box plots of bias.

## Usage

``` r
bias_plot(
  hiv_list,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  hiv_list_2 = NULL,
  model_1_name = "Model 1",
  model_2_name = "Model 2",
  sim_comparison = c("sampled", "generative")
)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  The quantity to plot bias for (default is `"incidence"`).
  Migration-model options are `"undiag_migration"`, `"ratio_abroad_uk"`,
  and `"mp"`.

- offset:

  The offset for text labels.

- lab:

  A logical indicating whether to add labels.

- hiv_list_2:

  A second list returned from the `postproc` function for comparison
  (default is NULL).

- model_1_name:

  The name of the first model for labeling (default is "Model 1").

- model_2_name:

  The name of the second model for labeling (default is "Model 2").

- sim_comparison:

  A string specifying which simulated values to compare against. Options
  are `"sampled"` (the stochastic realisation for the particular seed,
  the default) or `"generative"` (deterministic expected values).

## Value

A ggplot object of the bias box plot.
