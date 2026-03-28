# Box plots of mean squared error (MSE).

Box plots of mean squared error (MSE).

## Usage

``` r
mse_plot(
  hiv_list,
  hiv_list_2 = NULL,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  model_1_name = "Model 1",
  model_2_name = "Model 2"
)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- hiv_list_2:

  A second list returned from the `postproc` function for comparison
  (default is NULL).

- quantity:

  The quantity to plot MSE for (default is "incidence", other options
  are "undiag_migration", "ratio_abroad_uk", and "mp").

- offset:

  The offset for text labels.

- lab:

  A logical indicating whether to add labels.

- model_1_name:

  The name of the first model for labeling (default is "Model 1").

- model_2_name:

  The name of the second model for labeling (default is "Model 2").

## Value

A ggplot object of the MSE box plot.
