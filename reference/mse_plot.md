# Box plots of mean squared error (MSE).

Box plots of mean squared error (MSE).

## Usage

``` r
mse_plot(hiv_list, quantity = "incidence", offset = 0.1, lab = TRUE)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  The quantity to plot MSE for (default is "incidence", other options
  are "undiag_migration", "ratio_abroad_uk", and "mp").

- offset:

  The offset for text labels.

- lab:

  A logical indicating whether to add labels.

## Value

A ggplot object of the MSE box plot.
