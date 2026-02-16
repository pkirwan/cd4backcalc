# Box plots of bias.

Box plots of bias.

## Usage

``` r
bias_plot(hiv_list, quantity = "incidence", offset = 0.1, lab = TRUE)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  The quantity to plot bias for (default is "incidence", other options
  are "undiag_migration", "ratio_abroad_uk", and "mp").

- offset:

  The offset for text labels.

- lab:

  A logical indicating whether to add labels.

## Value

A ggplot object of the bias box plot.
