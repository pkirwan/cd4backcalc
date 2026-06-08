# Detected transfers dataframe

Calculates detected transfers dataframe from prevalence distribution and
detection probabilities.

## Usage

``` r
transfers_detected.df(
  stan_data,
  prev_dist,
  fit_mat,
  start.yr = 1978,
  annual = FALSE
)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- prev_dist:

  Prevalence distribution by year of infection.

- fit_mat:

  Matrix of posterior samples from stan model.

- start.yr:

  Starting year for the dataframe.

- annual:

  Logical indicating if annual detected transfers is required.

## Value

Dataframe of detected transfer estimates.
