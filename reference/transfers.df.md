# Transfers dataframe

Calculates transfers dataframe from prevalence distribution.

## Usage

``` r
transfers.df(stan_data, prev_dist, start.yr = 1978, annual = FALSE)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- prev_dist:

  Prevalence distribution by year of infection.

- start.yr:

  Starting year for the dataframe.

- annual:

  Logical indicating if annual transfers is required.

## Value

Dataframe of transfer estimates.
