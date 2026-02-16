# Migration dataframe

Calculates migration dataframe from prevalence distribution.

## Usage

``` r
total.migration.df(stan_data, prev_dist, start.yr = 1978, annual = FALSE)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- prev_dist:

  Prevalence distribution by year of infection.

- start.yr:

  Starting year for the dataframe.

- annual:

  Logical indicating if annual migration is required.

## Value

Dataframe of migration estimates.
