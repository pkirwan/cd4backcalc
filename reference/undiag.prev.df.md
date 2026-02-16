# Undiagnosed prevalence dataframe

Calculates undiagnosed prevalence dataframe from prevalence and
infection distributions.

## Usage

``` r
undiag.prev.df(stan_data, prev_dist, inf_dist, start.yr = 1978, annual = FALSE)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- prev_dist:

  Prevalence distribution by year of infection.

- inf_dist:

  Infection distribution by diagnosis.

- start.yr:

  Starting year for the dataframe.

- annual:

  Logical indicating if annual undiagnosed prevalence is required.

## Value

Dataframe of undiagnosed prevalence estimates.
