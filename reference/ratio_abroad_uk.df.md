# UK-acquired incidence-ratio dataframe

Calculates the born-abroad:UK-born ratio of UK-acquired incidence from
posterior samples.

## Usage

``` r
ratio_abroad_uk.df(stan_data, fit_mat, start.yr = 1978)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- fit_mat:

  Matrix of posterior samples from stan model.

- start.yr:

  Starting year for the dataframe.

## Value

Data frame of born-abroad:UK-born UK-acquired incidence-ratio estimates.
