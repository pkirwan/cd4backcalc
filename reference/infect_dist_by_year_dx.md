# Infection distribution by year of diagnosis

Calculates the distribution of infections by year of diagnosis.

## Usage

``` r
infect_dist_by_year_dx(stan_data, fit_mat, i, idx)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- fit_mat:

  Matrix of posterior samples from stan model.

- i:

  Index of the posterior sample to process.

- idx:

  Additional index parameter for processing.

## Value

List containing HIV and AIDS diagnosis infection distributions.
