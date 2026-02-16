# Prevalence distribution by year of infection

Calculates the prevalence distribution by year of infection, as
estimated in the final quarter of the model.

## Usage

``` r
prev_by_year_infection(stan_data, fit_mat, i, idx)
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

List containing posterior predictive quantities.
