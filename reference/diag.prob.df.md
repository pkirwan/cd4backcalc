# Diagnosis probabilities dataframe

Calculates diagnosis probabilities dataframe from posterior samples.

## Usage

``` r
diag.prob.df(stan_data, fit_mat, start.yr = 1978, group = c("nat", "mig"))
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- fit_mat:

  Matrix of posterior samples from stan model.

- start.yr:

  Starting year for the dataframe.

- group:

  Character indicating which diagnosis probabilities to extract (`"nat"`
  for UK-born, `"mig"` for born abroad).

## Value

Dataframe of diagnosis probabilities estimates.
