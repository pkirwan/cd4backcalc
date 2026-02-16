# Detection probabilities dataframe

Calculates detection probabilities dataframe from posterior samples.

## Usage

``` r
detect.prob.df(stan_data, fit_mat, start.yr = 1978)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- fit_mat:

  Matrix of posterior samples from stan model.

- start.yr:

  Starting year for the dataframe.

## Value

Dataframe of detection probabilities estimates.
