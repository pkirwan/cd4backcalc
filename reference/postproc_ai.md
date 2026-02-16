# Post-process age-independent model fit

Post-process age-independent model fit

## Usage

``` r
postproc_ai(
  fit,
  stan_data,
  start_yr = 1995,
  nworkers = 1,
  chain_results = FALSE
)
```

## Arguments

- fit:

  Fitted stan model object.

- stan_data:

  Data list used in stan model fitting.

- start_yr:

  Starting year for posterior predictive quantities.

- nworkers:

  Number of parallel workers to use.

- chain_results:

  Logical indicating whether to return results by chain.

## Value

List of posterior predictive quantities.
