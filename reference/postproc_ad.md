# Post-process the age-dependent model fit

Post-processes the fitted age-dependent back-calculation model to
extract various summaries and diagnostics.

## Usage

``` r
postproc_ad(fit, stan_data, start_yr = 0, nworkers = 1, chain_results = FALSE)
```

## Arguments

- fit:

  Fitted stan model object.

- stan_data:

  Data list used in stan model fitting.

- start_yr:

  Starting year for posterior predictive quantities.

- nworkers:

  Number of workers for parallel processing.

- chain_results:

  Logical indicating if chain-specific results are required.

## Value

List of post-processed results including summaries and diagnostics.
