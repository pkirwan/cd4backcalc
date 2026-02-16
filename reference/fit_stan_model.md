# Run CD4 back-calculation model

This function runs the CD4 back-calculation model using Stan.

## Usage

``` r
fit_stan_model(
  stan_data,
  model_seed = NULL,
  iter = 2000,
  chains = 4,
  refresh = 100,
  messages = TRUE
)
```

## Arguments

- stan_data:

  A list of data formatted for Stan.

- model_seed:

  Seed for model fitting (will be set to
  stan_data\$model_config\$model_seed if NULL).

- iter:

  Number of iterations for Stan sampling.

- chains:

  Number of chains for Stan sampling.

- refresh:

  Number of iterations between progress messages.

- messages:

  Logical indicating whether to show messages during sampling.

## Value

A stanfit object containing the results of the model fitting.

## Examples

``` r
sim_diags <- simulate_diagnoses(sim_seed = 123)
stan_data <- model_type_prep(sim_diags)
fit <- fit_stan_model(stan_data, iter = 10, chains = 1, refresh = 0, messages = FALSE)
#> Warning: There were 5 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
```
