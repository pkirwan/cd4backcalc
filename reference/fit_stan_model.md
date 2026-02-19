# Run CD4 back-calculation model

This function runs the CD4 back-calculation model using Stan.

## Usage

``` r
fit_stan_model(
  stan_data,
  model_seed = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
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

- iter_warmup:

  Number of warmup iterations for Stan sampling.

- iter_sampling:

  Number of sampling iterations for Stan sampling.

- chains:

  Number of chains for Stan sampling.

- adapt_delta:

  Adapt delta parameter for Stan sampling (default is 0.95).

- max_treedepth:

  Max tree depth for Stan sampling (default is 12).

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
fit <- fit_stan_model(
 stan_data, iter_warmup = 10, iter_sampling = 10, chains = 1, refresh = 0, messages = FALSE
)
#> Warning: 9 of 10 (90.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
```
