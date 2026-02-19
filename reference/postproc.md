# Post-processing model fit

This function decides which post-processing routine needs to be run and
passes the relevant data and arguments.

## Usage

``` r
postproc(fit, stan_data, start_yr = 0, nworkers = 1, chain_results = FALSE)
```

## Arguments

- fit:

  Fitted stan model object and stan_data list.

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

## Examples

``` r
sim_diags <- simulate_diagnoses(sim_seed = 123)
stan_data <- model_type_prep(sim_diags, model_seed = 123)
fit <- fit_stan_model(
   stan_data,
   iter_warmup = 10,
   iter_sampling = 10,
   chains = 1,
   refresh = 0,
   messages = FALSE
)
#> Warning: 1 of 1 chains had an E-BFMI less than 0.3.
#> See https://mc-stan.org/misc/warnings for details.
hiv_list <- postproc(fit, stan_data, start_yr = 0, nworkers = 1)
#> Iteration Number =  10 
#> Iteration Number =  10 
```
