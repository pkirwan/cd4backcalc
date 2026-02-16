# Run CD4 back-calculation model

This function runs the CD4 back-calculation model using Stan and runs
the post-processing.

## Usage

``` r
run_backcalc(
  stan_data,
  migration = NULL,
  rita = NULL,
  age = NULL,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  iter = 2000,
  chains = 4,
  refresh = 100,
  messages = TRUE,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE,
  return_fit = TRUE
)
```

## Arguments

- stan_data:

  A list of data formatted for Stan.

- migration:

  Logical indicating whether to include migration in the model.

- rita:

  Logical indicating whether to include RITA data in the model.

- age:

  Logical indicating whether to include age-specific data in the model.

- inf_model:

  Infection model type (1 = spline, 2 = random walk, 3 = GP).

- diag_model:

  Diagnosis model type (1 = spline, 2 = random walk).

- model_seed:

  Seed for model fitting (will be set randomly if NULL).

- iter:

  Number of iterations for Stan sampling.

- chains:

  Number of chains for Stan sampling.

- refresh:

  Number of iterations between progress messages.

- messages:

  Logical indicating whether to show messages during sampling.

- start_yr:

  Starting year for posterior predictive quantities.

- nworkers:

  Number of parallel workers to use.

- chain_results:

  Logical indicating whether to return results by chain.

- return_fit:

  Logical indicating whether to return the fit object

## Value

List of posterior predictive quantities and (optionlly) the stanfit
object.

## Examples

``` r
sim_diags <- simulate_diagnoses(sim_seed = 123)
hiv_list <- run_backcalc(sim_diags, iter = 5, chains = 1, messages = FALSE, nworkers = 1)
#> 
#> SAMPLING FOR MODEL 'AI_spl' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000573 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 5.73 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration: 1 / 5 [ 20%]  (Warmup)
#> Chain 1: Iteration: 3 / 5 [ 60%]  (Sampling)
#> Chain 1: Iteration: 5 / 5 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.008 seconds (Total)
#> Chain 1: 
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
```
