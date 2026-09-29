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
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 100,
  messages = TRUE,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE,
  return_fit = TRUE,
  checkpoint_dir = NULL,
  iter_per_chunk = NULL,
  init = NULL
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

  Diagnosis model type (1 = spline, 2 = AR(1)).

- model_seed:

  Seed for model fitting (will be set randomly if NULL).

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

- start_yr:

  Starting year for posterior predictive quantities.

- nworkers:

  Number of parallel workers to use.

- chain_results:

  Logical indicating whether to return results by chain.

- return_fit:

  Logical indicating whether to return the fit object

- checkpoint_dir:

  Directory for saving/loading checkpoints. If `NULL` (default), no
  checkpointing is performed.

- iter_per_chunk:

  Number of sampling iterations per chunk when checkpointing. If `NULL`,
  all sampling iterations are run in a single chunk.

- init:

  Initial values passed to Stan. When checkpointing, these are only used
  when starting a new run; resumed checkpointed fits use the saved
  checkpoint state.

## Value

List of posterior predictive quantities and (optionally) the stanfit
object.

## Examples

``` r
sim_diags <- simulate_diagnoses(sim_seed = 123)
hiv_list <- run_backcalc(
 sim_diags, iter_warmup = 10, iter_sampling = 10, chains = 1, messages = FALSE, nworkers = 1
)
#> CmdStan model cache dir: /tmp/RtmpgMTpkF/cd4backcalc-cmdstan
#> R tempdir(): /tmp/RtmpgMTpkF
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
