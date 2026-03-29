# Run CD4 back-calculation model

This function runs the CD4 back-calculation model using Stan. When
`checkpoint_dir` is provided, fitting is split into a warmup phase
followed by one or more sampling chunks of size `iter_per_chunk`. After
each chunk a checkpoint is written so that a timed-out SLURM job can be
resubmitted and will resume from the last completed chunk.

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
  messages = TRUE,
  checkpoint_dir = NULL,
  iter_per_chunk = NULL
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

- checkpoint_dir:

  Directory for saving/loading checkpoints. If `NULL` (default), no
  checkpointing is performed.

- iter_per_chunk:

  Number of sampling iterations per chunk when checkpointing. If `NULL`,
  all sampling iterations are run in a single chunk.

## Value

A CmdStanFit object (without checkpointing) or a
`cd4backcalc_combined_fit` object (with checkpointing).

## Examples

``` r
sim_diags <- simulate_diagnoses(sim_seed = 123)
stan_data <- model_type_prep(sim_diags)
fit <- fit_stan_model(
 stan_data, iter_warmup = 10, iter_sampling = 10, chains = 1, refresh = 0, messages = FALSE
)
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
