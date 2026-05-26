# Checkpointing and long runs

Age-dependent back-calculation models can take a long time to fit. On
resource-constrained systems (e.g. SLURM clusters with time limits),
checkpointing allows you to fit the model in stages, saving intermediate
results that can be resumed if the job is interrupted.

``` r

library(cd4backcalc)
```

## Basic checkpointing

Pass `checkpoint_dir` and `iter_per_chunk` to
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
to enable checkpointing. The model will fit in chunks of
`iter_per_chunk` sampling iterations, saving a checkpoint after each
chunk.

``` r

sim_diags <- simulate_diagnoses(age = TRUE, h_pattern = "increasing")

hiv_list <- run_backcalc(
  sim_diags,
  iter_warmup = 1000,
  iter_sampling = 2000,
  checkpoint_dir = "checkpoints",
  iter_per_chunk = 500
)
```

This will run:

1.  2 chunks of 500 warmup iterations each, with a checkpoint saved
    after each chunk
2.  4 chunks of 500 sampling iterations each, with a checkpoint saved
    after each chunk

## Resuming interrupted fits

If a job is interrupted, simply re-run the same code.
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
will detect existing checkpoints and resume from the last completed
chunk:

``` r

# re-run the same code — automatically resumes from checkpoint
hiv_list <- run_backcalc(
  sim_diags,
  iter_warmup = 1000,
  iter_sampling = 2000,
  checkpoint_dir = "checkpoints",
  iter_per_chunk = 500
)
```

## Creating initial values from a previous fit

For fine-tuning or extending a previous fit, use
[`create_stan_inits()`](https://pkirwan.github.io/cd4backcalc/reference/create_stan_inits.md)
to extract the last iteration’s parameter values and pass them to a new
run:

``` r

# fit a short initial run
short_fit <- run_backcalc(
  sim_diags,
  iter_warmup = 500,
  iter_sampling = 500
)

# extract initial values
inits <- create_stan_inits(short_fit$fit)

# reuse as starting point for a longer run
long_fit <- run_backcalc(
  sim_diags,
  iter_warmup = 500,
  iter_sampling = 1000,
  init = inits
)
```

The same `inits` object can also be passed directly to a custom
`cmdstanr` workflow. For interrupted chunked runs within `cd4backcalc`,
restarting via `checkpoint_dir` is usually simpler because it also
preserves the checkpointed adaptation state.

## When to use checkpointing

Checkpointing is most useful for:

- **Age-dependent models** which can take hours or days to fit;
- **SLURM or HPC environments** with wall-time limits;
- **Iterative model development** where you want to extend a run without
  starting from scratch.

For age-independent models, fitting is typically fast enough that
checkpointing is unnecessary.

## Single-chunk behaviour

If `iter_per_chunk` is greater than or equal to `iter_sampling`, no
chunking occurs and the model fits in a single run.
