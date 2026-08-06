# Create initial values for Stan from a fitted model

This function extracts the last iteration's parameter values from a
fitted Stan model. It returns a list of initial values for each chain,
which can be used to initialize a new Stan model.

## Usage

``` r
create_stan_inits(fit)
```

## Arguments

- fit:

  A fitted Stan model object (e.g., from rstan or cmdstanr).

## Value

A list of initial values for each chain.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
