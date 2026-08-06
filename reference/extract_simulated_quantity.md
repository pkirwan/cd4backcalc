# Extract simulated quantities based on the provided configuration.

Extract simulated quantities based on the provided configuration.

## Usage

``` r
extract_simulated_quantity(sim_diags, quantity, deterministic = FALSE)
```

## Arguments

- sim_diags:

  The simulated diagnosis data containing stan_data and sim_data.

- quantity:

  The quantity to extract.

- deterministic:

  Logical indicating whether to re-run the simulation deterministically
  (dist = FALSE) to obtain the expected values (lambda) rather than the
  stochastic realisation. Useful for bias/MSE calculations across
  multiple seeds (default is FALSE).

## Value

A list of model estimates including undiagnosed prevalence, incidence,
migration, diagnosis probabilities, diagnoses, migration probabilities,
the born-abroad:UK-born ratio of UK-acquired incidence, transfers of
care, and CD4 strata.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
