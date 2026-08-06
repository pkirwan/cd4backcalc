# Make a list of mean squared error (MSE) values for HIV model estimates.

Make a list of mean squared error (MSE) values for HIV model estimates.

## Usage

``` r
make_mse_list(
  hiv_list,
  sim_comparison = c("sampled", "generative"),
  years = NULL
)
```

## Arguments

- hiv_list:

  A list of HIV model results.

- sim_comparison:

  A string specifying which simulated values to compare against. Options
  are `"sampled"` (the stochastic realisation for the particular seed,
  the default) or `"generative"` (deterministic expected values).

- years:

  Optional numeric/integer vector specifying which years to include in
  bias/MSE calculations. If `NULL` (default), all available years are
  used.

## Value

A list of bias and MSE values for incidence, migration, the
born-abroad:UK-born ratio of UK-acquired incidence, and migration
probabilities.
