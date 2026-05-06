# Make a list of mean squared error (MSE) values for HIV model estimates.

Make a list of mean squared error (MSE) values for HIV model estimates.

## Usage

``` r
make_mse_list(hiv_list, sim_comparison = c("generative", "sampled"))
```

## Arguments

- hiv_list:

  A list of HIV model results.

- sim_comparison:

  A string specifying which simulated values to compare against. Options
  are `"generative"` (deterministic expected values, the default) or
  `"sampled"` (the stochastic realisation for the particular seed).

## Value

A list of bias and MSE values for incidence, migration, proportion of
UK-acquired infections that are migrants, and migration probabilities.
