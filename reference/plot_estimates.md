# Plot estimated model quantities

Generate a ggplot2 object plotting estimated quantities over time,
optionally compared to the simulated values.

## Usage

``` r
plot_estimates(hiv_list, quantity, start_year = 0, chain_results = FALSE)
```

## Arguments

- hiv_list:

  A list returned from the `postproc` function containing model
  estimates.

- quantity:

  A string specifying the quantity to plot. Options are "incidence",
  "undiag_prev", "undiag_migration", "all_migration", "diag_prob",
  "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and
  "detect_prob".

- start_year:

  An integer specifying the starting year for the x-axis (default is 0).

- chain_results:

  A logical indicating whether to plot chain-by-chain estimates (default
  is FALSE).

## Value

A ggplot2 object.
