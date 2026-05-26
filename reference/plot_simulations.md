# Plot simulated quantities

Generate a ggplot2 object plotting simulated quantities over time.

## Usage

``` r
plot_simulations(sim_diags, quantity)
```

## Arguments

- sim_diags:

  A list containing the simulated diagnosis data objects sim_data and
  stan_data.

- quantity:

  A string specifying the quantity to plot. Options are "incidence",
  "undiag_prev", "undiag_migration", "all_migration", "diag_prob",
  "diagnoses", "migration_prob", "transfers", "transfers_detected",
  "ratio_abroad_uk", and "detect_prob".

## Value

A ggplot2 object.
