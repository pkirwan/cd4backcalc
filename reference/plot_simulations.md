# Plot simulated or observed quantities

Generate a ggplot2 object plotting simulated quantities over time. When
`sim_data` is absent, quantities available directly from real-data Stan
inputs are plotted instead.

## Usage

``` r
plot_simulations(sim_diags, quantity)
```

## Arguments

- sim_diags:

  A list containing simulated diagnosis data objects `sim_data` and
  `stan_data`, or a real-data result containing `stan_data` and no
  `sim_data`.

- quantity:

  A string specifying the quantity to plot. Options include "incidence",
  "incidence_cob", "undiag_prev", "undiag_migration", "all_migration",
  "diag_prob", "diagnoses", "migration_prob", "transfers",
  "transfers_detected", "ratio_abroad_uk", "detect_prob",
  "migrant_hiv_prevalence", "aids_diagnoses", "aids_diagnoses_age",
  "cd4_diagnoses", "cd4_diagnoses_age", "transfers_detected_age",
  "gamma_obs", "gamma_obs_group", and "recent_arrival_cd4".

## Value

A ggplot2 object.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
