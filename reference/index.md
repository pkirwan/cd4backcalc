# Package index

## Simulation

Simulate HIV diagnosis data for model fitting

- [`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md)
  : Simulate HIV Epidemic Data
- [`sim_remodel()`](https://pkirwan.github.io/cd4backcalc/reference/sim_remodel.md)
  : Remodel the simulated diagnosis data for fitting a different model
  type.

## Model fitting

Fit and configure CD4-staged back-calculation models

- [`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
  : Run CD4 back-calculation model
- [`create_stan_inits()`](https://pkirwan.github.io/cd4backcalc/reference/create_stan_inits.md)
  : Create initial values for Stan from a fitted model
- [`calibrate_init_prev()`](https://pkirwan.github.io/cd4backcalc/reference/calibrate_init_prev.md)
  : Calibrate initial prevalence from fitted draws

## Visualisation

Plot model estimates, diagnostics, and goodness of fit

- [`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md)
  : Plot estimated model quantities
- [`plot_simulations()`](https://pkirwan.github.io/cd4backcalc/reference/plot_simulations.md)
  : Plot simulated quantities
- [`plot_diagnostics()`](https://pkirwan.github.io/cd4backcalc/reference/plot_diagnostics.md)
  : Plot diagnostics for fitted model
- [`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
  : Box plots of bias.
- [`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
  : Box plots of mean squared error (MSE).

## Data extraction

Extract model estimates and data for custom analysis

- [`get_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/get_estimates.md)
  : Get estimates from HIV model results for plotting.
- [`extract_simulated_quantity()`](https://pkirwan.github.io/cd4backcalc/reference/extract_simulated_quantity.md)
  : Extract simulated quantities based on the provided configuration.
- [`extract_real_data()`](https://pkirwan.github.io/cd4backcalc/reference/extract_real_data.md)
  : Extract real data for plotting
