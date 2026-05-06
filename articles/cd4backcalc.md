# Getting started

This vignette walks through a basic workflow using `cd4backcalc`:
simulate HIV diagnosis data, fit a CD4-staged back-calculation model,
assess convergence and goodness of fit, and plot the model estimates.

``` r

library(cd4backcalc)
```

## Simulate data

We first simulate HIV diagnosis data using
[`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md).
This function propagates a specified number of HIV infections through
the back-calculation model to generate diagnoses over time. Key
arguments control the incidence pattern (`h_pattern`), whether to
include migration, RITA, and/or age structure, and a seed for
reproducibility.

``` r

sim_diags <- simulate_diagnoses(h_pattern = "increasing", sim_seed = 123)
```

If you are working with real data rather than simulated data, you can
skip this step and pass a pre-formatted `stan_data` list directly to
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md).

We can visualise the simulated data using
[`plot_simulations()`](https://pkirwan.github.io/cd4backcalc/reference/plot_simulations.md):

``` r

plot_simulations(sim_diags, quantity = "incidence")
```

## Fit the model

The
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
function fits the CD4-staged back-calculation model using `cmdstanr`. By
default, each chain runs 1000 warmup and 1000 sampling iterations across
4 chains.

``` r

hiv_list <- run_backcalc(sim_diags)
```

For faster iteration during development, you can reduce the number of
iterations and chains:

``` r

hiv_list <- run_backcalc(
  sim_diags,
  iter_warmup = 500,
  iter_sampling = 500,
  chains = 2
)
```

## Assess convergence

Use
[`plot_diagnostics()`](https://pkirwan.github.io/cd4backcalc/reference/plot_diagnostics.md)
to inspect traceplots and posterior density plots for key parameters:

``` r

plot_diagnostics(hiv_list, ntraces = 4)
```

You can also access the underlying `cmdstanr` fit object for additional
diagnostics:

``` r

hiv_list$fit$summary()
hiv_list$fit$cmdstan_diagnose()
```

## Goodness of fit

When working with simulated data,
[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
compares the model estimates against the known simulated values:

``` r

bias_plot(hiv_list, quantity = "incidence")
```

## Plot estimates

[`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md)
produces plots of the posterior estimates with credible intervals. When
simulated data is available, the true values are overlaid for
comparison.

``` r

# HIV incidence over time
plot_estimates(hiv_list, quantity = "incidence")

# undiagnosed prevalence over time
plot_estimates(hiv_list, quantity = "undiag_prev")

# diagnosis probabilities
plot_estimates(hiv_list, quantity = "diag_prob")

# expected diagnoses vs observed
plot_estimates(hiv_list, quantity = "diagnoses")
```

All plotting functions return `ggplot2` objects, so they can be further
customised:

``` r

library(ggplot2)
plot_estimates(hiv_list, quantity = "incidence") +
  labs(title = "Estimated HIV incidence") +
  theme_minimal()
```

## Next steps

- See [Model
  types](https://pkirwan.github.io/cd4backcalc/articles/model-types.html)
  for fitting models with RITA evidence, migration, and age structure.
- See [Analysing
  results](https://pkirwan.github.io/cd4backcalc/articles/analysing-results.html)
  for extracting estimates, comparing models, and the full set of
  available quantities.
- See
  [Checkpointing](https://pkirwan.github.io/cd4backcalc/articles/checkpointing.html)
  for fitting large models in stages on resource-constrained systems.
