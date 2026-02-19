# cd4backcalc

This vignette demonstrates how to use the `cd4backcalc` package to
simulate HIV diagnosis data and fit a CD4-staged back-calculation model.

``` r
library(cd4backcalc)
```

## Simulate data

To demonstrate the process of fitting a back-calculation model we need
to first simulate some diagnosis data. Data are simulated by propagating
a number of HIV infections through the back-calculation model to
generate diagnoses over time.

We can define a simulation configuration using the `define_sim_config`
function. This function allows us to specify whether to include
migration status, RITA status, and/or age groups in the simulation, as
well as varying the patterns for each of the model parameters,
e.g. `h_pattern` which controls the pattern of HIV incidence over time
can be set to “constant”, “linear”, or “quadratic”.

In case real data are being used, then the simulation step can be
skipped, and the migration, rita, and age parameters can be passed
directly to the `run_backcalc` function with the formatted `stan_data`.

The `simulate_diagnoses` function then uses this configuration to
simulate HIV diagnoses over time.

``` r
## simulate data
stan_data <- simulate_diagnoses(migration = FALSE, rita = FALSE, age = FALSE)
```

## Fit model

Next, the `run_backcalc` function is used to fit the back-calculation
model to the simulated diagnosis data using the specified model
configuration. The model fitting uses the cmdstanr package and
additional options can be passed to control the MCMC sampling process
(e.g. number of iterations, number of chains, etc).

``` r
## fit model
hiv_list <- run_backcalc(stan_data)
```

## Assessing convergence

Convergence of the MCMC chains can be assessed using the `cmdstanr`
package functions such as `traceplot` and `stan_rhat`. Additionally, the
`cd4backcalc` package provides a convenience function `plot_diagnostics`
to plot key diagnostics for the fitted model.

``` r
## plot diagnostics
plot_diagnostics(hiv_list)
```

## Assessing goodness of fit

The `hiv_list` object returned from the post-processing can now be
compared to the simulated data using goodness of fit tests such as bias.
When multiple models have been fitted a complementary mean squared error
(MSE) plot can be used. This step can be omitted if models were not
fitted to simulated data.

``` r
## generate bias and MSE plots
bias_plot(hiv_list, quantity = "incidence")
```

## Model estimates

A selection of pre-written plotting functions are available to plot the
model estimates. These functions use `ggplot2` so the plots may be
modified by adding additional options.

``` r
## plot incidence estimates compared to simulated values
plot_estimates(hiv_list, quantity = "incidence")

## plot undiagnosed prevalence estimates compared to simulated values
plot_estimates(hiv_list, quantity = "undiag_prev")
```

## References

- Birrell PJ, Chadborn TR, Gill ON, Delpech VC, et al. (2012).
  Estimating trends in incidence, time-to-diagnosis and undiagnosed
  prevalence using a CD4-based Bayesian back-calculation. *Stat. Commun.
  Infect. Dis.* 4(1). doi:
  [10.1515/1948-4690.1055](https://doi.org/10.1515/1948-4690.1055).
- Brizzi F, Birrell PJ, Plummer MT, Kirwan PD, et al. (2019). Extending
  Bayesian back-calculation to estimate age and time specific HIV
  incidence. *Lifetime Data Anal.* 25(4), pp.757-780. doi:
  [10.1007/s10985-019-09465-1](https://doi.org/10.1007/s10985-019-09465-1).
