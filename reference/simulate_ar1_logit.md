# Simulate AR(1) trajectory on the logit scale

Produces a time series with the same correlation structure as the Stan
AR(1) random-walk prior used for diagnosis probabilities. AR(1)
residuals are added on top of a deterministic trend and then
mean-centred so that the stochastic trajectory does not drift
systematically from the deterministic trend. This ensures that
`dist = TRUE` and `dist = FALSE` simulations produce comparable epidemic
trajectories (e.g. undiagnosed prevalence), while still exhibiting
realistic smooth quarter-to-quarter variation.

## Usage

``` r
simulate_ar1_logit(nquar, trend_logit, rho = 0.86, vardelta = exp(-4))
```

## Arguments

- nquar:

  Number of time quarters.

- trend_logit:

  Numeric vector of length `nquar` giving the deterministic mean
  trajectory on the logit scale.

- rho:

  AR(1) autocorrelation coefficient.

- vardelta:

  Innovation variance (step variance per quarter).

## Value

Numeric vector of length `nquar` on the logit scale.
