
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cd4backcalc

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pkirwan/cd4backcalc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pkirwan/cd4backcalc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The cd4backcalc package provides a user-friendly way to perform
Bayesian, discrete-time, multi-state modelling to estimate HIV incidence
and undiagnosed HIV prevalence from diagnosis data.

## Description

CD4 back-calculation is a discrete-time, multi-state Markov model used
to represent HIV progression through different states until a diagnosis.
The undiagnosed (latent) states in the model represent decreasing CD4
strata, and may include recent HIV acquisition (defined using RITA
evidence), migration status, and age group. An illustration of the
latent and diagnosed states in a CD4-staged back-calculation model is
shown below.

Models are fitted within a Bayesian framework (using rstan), with
suitable prior distributions, to estimate HIV incidence and prevalence
of undiagnosed HIV over time, as well as the probabilities of diagnosis
in each CD4 state.

<div class="figure">

<img src="man/figures/README-backcalc-model-1.png" alt="CD4-staged back-calculation model. The model compartments represent undiagnosed HIV (round boxes) and diagnosed HIV or AIDS (square boxes). The $h_j$ arrow represents HIV acquired at each time point ($j$). HIV progression probabilities ($q_i$) vary by CD4 stratum. HIV diagnosis probabilities ($d_{i,j}$) vary by CD4 stratum and time." width="100%" />
<p class="caption">

CD4-staged back-calculation model. The model compartments represent
undiagnosed HIV (round boxes) and diagnosed HIV or AIDS (square boxes).
The $h_j$ arrow represents HIV acquired at each time point ($j$). HIV
progression probabilities ($q_i$) vary by CD4 stratum. HIV diagnosis
probabilities ($d_{i,j}$) vary by CD4 stratum and time.
</p>

</div>

## Installation

You can install the development version of cd4backcalc from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("pkirwan/cd4backcalc")
```

### System requirements for installation

- Windows users will require Rtools to compile the Stan models, see:
  <https://cran.r-project.org/bin/windows/Rtools/>.
- MacOS users will require the correct GNU Fortran compiler, see:
  <https://mac.r-project.org/tools/>.

## Supported models

The cd4backcalc package includes several pre-compiled models which can
be used to fit the back-calculation model to real or simulated diagnosis
data under different scenarios:

- CD4-staged back-calculation
- Age-dependent CD4-staged back-calculation
- Dual-biomarker back-calculation with CD4 and RITA data
- CD4-staged back-calculation with migration status

Most of the functions in the package are designed to work with all model
variants, with the appropriate model selected based on user-specified
options. Model variants are designed to be modular, e.g. an
age-dependent dual-biomarker model.

## Usage

A basic workflow to simulate diagnosis data and fit the CD4-staged
back-calculation model, check convergence and goodness of fit to the
simulated data, and plot model estimates is shown in the example below.

``` r
library(cd4backcalc)

# simulate data
sim_diags <- simulate_diagnoses(migration = FALSE, rita = FALSE, age = FALSE)

# run back-calculation model
hiv_list <- run_backcalc(stan_data)

# check convergence
plot_diagnostics(hiv_list)

# check goodness of fit
bias_plot(hiv_list, quantity = "incidence")

# plot estimates compared to simulated values
plot_estimates(hiv_list, quantity = "incidence")
plot_estimates(hiv_list, quantity = "undiag_prev")
```
