# Changelog

## Changes in v0.2.4

#### Improvements

- Optimised Stan models by moving migration probability calculations
  outside of loop.
- Matched Stan models by removing AR(1) process from simulated migration
  probabilities.
- `cmdstanr` is now a suggested rather than required package, as it is
  not on CRAN.
- Package size reduced to under 5 MB and vignette figures regenerated.

## Changes in v0.2.3

#### Improvements

- New \`quantity=incidence_cob’ option to plot estimated HIV incidence
  by country of birth.
- Migration-model post-processing now returns `alpha_diag_mig`, the
  posterior summary of the shared born-abroad diagnosis offset (logit
  scale and odds ratio), which can be plotted with
  `plot_estimates(quantity = "alpha_diag_mig")`.

#### Bug fixes

- Plots will now use the stored `start_year` value, as defined when
  running the model fit, if available.

## Changes in v0.2.2

#### Improvements

- Emit undiag_prev_matrix from
  [`postproc_ad()`](https://pkirwan.github.io/cd4backcalc/reference/postproc_ad.md)
  to simplify defining initial prevalence from existing runs.
- Fix a failing text on macOS.

## Changes in v0.2.1

#### Improvements

- Plotting and labelling improvements for
  [`plot_simulations()`](https://pkirwan.github.io/cd4backcalc/reference/plot_simulations.md)
  and
  [`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md).
- Updated target_inf and target_mig in \`simulate_diagnoses()’ to be
  consistent.
- Updated vignettes to show latest model DAGs.

## Changes in v0.2.0

#### Improvements

- The simplified age-dependent migration model is now the default, with
  a re-parameterisation of transfers of care to align to
  epidemiologically relevant quantities.
- Migration probability simulation patterns now mirror the Stan model
  and the redundant `m_pattern="reduced"` option has been removed.
- Ported the optional `gamma_obs` into the age-independent migration
  models.
- Incorporated efficiency improvements from age-dependent migration
  models into other age-dependent Stan models.
- New `gamma_obs_drift_sd` option in
  [`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md)
  to control the log-scale random-walk variation on
  `migrant_hiv_prevalence`.
- [`plot_simulations()`](https://pkirwan.github.io/cd4backcalc/reference/plot_simulations.md)
  now also plots observed quantities from real-data results without
  `sim_config`, e.g. `aids_diagnoses_age`, `gamma_obs_group`, and
  `cd4_diagnoses_age`.
- [`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md)
  now accepts a `colours` argument with the model or chain colours.
- Re-ran R vignettes with improved plotting and model logic.

#### Bug fixes

- Fixed the simulation of old-age AIDS cases using a multinomial draw
  rather than simple rounding logic.
- Fixed the reconstruction of age-independent initial prevalence in C++
  code.

## Changes in v0.1.5

#### Improvements

- This version of the package implements a simplified age-dependent
  migration Stan model, to overcome identifiability issues with the more
  complex models.
- New simulation option `recent_arrival_completeness` controls the
  proportion of recently arrived migrants with non-recently-acquired HIV
  for fitting the reduced model.
- New simulation option `m_pattern="reduced"` chooses the simplified
  migration pattern, which may be adopted for future releases.

#### Bug fixes

- Addressed a rounding issue in the construction of ABROAD_age_mat in
  simulated data.

## Changes in v0.1.4

#### Improvements

- Updates to migration models to improve fit against real-world data.
- Additional parameters for
  [`define_sim_config()`](https://pkirwan.github.io/cd4backcalc/reference/define_sim_config.md)
  to control:
  - change in born-abroad arrival age composition over time
    (`abroad_age_trend`);
  - ascertainment step change in `detect_prob` (`detect_step_quarter`
    and `detect_step_size`) and;
  - drift in the HIV-positive share of the external migration
    (`gamma_obs_scale_trend`).
- `years` option added to
  [`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
  and
  [`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
  functions to restrict to a subset of the data.
- Added News section to the pkgdown website and various improvements to
  comments and documentation.

#### Bug fixes

- Correct `init_prev_calibration()` to prevent zeros in the calibration
  of initial prevalence.

## Changes in v0.1.3

#### Improvements

- An optional rita_completeness_by_state vector can be passed to
  simulate_diagnoses() to allow modelling of informative missingness in
  the rita compartment.
- The age profile of migrants was previously modelled deterministically
  from the data, now age_group_prob will be modelled from the data.
- Update gamma_obs_group to feed in age-group-specific migration to the
  Stan model, where available, or gamma_obs to provide
  non-age-group-specific data. The logical scalar/vector
  gamma_obs_grouped is provided to simulate_diagnoses() to determine
  when to use grouped vs. ungrouped data.

#### Bug fixes

- Fix legacy scaling factor in derive_init_props which was applying an
  arbitrary multiplier to initial probabilities.

## Changes in v0.1.2

#### Improvements

- Derive initial state conditions analytically to aid fit during initial
  time points.
- Fit migration parameters using lower-dimensional splines/curves to
  improve identifiability.
- Allow an external overall migration time-series with rolling-annual
  counts to be included to anchor the migration data (via gamma_obs_sd
  and gamma_obs_scale arguments).
- Print the name of the checkpoint directory when checkpointing is
  enabled.

## Changes in v0.1.1

#### Improvements

- Added linetype to plots to distinguish between models.
- Added transfers_detected to the age-independent migration-specific
  output (see v0.1.0 for details on usage).

#### Bug fixes

- Fixed the calculation of AIDS counts in simulated age-specific data
  with migration.

## Changes in v0.1.0

#### Improvements

- New calibrate_init_prev() function which uses the fitted first quarter
  diagnosis probabilities and observed first quarter diagnoses to
  calculate improved initial prevalence values in case of poor fit to
  the data in the first few quarters.
- Added quantity = “transfers_detected” option for plot_estimates() to
  allow comparison against detected transfers of care, rather than all
  transfers of care. When comparing to observed data, the posterior
  predictive distribution is shown (as for HIV/AIDS diagnoses).
- Allow multiple fits to be compared in bias and MSE plots, with
  age-specific calculations supported using quantity = “incidence_age”.

## Changes in v0.0.8

#### Improvements

- Added detect_obs_sd, detect_obs_times, and eta_sd as parameters for
  the simulate_diagnoses() function, allowing control of the variance
  and sparsity in the detection probability and (native born:born
  abroad) ratio which is passed to the Stan migration models.
- Updated prior specification in the age-dependent Stan models.
- Revised incidence and migration patterns for combo_4 to provide a
  qualitatively distinct pattern to combo_3.
- Updated figures in the README and added figures to the vignettes.

#### Bug fixes

- Fixed a trace plot display issue in plot_diagnostics()
- Fixed the calculation of age-specific incidence in post-processing
  which was previously excluding born-abroad, UK-acquired.

## Changes in v0.0.7

#### Improvements

- New init option for run_backcalc() allows parameter values from a
  previous run to be passed to a new run for fine-tuning. Use inits \<-
  create_stan_inits(fit) to extract the parameter values, see
  vignette(“checkpointing”) for further details.
- Updated diagnosis probability age-effects in age-dependent Stan models
  to use zero-sum constraints and shared shrinkage to aid
  identifiability.
- Vignettes updated to provide more comprehensive details about model
  fitting and the quantities available.

## Changes in v0.0.6

#### Improvements

- Switched to an auto-regressive, AR(1), process for the diagnosis
  probabilities across all Stan models to mitigate identifiability
  issues in the most recent time-points.
- Simulation of diagnosis probabilities is now also specified according
  to a comparable AR(1) process.
- Updated the diagnosis probability and migration probability priors to
  use a non-centered paramaterisation to improve model fit.
- Revised several other priors to improve consistency and
  identifiability across models and to more closely align with
  real-world values.

#### Bug fixes

- Fixed the implementation of the migration-adjusted Stan models so that
  misclassified migration is correctly accounted for within the
  multinomial likelihood.
- Fixed post-processing to use migration-specific diagnosis
  probabilities where appropriate.
- Added posterior quantiles for undiagnosed migration by age band to the
  post-processing output.
- For the migration models, age-independent post-processing now uses the
  latent prevalence directly from the model (as in the age-dependent
  model). This resolved an issue with inflated undiagnosed prevalence
  due to misclassified migration arrivals.

## Changes in v0.0.5

#### Improvements

- Stan model updates:
- Replaced the random walk with a spline in age-dependent Stan models
  for HIV diagnosis probabilities for consistency with age-independent
  models.
- Added false recency rate to RITA-specific Stan models with an
  informative prior.
- Added a diagnosis state-specific intercept to migration-adjusted Stan
  models to account for differing probability of diagnosis in the
  born-abroad group.
- Updated priors across all Stan models for consistency.
- Changed expected-value simulation from piecewise linear to a smooth
  curve.
- Added realistic noise to all simulated values.
- Added new age_diag_effect, frr, and mig_diag_offset arguments to
  simulate_diagnoses:
- age_diag_effect applies a logit-additive offset to control the
  strength of age variation in diagnosis probabilities. The default
  (age_diag_effect = 0) applies no age variation.
- frr applies a false recency rate for RITA models, with a default of
  1.9%.
- mig_diag_offset applies a state-specific logit-scale offset to the
  diagnosis probabilities for the born-abroad group. The default
  (mig_diag_offset = NULL) applies no offset, offsets should be provided
  as a vector of values, e.g. mig_diag_offset = c(-0.4, -0.3, -0.2,
  -0.1).

#### Bug fixes

- Improved the simulation of misclassified transfers from abroad in the
  simulate_diagnoses() function to better align with the Stan model
  assumptions.
- Reduced simulated quarterly diagnosis probabilities to more realistic
  values which are in line with the model priors.

## Changes in v0.0.4

#### Improvements

- Parameterisation of migration probabilities in migration-adjusted Stan
  models now uses a softmax transformation instead of a sequential logit
  (stick-breaking) approach.
- Numerical safeguard added to age-dependent Stan models to reduce the
  chance of NaN being passed to the multinomial likelihood.
- Improved plotting defaults (colours, axis labels, facet layout).
- Expanded vignettes and added roxygen documentation for key exported
  functions on the pkgdown site.

#### Bug fixes

- Checkpoint functions now handle chain dropout that may occur during
  long MCMC runs.
- Plotting functions no longer error when list elements are NULL.
- Corrected age-group redistribution logic in simulate_diagnoses()
  function.

## Changes in v0.0.3

#### Major changes

- Inspired by the (archived) <https://venpopov.com/chkptstanr/>, this
  release provides checkpointing functionality to allow MCMC fits to be
  resumed when working on a cluster which imposes time limits on single
  jobs. To enable checkpointing the following options should be provided
  to run_backcalc:
- checkpoint_dir to define the location to save intermediate model
  states, and;
- iter_per_chunk to define the number of sampling iterations to perform
  for each chunk (applies to warmup and sampling phases.)

#### Bug fixes

- Fixed a crash when specifying nage smaller than 52.
- Stan should issue fewer complaints about non-symmetric precision
  matrices.

#### Improvements

- New combo_4 simulation pattern with qualitatively different incidence
  and migration patterns.
- Plotting functions now accept a second model fit (hiv_list_2) and
  model labels (hiv_name_1 and hiv_list_2) for comparison.
- A writable location will be chosen for compiled CmdStan executables.
  By default a job-local temporary directory will be set to avoid
  failures from shared home-directory caches on batch systems, but this
  can be changed via the CD4BACKCALC_CMDSTAN_DIR environment variable.

#### Dependencies

- Removed reliance on specific package versions where possible to
  improve backwards compatibility to R 4.1.0.
- Install instructions now recommend pak::pak() following the
  deprication of devtools::install_github().

## Changes in v0.0.2

#### Breaking changes

- run_backcalc() and fit_stan_model() now use iter_warmup and
  iter_sampling arguments (replacing the single iter argument), in line
  with the cmdstanr interface. Both default to 1000.

#### Major changes

- Stan backend migrated from rstan to cmdstanr. This brings compiler
  improvements and better performance.
- Stan models rewritten to reduce computational overhead and improve
  sampling speed across all model variants.
- Models are now compiled at first use rather than at package
  installation. This substantially shortens install time and sidesteps
  an execute permissions issue on some systems. The compiled binaries
  are saved to cmdstanr::cmdstan_default_install_path() for faster
  subsequent loads.

#### Dependencies

- Added: cmdstanr (\>= 0.9.0), sourced from the Stan R packages
  repository (<https://stan-dev.r-universe.dev>).
- Removed: rstan and rstantools.

## Changes in v0.0.1

#### Initial release

- Initial release of cd4backcalc including age-independent and
  age-specific models for CD4-only, CD4+RITA, and migration model
  variants.
- Expected workflow:
- Simulate data using simulate_diagnoses() (or format real data for
  model fitting);
- Fit appropriate model using run_backcalc() (if running multiple fits
  store these in a list);
- Examine traceplots and model fit using plot_diagnostics() and
  bias_plot();
- Plot estimates for specific quantities using plot_estimates().

#### Dependencies

- rstantools and Rcpp provide the skeleton for including Stan models and
  C++ code within an R package.
