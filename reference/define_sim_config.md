# Define simulation configuration

Defines the patterns for various parameters in the epidemic simulation.
The number of latent and diagnosis states will be inferred from the
simulation specification.

## Usage

``` r
define_sim_config(
  h_pattern = "constant",
  o_pattern = "constant",
  d_pattern = "constant",
  p_pattern = "constant",
  m_pattern = "constant",
  r_pattern = "constant",
  a_pattern = "constant",
  sim_type = NULL,
  sim_seed = NULL,
  dist = TRUE,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  nquar = 40,
  nage = 52,
  q = c(0.0976, 0.115, 0.116, 0.148),
  target_inf = 30000,
  target_mig = 6000,
  abroad_prop = 0.2,
  cd4_completeness = 1,
  rita_completeness = 1,
  rita_completeness_by_state = NULL,
  age_diag_effect = 0,
  mig_diag_offset = NULL,
  frr = 0.019,
  detect_obs_times = NULL,
  gamma_obs_grouped = FALSE,
  gamma_obs_sd = 0,
  gamma_obs_scale = 500,
  abroad_age_trend = 0,
  detect_step_quarter = NULL,
  detect_step_size = 0,
  gamma_obs_scale_trend = 0,
  gamma_obs_drift_sd = 0.02,
  recent_arrival_completeness = 0
)
```

## Arguments

- h_pattern:

  Pattern type for infections ("constant", "increasing", "varying",
  "varying_2").

- o_pattern:

  Pattern type for omega ("constant", "increasing", "step_change_1",
  "step_change_2", "none").

- d_pattern:

  Pattern type for diagnosis probabilities ("constant", "varying").

- p_pattern:

  Pattern for the born-abroad:UK-born ratio of UK-acquired incidence
  ("constant", "increasing").

- m_pattern:

  Pattern type for mu ("constant", "differing", or "varying").

- r_pattern:

  Pattern type for RITA data ("constant", "increasing", "none").

- a_pattern:

  Pattern type for age-based data ("constant", "wedge",
  "varying_wedge").

- sim_type:

  If provided, a preset type of simulation to define patterns for.

- sim_seed:

  Random seed for simulation. If NULL, a random seed is generated.

- dist:

  Logical indicating whether to add stochastic variation.

- migration:

  Logical indicating whether to include migration effects.

- rita:

  Logical indicating whether to simulate RITA (Recent Infection Testing
  Algorithm) data.

- age:

  Logical indicating whether to simulate age-stratified data.

- nquar:

  Number of time quarters.

- nage:

  Number of age groups.

- q:

  Vector of progression probabilities.

- target_inf:

  Target total number of HIV acquisitions.

- target_mig:

  Target total number of migrations.

- abroad_prop:

  Proportion of infections acquired abroad at the outset.

- cd4_completeness:

  CD4 completeness.

- rita_completeness:

  RITA completeness.

- rita_completeness_by_state:

  Optional numeric vector of length 5 giving state-specific completeness
  for RITA simulations. Elements correspond to RITA positive, CD4 500 or
  above, CD4 350–499, CD4 200–349, and CD4 below 200. Use lower values
  in lower CD4 states to simulate informative missingness.

- age_diag_effect:

  Numeric scalar controlling the strength of age variation in diagnosis
  probabilities. Applied as a logit-additive offset scaled linearly
  across ages.

- mig_diag_offset:

  Numeric vector of per-state logit-scale offsets applied to diagnosis
  probabilities for people born abroad (`y = 2`). Length must equal the
  number of diagnosis states (4 for non-RITA, 5 for RITA). `NULL`
  (default) applies no offset.

- frr:

  False recency rate for RITA misclassification.

- detect_obs_times:

  Integer vector of quarter indices (1-based) at which external
  detection-probability observations are generated.

- gamma_obs_grouped:

  Logical scalar or logical vector of length `nquar` controlling which
  quarters receive age-stratified external migration observations on the
  same 4 age groups as `ABROAD_group`. The scalar `gamma_obs` series is
  always simulated; quarters where `gamma_obs_grouped` is `FALSE` have
  `-1` rows in `gamma_obs_group`.

- abroad_age_trend:

  Total logit-scale change (start to end) in the 15-24/25-34 age-group
  shares of born-abroad arrivals relative to the 45+ reference group; 0
  = time-constant age composition.

- detect_step_quarter:

  Quarter index (1-based) of an ascertainment step in detect_prob (e.g.
  a surveillance reporting change). `NULL` (default) means no step.

- detect_step_size:

  Logit-scale jump in detect_prob applied from `detect_step_quarter`
  onward.

- gamma_obs_scale_trend:

  Total log-scale change (start to end) in the HIV-positive share of the
  external migration index (gamma_obs_scale); 0 = constant share.

- gamma_obs_drift_sd:

  Quarterly random-walk innovation SD on the log ratio of the external
  migration index to HIV-positive migration. This controls stochastic
  variation in migrant HIV prevalence and is passed to the matching Stan
  model. Set to 0 to disable.

- recent_arrival_completeness:

  Proportion of eligible diagnoses with avidity testing done and a
  recorded year of arrival. Combined multiplicatively with
  `cd4_completeness` (a recorded CD4 value is also required) to give the
  overall completeness of the recent-arrival CD4 auxiliary data. Set to
  0 to disable.

## Value

Configuration list.
