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
  nquar = 120,
  nage = 52,
  q = c(0.0976, 0.115, 0.116, 0.148),
  target_inf = 90000,
  target_mig = 12000,
  init_props = c(7, 3.5, 1.75, 0.5),
  init_props_abroad = c(1.75, 2.25, 2.25, 1.25),
  abroad_prop = 0.2,
  cd4_completeness = 1,
  rita_completeness = 1,
  age_diag_effect = 0,
  mig_diag_offset = NULL,
  frr = 0.019,
  detect_obs_times = NULL
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

  Pattern type for proportion UK ("constant", "increasing").

- m_pattern:

  Pattern type for mu ("constant", "differing", "varying").

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

- init_props:

  Initial state proportions.

- init_props_abroad:

  Initial state proportions for infections acquired abroad.

- abroad_prop:

  Proportion of infections acquired abroad at the outset.

- cd4_completeness:

  CD4 completeness.

- rita_completeness:

  RITA completeness.

- age_diag_effect:

  Numeric scalar controlling the strength of age variation in diagnosis
  probabilities. Applied as a logit-additive offset scaled linearly
  across ages.

- mig_diag_offset:

  Numeric vector of per-state logit-scale offsets applied to diagnosis
  probabilities for the migrant group (y=2). Length must equal the
  number of diagnosis states (4 for non-RITA, 5 for RITA). NULL
  (default) applies no offset.

- frr:

  False recency rate for RITA misclassification.

- detect_obs_times:

  Integer vector of quarter indices (1-based) at which external
  detection-probability observations are generated.

## Value

Configuration list.
