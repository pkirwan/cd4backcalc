# Define simulation type

This is a helper function to define various pre-specified simulation
types for epidemic simulations. Each configuration specifies patterns
for infections, migration, diagnosis probabilities, and other
parameters.

## Usage

``` r
define_simulation_type(config, sim_type = "constant")
```

## Arguments

- config:

  Configuration list to modify.

- sim_type:

  Type of simulation to define patterns for.

## Value

Configuration list for the specified simulation type.

## Details

List of available simulation types:

- "constant": All parameters constant over time.

- "h_increasing": Infections increase over time.

- "h_varying": Infections vary over time.

- "o_increasing": Migration (omega) increases over time.

- "o_step_change": Migration (omega) has step changes over time.

- "d_varying": Diagnosis probabilities vary over time.

- "p_increasing": Proportion UK increases over time.

- "m_differing": Migration probabilities (mu) differ by diagnosed state.

- "m_varying": Migration probabilities (mu) vary over time.

- "combo_1": Combination of varying infections, step change migration,
  and varying diagnosis probabilities.

- "combo_2": Combination of varying infections, step change migration,
  varying diagnosis probabilities, and increasing proportion UK.

- "combo_3": Combination of varying infections, step change migration,
  varying diagnosis probabilities, increasing proportion UK, and varying
  migration probabilities.

- "combo_3_no_migration": Similar to "combo_3" but without migration
  effects.

- "combo_4": As "combo_3" but with different incidence and migration
  patterns.
