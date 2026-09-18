# Generate recent-arrival CD4 counts

Generates annual age-group by CD4 counts among diagnoses of infections
present on arrival, restricted to arrival in the diagnosis year or the
preceding calendar year.

## Usage

``` r
generate_recent_arrival_cd4(m_age, d_probs_mig, q_age, completeness, dist)
```

## Arguments

- m_age:

  Undiagnosed migrant arrivals by age, quarter, and the four non-recent
  CD4 states.

- d_probs_mig:

  Migrant diagnosis probabilities by quarter, age, and the four non-AIDS
  CD4 states.

- q_age:

  Progression probabilities by age and state.

- completeness:

  Proportion of eligible diagnoses observed.

- dist:

  Logical indicating whether to add stochastic variation.

## Value

A list containing observed and expected annual count arrays.
