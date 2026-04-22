# Generate diagnosis probability patterns

Generate diagnosis probability patterns

## Usage

``` r
generate_diagnosis_probs(
  nquar,
  nage = 52,
  pattern,
  dist = TRUE,
  rita = FALSE,
  age_diag_effect = 0
)
```

## Arguments

- nquar:

  Number of time quarters.

- nage:

  Number of age groups.

- pattern:

  Pattern type for diagnosis probabilities ("constant", "varying").

- dist:

  Logical indicating whether to add normal variation.

- rita:

  Logical indicating whether to generate an additional diagnosis
  probability for recently acquired HIV.

- age_diag_effect:

  Numeric scalar controlling the strength of age variation in diagnosis
  probabilities. A value of 0 means no age effect; positive values
  increase diagnosis probability with age. Applied as a multiplicative
  factor on the logit scale.

## Value

A 3D array of diagnosis probabilities with dimensions
`[nquar, nage, nstates]`.
