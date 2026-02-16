# Generate diagnosis probability patterns

Generate diagnosis probability patterns

## Usage

``` r
generate_diagnosis_probs(nquar, pattern, dist = TRUE, rita = FALSE)
```

## Arguments

- nquar:

  Number of time quarters.

- pattern:

  Pattern type for diagnosis probabilities ("constant", "varying").

- dist:

  Logical indicating whether to add normal variation.

- rita:

  Logical indicating whether to generate an additional diagnosis
  probability for recently acquired HIV.

## Value

Matrix of diagnosis probabilities over time.
