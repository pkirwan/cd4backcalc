# Generate migration probability patterns

Generate migration probability patterns

## Usage

``` r
generate_mu(nquar, nlatst, pattern, detect_prob)
```

## Arguments

- nquar:

  Number of time quarters.

- nlatst:

  Number of latent states.

- pattern:

  Pattern type for mu ("constant", "differing", "varying").

- detect_prob:

  Vector of detection probabilities for those migrating already
  diagnosed.

## Value

Matrix of mu values over time.
