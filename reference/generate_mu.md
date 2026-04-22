# Generate migration probability patterns

Generate migration probability patterns

## Usage

``` r
generate_mu(nquar, nlatst, pattern, dist = FALSE)
```

## Arguments

- nquar:

  Number of time quarters.

- nlatst:

  Number of latent states.

- pattern:

  Pattern type for mu ("constant", "differing", "varying").

- dist:

  Logical indicating whether to add logit-scale noise.

## Value

Matrix of mu values over time.
