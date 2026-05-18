# Generate infections according to specified pattern

Generate infections according to specified pattern

## Usage

``` r
generate_infections(nquar, nlatst, pattern, target_inf = 90000, dist = TRUE)
```

## Arguments

- nquar:

  Number of time quarters.

- nlatst:

  Number of latent states.

- pattern:

  Pattern type for infections ("varying", "varying_2", "constant",
  "increasing").

- target_inf:

  Target total number of infections to simulate (may be higher/lower
  with Poisson variation).

- dist:

  Logical indicating whether to add Poisson variation.

## Value

Matrix of new infections over time.
