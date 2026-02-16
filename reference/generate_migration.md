# Generate migration combining omega and mu

Generate migration combining omega and mu

## Usage

``` r
generate_migration(nquar, nlatst, omega, mu, dist = TRUE)
```

## Arguments

- nquar:

  Number of time quarters.

- nlatst:

  Number of latent states.

- omega:

  Vector of omega values over time.

- mu:

  Matrix of mu values over time.

- dist:

  Logical indicating whether to use multinomial distribution.

## Value

Matrix of migration counts over time.
