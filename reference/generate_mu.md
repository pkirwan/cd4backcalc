# Generate migration probability patterns

Each pattern combines a time-varying transfer share with a time-constant
conditional CD4 composition among undiagnosed arrivals, matching the
Stan model parameterisation.

## Usage

``` r
generate_mu(nquar, nlatst, pattern)
```

## Arguments

- nquar:

  Number of time quarters.

- nlatst:

  Number of latent states.

- pattern:

  Pattern type for mu ("constant", "differing", or "varying").

## Value

Matrix of mu values over time.
