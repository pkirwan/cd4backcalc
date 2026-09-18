# Generate migration probability patterns

Each pattern combines a time-varying transfer share with a time-constant
conditional CD4 composition among undiagnosed arrivals, matching the
Stan model parameterisation. With `dist = TRUE`, smooth temporal
variation is added to the transfer share on the logit scale.

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

  Pattern type for mu ("constant", "differing", or "varying").

- dist:

  Logical indicating whether to add logit-scale temporal variation.

## Value

Matrix of mu values over time.
