# Round specified fields in an epidemic object

Round specified fields in an epidemic object

## Usage

``` r
round_fields(
  epi_obj,
  fields = c("lat", "diags", "lat_age", "diags_age"),
  stochastic_diags_age = FALSE
)
```

## Arguments

- epi_obj:

  Epidemic object containing fields to round.

- fields:

  Character vector of field names to round.

- stochastic_diags_age:

  Logical; if TRUE, allocate diagnosis totals across ages with
  multinomial draws.

## Value

Epidemic object with specified fields rounded.
