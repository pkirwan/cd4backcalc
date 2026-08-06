# Get bias between HIV model estimates and simulation results.

Get bias between HIV model estimates and simulation results.

## Usage

``` r
get_bias(hiv_list, type, deterministic = FALSE, years = NULL)
```

## Arguments

- hiv_list:

  A list of HIV model results.

- type:

  The variable type to calculate bias for.

- deterministic:

  Logical indicating whether to compare against deterministic generative
  values.

- years:

  Optional numeric vector specifying which years to include. If `NULL`,
  all available years after year 0 are used.

## Value

A data frame of bias values.
