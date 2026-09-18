# Round age-stratified arrays to match pre-rounded quarter/state totals

Round age-stratified arrays to match pre-rounded quarter/state totals

## Usage

``` r
round_age_array_to_totals(arr, totals, stochastic = FALSE)
```

## Arguments

- arr:

  Numeric array with age in the first dimension.

- totals:

  Numeric array giving the target sums for each combination of the
  remaining dimensions.

- stochastic:

  Logical; if TRUE, allocate each total across ages with a multinomial
  draw instead of largest-remainder rounding.

## Value

Rounded array whose sums over age match totals exactly.
