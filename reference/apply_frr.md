# Apply false recency rate misclassification to a RITA CD4 array

Moves a fraction `frr` of long-standing (non-RITA-recent, states 2+)
diagnoses to the RITA-recent state (state 1), simulating false recency
classification errors in the RITA assay. Works on arrays of any rank
whose *last* dimension indexes the RITA/CD4 states.

## Usage

``` r
apply_frr(cd4_arr, frr, dist = TRUE)
```

## Arguments

- cd4_arr:

  Integer array whose last dimension indexes RITA/CD4 states

- frr:

  Scalar false recency rate in \\\[0, 1\]\\

- dist:

  Logical; if `TRUE` the number of misclassified cases per cell is drawn
  from a Binomial distribution, otherwise the deterministic expectation
  (rounded, preserving the total) is used.

## Value

An integer array with the same dimensions as `cd4_arr`.
