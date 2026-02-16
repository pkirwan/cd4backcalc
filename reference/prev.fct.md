# Function to get posterior functions of prevalence

This function generates prevalence by age group from prevalence arrays
(per each posterior iteration).

## Usage

``` r
prev.fct(prev.arr, state, age.st, age.end, rita.flag, nt, start_yr)
```

## Arguments

- prev.arr:

  A 4D array with prevalence for each posterior sample.

- state:

  State for prevalence calculation, either "all" or 1 to 7.

- age.st:

  Starting age for prevalence calculation.

- age.end:

  Ending age for prevalence calculation.

- rita.flag:

  Logical flag indicating if RITA diagnoses are considered.

- nt:

  Total number of time points.

- start_yr:

  Starting year for the time series.

## Value

A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time
point.
