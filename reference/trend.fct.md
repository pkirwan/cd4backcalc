# Function to identify the trend in infection in the last x years

Function to identify the trend in infection in the last x years

## Usage

``` r
trend.fct(infs, age.st, age.end, nt, last.yr = 3, nage = 52)
```

## Arguments

- infs:

  A matrix with quarterly log-infections for each posterior sample in
  each row.

- age.st:

  Starting age for trend identification.

- age.end:

  Ending age for trend identification.

- nt:

  Total number of time points (quarters).

- last.yr:

  Number of last years to consider for trend identification.

- nage:

  Number of age groups.

## Value

A character vector indicating the trend for each posterior sample.
