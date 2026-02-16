# Function to identify the peak in infection in the last x years

Function to identify the peak in infection in the last x years

## Usage

``` r
peak.fct(infs, age.st, age.end, nt, last.yr = 10, nage = 52)
```

## Arguments

- infs:

  A matrix with quarterly log-infections for each posterior sample in
  each row.

- age.st:

  Starting age for peak identification.

- age.end:

  Ending age for peak identification.

- nt:

  Total number of time points (quarters).

- last.yr:

  Number of last years to consider for peak identification.

- nage:

  Number of age groups.

## Value

A numeric vector indicating the peak position for each posterior sample.
