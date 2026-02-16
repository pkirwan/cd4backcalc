# Function to get posterior functions of incidence

This function generates incidence by age group from log-infections (per
each posterior iteration).

## Usage

``` r
infs.fct(infs, age.st, age.end, nt, type = "time", nage = 52, start_yr = 0)
```

## Arguments

- infs:

  A matrix with quarterly log-infections for each posterior sample in
  each row.

- age.st:

  Starting age for incidence calculation.

- age.end:

  Ending age for incidence calculation.

- nt:

  Total number of time points (quarters).

- type:

  Type of incidence calculation, either "time" or "age".

- nage:

  Number of age groups.

- start_yr:

  Starting year for the time series.

## Value

A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time
point or age group.
