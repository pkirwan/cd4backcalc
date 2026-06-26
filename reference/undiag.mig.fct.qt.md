# Function to get posterior quantiles of undiagnosed migration by age band

Function to get posterior quantiles of undiagnosed migration by age band

## Usage

``` r
undiag.mig.fct.qt(draws, abroad_age_prop, age.st, age.end, nage, nt, start_yr)
```

## Arguments

- draws:

  A matrix of posterior draws for undiagnosed migration.

- abroad_age_prop:

  An array of posterior migrant age proportions indexed by draw,
  quarter, and age (dimensions: `[draw, quarter, age]`).

- age.st:

  Starting age for the band.

- age.end:

  Ending age for the band.

- nage:

  Number of single-year age groups.

- nt:

  Total number of time points (quarters).

- start_yr:

  Starting year for the time series.

## Value

A tibble with summary statistics (2.5%, 50%, mean, 97.5%) for each time
point.
