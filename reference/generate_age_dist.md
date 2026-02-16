# Generate age distributions

Generates age distributions for a given number of age categories over
the epidemic.

## Usage

``` r
generate_age_dist(nquar, nage = 52, pattern, peak_start = 10, peak_end = 30)
```

## Arguments

- nquar:

  Number of time quarters.

- nage:

  Number of age groups.

- pattern:

  Pattern type for age distribution ("constant", "wedge",
  "varying_wedge").

- peak_start:

  Peak age during first time quarter.

- peak_end:

  Peak age during final time quarter.

## Value

Matrix of age distribution over time.
