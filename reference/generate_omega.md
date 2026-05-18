# Generate migration patterns

Generate migration patterns

## Usage

``` r
generate_omega(nquar, pattern, target_mig = 12000, dist = TRUE)
```

## Arguments

- nquar:

  Number of time quarters.

- pattern:

  Pattern type for omega ("constant", "increasing", "step_change_1",
  "step_change_2", "none").

- target_mig:

  Target total number of migrations to simulate (may be higher/lower
  with Poisson variation).

- dist:

  Logical indicating whether to add Poisson variation.

## Value

Vector of omega values over time.
