# Helper function to bind tidy data frames

Helper function to bind tidy data frames

## Usage

``` r
tidy_bind(
  list,
  age_bands,
  q_levels,
  nquar,
  start_yr,
  qt.flag = TRUE,
  CD4 = FALSE
)
```

## Arguments

- list:

  List of HIV posterior predictive matrices.

- age_bands:

  Age bands list.

- q_levels:

  Quartile levels.

- nquar:

  Number of quarters.

- start_yr:

  Starting year.

- qt.flag:

  Logical indicating if quarterly output is required.

- CD4:

  Logical indicating if CD4 data is being processed.

## Value

Combined tidy data frame.
